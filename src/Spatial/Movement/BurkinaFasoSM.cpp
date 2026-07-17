#include "Spatial/Movement/BurkinaFasoSM.h"

#include "Simulation/Model.h"

Spatial::BurkinaFasoSM::BurkinaFasoSM(
    double tau, double alpha, double rho, double capital, double penalty,
    int number_of_locations, std::vector<std::vector<double>> spatial_distance_matrix)
    : tau_(tau),
      alpha_(alpha),
      rho_(rho),
      capital_(capital),
      penalty_(penalty),
      number_of_locations_(number_of_locations),
      spatial_distance_matrix_(std::move(spatial_distance_matrix)) {
#ifdef USE_DISTANCE_LUT
  if (!spatial_distance_matrix_.empty()) {
    constructor_distance_lut_ = LocationPairTable::make_dense(spatial_distance_matrix_);
  }
#endif
}

void Spatial::BurkinaFasoSM::prepare() {
  prepare_kernel();
#ifdef USE_DISTANCE_LUT
  spdlog::info("Kernel prepared for BurkinaFasoSM, {} locations, {:.1f} MB",
               kernel_lut_.size(), kernel_lut_.memory_bytes() / 1048576.0);
#else
  spdlog::info("Kernel prepared for BurkinaFasoSM, kernel size x,y: {} - {}", kernel_.size(),
               kernel_.empty() ? 0 : kernel_[0].size());
#endif

  travel_.clear();
  if (Model::get_spatial_data() != nullptr) {
    AscFile* travel_raster =
        Model::get_spatial_data()->get_raster(SpatialData::SpatialFileType::TRAVEL);
    if (travel_raster == nullptr) {
      spdlog::warn("BurkinaFasoSM: travel raster not found, surface travel not prepared.");
    } else {
      travel_ = std::move(prepare_surface(travel_raster, number_of_locations_));
      spdlog::info("BurkinaFasoSM: surface travel prepared, size: {}", travel_.size());
    }
  } else {
    spdlog::warn("BurkinaFasoSM: no spatial data found, surface travel not prepared.");
  }

#ifdef USE_DISTANCE_LUT
  prepare_districts();
#endif
}

void Spatial::BurkinaFasoSM::prepare_kernel() {
#ifdef USE_DISTANCE_LUT
  kernel_lut_ = LocationPairTable{};
  const LocationPairTable* distances = nullptr;

  if (!constructor_distance_lut_.empty()) {
    distances = &constructor_distance_lut_;
  } else if (Model::get_config() != nullptr) {
    const auto &configured =
        Model::get_config()->get_spatial_settings().get_spatial_distance_lut();
    if (!configured.empty()) { distances = &configured; }
  }

  if (distances == nullptr) { return; }

  const double rho = rho_;
  const double alpha = alpha_;
  kernel_lut_ = distances->map_with_zero_sentinel(
      [rho, alpha](double distance) {
        return std::pow(1.0 + (distance / rho), -alpha);
      },
      "BurkinaFasoSM kernel");
#else
  kernel_.assign(number_of_locations_, std::vector<double>(number_of_locations_, 0.0));
  for (uint64_t source = 0; source < number_of_locations_; ++source) {
    for (uint64_t destination = 0; destination < number_of_locations_; ++destination) {
      kernel_[source][destination] = std::pow(
          1.0 + (spatial_distance_matrix_[source][destination] / rho_), -alpha_);
    }
  }
#endif
}

#ifdef USE_DISTANCE_LUT
void Spatial::BurkinaFasoSM::prepare_districts() {
  has_district_level_ = false;
  district_by_location_.clear();

  if (Model::get_spatial_data() == nullptr
      || !Model::get_spatial_data()->has_admin_level("district")) {
    return;
  }

  district_by_location_.resize(number_of_locations_);
  for (uint64_t location = 0; location < number_of_locations_; ++location) {
    district_by_location_[location] =
        Model::get_spatial_data()->get_admin_unit("district", static_cast<int>(location));
  }
  has_district_level_ = true;
}
#endif

DoubleVector Spatial::BurkinaFasoSM::get_v_relative_out_movement_to_destination(
    const int &from_location, const int &number_of_locations,
    const DoubleVector &relative_distance_vector,
    const IntVector &v_number_of_residents_by_location) const {
  const double source_population_power =
      std::pow(v_number_of_residents_by_location[from_location], tau_);
  const bool use_travel = travel_.size() == static_cast<size_t>(number_of_locations);
  const double source_travel = use_travel ? travel_[from_location] : 0.0;
  DoubleVector results(number_of_locations, 0.0);

#ifdef USE_DISTANCE_LUT
  const bool use_kernel = kernel_lut_.size() == static_cast<size_t>(number_of_locations);
  if (!use_kernel
      && relative_distance_vector.size() < static_cast<size_t>(number_of_locations)) {
    throw std::runtime_error(fmt::format(
        "BurkinaFasoSM called without a prepared distance LUT or compatibility distance row"));
  }
  const auto kernel_row = use_kernel ? kernel_lut_.row_view(from_location)
                                     : LocationPairTable::RowView{};
  const bool source_in_capital =
      has_district_level_ && district_by_location_[from_location] == capital_;
#else
  if (kernel_.empty()) {
    throw std::runtime_error(fmt::format("{} called without kernel prepared", __FUNCTION__));
  }
#endif

  for (int destination = 0; destination < number_of_locations; ++destination) {
    double kernel_value = 0.0;
#ifdef USE_DISTANCE_LUT
    if (use_kernel) {
      kernel_value = kernel_row[static_cast<size_t>(destination)];
    } else {
      const double distance = relative_distance_vector[destination];
      if (NumberHelpers::is_zero(distance)) { continue; }
      kernel_value = std::pow(1.0 + (distance / rho_), -alpha_);
    }
#else
    if (NumberHelpers::is_zero(relative_distance_vector[destination])) { continue; }
    kernel_value = kernel_[from_location][destination];
#endif
    if (NumberHelpers::is_zero(kernel_value)) { continue; }

    double probability = source_population_power * kernel_value;
    if (use_travel) {
      probability /= 1.0 + source_travel + travel_[destination];
    }

#ifdef USE_DISTANCE_LUT
    if (source_in_capital && district_by_location_[destination] == capital_) {
      probability /= penalty_;
    }
#else
    if (Model::get_spatial_data() != nullptr
        && Model::get_spatial_data()->has_admin_level("district")) {
      const auto source_district =
          Model::get_spatial_data()->get_admin_unit("district", from_location);
      if (source_district == capital_
          && Model::get_spatial_data()->get_admin_unit("district", destination) == capital_) {
        probability /= penalty_;
      }
    }
#endif

    results[destination] = probability;
  }

  return results;
}
