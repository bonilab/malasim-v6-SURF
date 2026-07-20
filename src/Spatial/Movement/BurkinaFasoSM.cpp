#include "Spatial/Movement/BurkinaFasoSM.h"

#include "Simulation/Model.h"
#include "Spatial/Movement/MovementKernel.h"

Spatial::BurkinaFasoSM::BurkinaFasoSM(
    double tau, double alpha, double rho, double capital, double penalty, int number_of_locations)
    : tau_(tau),
      alpha_(alpha),
      rho_(rho),
      capital_(capital),
      penalty_(penalty),
      number_of_locations_(number_of_locations) {}

void Spatial::BurkinaFasoSM::prepare() {
  prepare_kernel();
  spdlog::info("Kernel prepared for BurkinaFasoSM, {} locations, {:.1f} MB", kernel_lut_.size(),
               kernel_lut_.memory_bytes() / 1048576.0);

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

  prepare_districts();
}

void Spatial::BurkinaFasoSM::prepare_kernel() {
  kernel_lut_ = GridPairTable{};
  const GridPairTable* distances = nullptr;
  if (Model::get_config() != nullptr) {
    const auto* provider = Model::get_config()->get_spatial_settings().get_distance_provider();
    distances = provider == nullptr ? nullptr : provider->grid_table();
  }

  if (distances == nullptr) { return; }

  const double rho = rho_;
  const double alpha = alpha_;
  kernel_lut_ = make_grid_movement_kernel(
      *distances,
      [rho, alpha](double distance) { return std::pow(1.0 + (distance / rho), -alpha); },
      "BurkinaFasoSM kernel");
}

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

DoubleVector Spatial::BurkinaFasoSM::get_v_relative_out_movement_to_destination(
    const int &from_location,
    const int &number_of_locations,
    const DoubleVector &relative_distance_vector,
    const IntVector &v_number_of_residents_by_location) const {
  const double source_population_power =
      std::pow(v_number_of_residents_by_location[from_location], tau_);
  const bool use_travel = travel_.size() == static_cast<size_t>(number_of_locations);
  const double source_travel = use_travel ? travel_[from_location] : 0.0;
  DoubleVector results(number_of_locations, 0.0);

  const bool use_kernel = kernel_lut_.size() == static_cast<size_t>(number_of_locations);
  if (!use_kernel && relative_distance_vector.size() < static_cast<size_t>(number_of_locations)) {
    throw std::runtime_error(fmt::format(
        "BurkinaFasoSM called without a prepared distance LUT or compatibility distance row"));
  }
  const auto kernel_values =
      use_kernel ? kernel_lut_.values_from(from_location) : GridPairTable::ValuesFromLocation{};
  const bool source_in_capital =
      has_district_level_ && district_by_location_[from_location] == capital_;

  for (int destination = 0; destination < number_of_locations; ++destination) {
    double kernel_value = 0.0;
    if (use_kernel) {
      kernel_value = kernel_values[static_cast<size_t>(destination)];
    } else {
      const double distance = relative_distance_vector[destination];
      if (NumberHelpers::is_zero(distance)) { continue; }
      kernel_value = std::pow(1.0 + (distance / rho_), -alpha_);
    }
    if (NumberHelpers::is_zero(kernel_value)) { continue; }

    double probability = source_population_power * kernel_value;
    if (use_travel) { probability /= 1.0 + source_travel + travel_[destination]; }

    if (source_in_capital && district_by_location_[destination] == capital_) {
      probability /= penalty_;
    }

    results[destination] = probability;
  }

  return results;
}
