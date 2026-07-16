#include "BurkinaFasoSM.h"

#include "Simulation/Model.h"

void Spatial::BurkinaFasoSM::prepare() {
  prepare_kernel();
  spdlog::info("Kernel prepared for BurkinaFasoSM, {} locations, {:.1f} MB", kernel_.size(),
               kernel_.memory_bytes() / 1048576.0);

  travel_.clear();
  if (Model::get_spatial_data() != nullptr) {
    AscFile* travel_raster =
        Model::get_spatial_data()->get_raster(SpatialData::SpatialFileType::TRAVEL);
    if (travel_raster != nullptr) {
      travel_ = std::move(prepare_surface(travel_raster, static_cast<int>(number_of_locations_)));
      spdlog::info("Travel raster prepared for BurkinaFasoSM, travel size: {}", travel_.size());
    } else {
      spdlog::warn(
          "Travel raster is not set for BurkinaFasoSM, no travel friction will be applied.");
    }
  } else {
    spdlog::warn("BurkinaFasoSM: no spatial data found, surface travel not prepared.");
  }

  prepare_districts();
}

void Spatial::BurkinaFasoSM::prepare_kernel() {
  if (spatial_distance_ == nullptr) {
    throw std::runtime_error(fmt::format("{} called without spatial distances", __FUNCTION__));
  }

  spdlog::info("Preparing kernel for BurkinaFasoSM, number of locations: {}",
               number_of_locations_);
  const double rho = rho_;
  const double alpha = alpha_;
  kernel_ = spatial_distance_->map_with_zero_sentinel(
      [rho, alpha](double distance) { return std::pow(1.0 + (distance / rho), -alpha); },
      "BurkinaFasoSM kernel");
}

void Spatial::BurkinaFasoSM::prepare_districts() {
  // The destination loop used to call get_admin_unit("district", destination),
  // which hashes a std::string and looks it up in a map, once per destination.
  // The mapping is fixed for the whole run, so resolve it here instead.
  has_district_level_ = false;
  district_by_location_.clear();

  if (Model::get_spatial_data() == nullptr
      || !Model::get_spatial_data()->has_admin_level("district")) {
    spdlog::info("BurkinaFasoSM: no district admin level, capital penalty will not be applied.");
    return;
  }

  district_by_location_.resize(number_of_locations_);
  for (uint64_t location = 0; location < number_of_locations_; location++) {
    district_by_location_[location] =
        Model::get_spatial_data()->get_admin_unit("district", static_cast<int>(location));
  }
  has_district_level_ = true;
  spdlog::info("BurkinaFasoSM: cached district for {} locations.", number_of_locations_);
}

DoubleVector Spatial::BurkinaFasoSM::get_v_relative_out_movement_to_destination(
    const int& from_location, const int& number_of_locations,
    const IntVector& v_number_of_residents_by_location) const {
  if (kernel_.empty()) {
    throw std::runtime_error(fmt::format("{} called without kernel prepared", __FUNCTION__));
  }

  const double source_population_power =
      std::pow(v_number_of_residents_by_location[from_location], tau_);
  const auto kernel_row = kernel_.row_view(from_location);
  const bool use_travel = travel_.size() == static_cast<size_t>(number_of_locations);
  const double source_travel = use_travel ? travel_[from_location] : 0.0;

  // If the source is not in the capital the penalty can never fire, so the
  // per-destination district read is not needed at all.
  const bool source_in_capital =
      has_district_level_ && (district_by_location_[from_location] == capital_);

  DoubleVector results(number_of_locations, 0.0);
  for (int destination = 0; destination < number_of_locations; ++destination) {
    const double kernel_value = kernel_row[static_cast<size_t>(destination)];
    if (NumberHelpers::is_zero(kernel_value)) { continue; }

    double probability = source_population_power * kernel_value;

    if (use_travel) {
      probability /= 1.0 + source_travel + travel_[destination];
    }

    if (source_in_capital && district_by_location_[destination] == capital_) {
      probability /= penalty_;
    }

    results[destination] = probability;
  }

  return results;
}
