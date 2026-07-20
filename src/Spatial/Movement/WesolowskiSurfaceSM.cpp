#include "Spatial/Movement/WesolowskiSurfaceSM.h"

#include "Simulation/Model.h"
#include "Spatial/Movement/MovementKernel.h"

void Spatial::WesolowskiSurfaceSM::prepare() {
  distance_power_ = GridPairTable{};
  if (Model::get_config() != nullptr) {
    const auto* provider = Model::get_config()->get_spatial_settings().get_distance_provider();
    const auto* distances = provider == nullptr ? nullptr : provider->grid_table();
    if (distances != nullptr) {
      const double gamma = gamma_;
      distance_power_ = make_grid_movement_kernel(
          *distances, [gamma](double distance) { return std::pow(distance, gamma); },
          "WesolowskiSurfaceSM distance powers");
    }
  }

  AscFile* travel_raster =
      Model::get_spatial_data()->get_raster(SpatialData::SpatialFileType::TRAVEL);
  travel = std::move(prepare_surface(travel_raster, number_of_locations_));
}

DoubleVector Spatial::WesolowskiSurfaceSM::get_v_relative_out_movement_to_destination(
    const int &from_location,
    const int &number_of_locations,
    const DoubleVector &relative_distance_vector,
    const IntVector &v_number_of_residents_by_location) const {
  if (travel.empty()) {
    throw std::runtime_error(
        fmt::format("{} called without travel surface prepared", __FUNCTION__));
  }

  DoubleVector results(number_of_locations, 0.0);
  const double source_population_power =
      std::pow(v_number_of_residents_by_location[from_location], alpha_);
  const double source_travel = travel[from_location];

  if (distance_power_.size() == static_cast<size_t>(number_of_locations)) {
    const auto distance_power = distance_power_.values_from(from_location);
    for (int destination = 0; destination < number_of_locations; ++destination) {
      const double denominator = distance_power[static_cast<size_t>(destination)];
      if (NumberHelpers::is_zero(denominator)) { continue; }

      const double probability =
          kappa_
          * (source_population_power
             * std::pow(v_number_of_residents_by_location[destination], beta_))
          / denominator;
      results[destination] = probability / (1.0 + source_travel + travel[destination]);
    }
    return results;
  }

  // Compatibility fallback for standalone/unit-test construction and old mode.
  if (relative_distance_vector.size() < static_cast<size_t>(number_of_locations)) {
    throw std::runtime_error(
        fmt::format("WesolowskiSurfaceSM called without a prepared distance LUT or compatibility "
                    "distance row"));
  }
  for (int destination = 0; destination < number_of_locations; ++destination) {
    const double distance = relative_distance_vector[destination];
    if (NumberHelpers::is_zero(distance)) { continue; }

    const double probability = kappa_
                               * (source_population_power
                                  * std::pow(v_number_of_residents_by_location[destination], beta_))
                               / std::pow(distance, gamma_);
    results[destination] = probability / (1.0 + source_travel + travel[destination]);
  }
  return results;
}
