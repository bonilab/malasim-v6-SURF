#include "Spatial/Movement/WesolowskiSurfaceSM.h"

#include "Simulation/Model.h"

void Spatial::WesolowskiSurfaceSM::prepare() {
#ifdef USE_DISTANCE_LUT
  distance_power_ = LocationPairTable{};
  if (Model::get_config() != nullptr) {
    const auto &distances =
        Model::get_config()->get_spatial_settings().get_spatial_distance_lut();
    if (!distances.empty()) {
      const double gamma = gamma_;
      distance_power_ = distances.map_with_zero_sentinel(
          [gamma](double distance) { return std::pow(distance, gamma); },
          "WesolowskiSurfaceSM distance powers");
    }
  }
#endif

  AscFile* travel_raster =
      Model::get_spatial_data()->get_raster(SpatialData::SpatialFileType::TRAVEL);
  travel = std::move(prepare_surface(travel_raster, number_of_locations_));
}

DoubleVector Spatial::WesolowskiSurfaceSM::get_v_relative_out_movement_to_destination(
    const int &from_location, const int &number_of_locations,
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

#ifdef USE_DISTANCE_LUT
  if (distance_power_.size() == static_cast<size_t>(number_of_locations)) {
    const auto distance_power = distance_power_.row_view(from_location);
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
#endif

  // Compatibility fallback for standalone/unit-test construction and old mode.
#ifdef USE_DISTANCE_LUT
  if (relative_distance_vector.size() < static_cast<size_t>(number_of_locations)) {
    throw std::runtime_error(fmt::format(
        "WesolowskiSurfaceSM called without a prepared distance LUT or compatibility distance row"));
  }
#endif
  for (int destination = 0; destination < number_of_locations; ++destination) {
    const double distance = relative_distance_vector[destination];
    if (NumberHelpers::is_zero(distance)) { continue; }

    const double probability =
        kappa_
        * (source_population_power
           * std::pow(v_number_of_residents_by_location[destination], beta_))
        / std::pow(distance, gamma_);
    results[destination] = probability / (1.0 + source_travel + travel[destination]);
  }
  return results;
}
