#include "Spatial/Movement/WesolowskiSM.hxx"

#ifdef USE_DISTANCE_LUT
#include "Simulation/Model.h"
#endif

void Spatial::WesolowskiSM::prepare() {
#ifdef USE_DISTANCE_LUT
  distance_power_ = LocationPairTable{};
  if (Model::get_config() == nullptr) { return; }

  const auto &distances = Model::get_config()->get_spatial_settings().get_spatial_distance_lut();
  if (distances.empty()) { return; }

  const double gamma = gamma_;
  distance_power_ = distances.map_with_zero_sentinel(
      [gamma](double distance) { return std::pow(distance, gamma); },
      "WesolowskiSM distance powers");
#endif
}

DoubleVector Spatial::WesolowskiSM::get_v_relative_out_movement_to_destination(
    const int &from_location, const int &number_of_locations,
    const DoubleVector &relative_distance_vector,
    const IntVector &v_number_of_residents_by_location) const {
  DoubleVector results(number_of_locations, 0.0);
  const double source_population_power =
      std::pow(v_number_of_residents_by_location[from_location], alpha_);

#ifdef USE_DISTANCE_LUT
  if (distance_power_.size() == static_cast<size_t>(number_of_locations)) {
    const auto distance_power = distance_power_.row_view(from_location);
    for (int target_location = 0; target_location < number_of_locations; ++target_location) {
      const double denominator = distance_power[static_cast<size_t>(target_location)];
      if (NumberHelpers::is_zero(denominator)) { continue; }
      results[target_location] =
          kappa_
          * (source_population_power
             * std::pow(v_number_of_residents_by_location[target_location], beta_))
          / denominator;
    }
    return results;
  }
#endif

  // Compatibility fallback for standalone/unit-test construction and old mode.
#ifdef USE_DISTANCE_LUT
  if (relative_distance_vector.size() < static_cast<size_t>(number_of_locations)) {
    throw std::runtime_error(fmt::format(
        "WesolowskiSM called without a prepared distance LUT or compatibility distance row"));
  }
#endif
  for (int target_location = 0; target_location < number_of_locations; ++target_location) {
    const double distance = relative_distance_vector[target_location];
    if (NumberHelpers::is_zero(distance)) { continue; }
    results[target_location] =
        kappa_
        * (source_population_power
           * std::pow(v_number_of_residents_by_location[target_location], beta_))
        / std::pow(distance, gamma_);
  }
  return results;
}
