#include "Spatial/Movement/WesolowskiSM.hxx"

#include "Simulation/Model.h"
#include "Spatial/Movement/MovementKernel.h"

void Spatial::WesolowskiSM::prepare() {
  distance_power_ = GridPairTable{};
  if (Model::get_config() == nullptr) { return; }

  const auto* provider = Model::get_config()->get_spatial_settings().get_distance_provider();
  const auto* distances = provider == nullptr ? nullptr : provider->grid_table();
  if (distances == nullptr) { return; }

  const double gamma = gamma_;
  distance_power_ = make_grid_movement_kernel(
      *distances, [gamma](double distance) { return std::pow(distance, gamma); },
      "WesolowskiSM distance powers");
}

DoubleVector Spatial::WesolowskiSM::get_v_relative_out_movement_to_destination(
    const int &from_location,
    const int &number_of_locations,
    const DoubleVector &relative_distance_vector,
    const IntVector &v_number_of_residents_by_location) const {
  DoubleVector results(number_of_locations, 0.0);
  const double source_population_power =
      std::pow(v_number_of_residents_by_location[from_location], alpha_);

  if (distance_power_.size() == static_cast<size_t>(number_of_locations)) {
    const auto distance_power = distance_power_.values_from(from_location);
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

  // Compatibility fallback for standalone/unit-test construction and old mode.
  if (relative_distance_vector.size() < static_cast<size_t>(number_of_locations)) {
    throw std::runtime_error(fmt::format(
        "WesolowskiSM called without a prepared distance LUT or compatibility distance row"));
  }
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
