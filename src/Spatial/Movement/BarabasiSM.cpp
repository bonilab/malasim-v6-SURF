#include "Spatial/Movement/BarabasiSM.hxx"

#include "Simulation/Model.h"
#include "Spatial/Movement/MovementKernel.h"

void Spatial::BarabasiSM::prepare() {
  movement_weight_ = GridPairTable{};
  if (Model::get_config() == nullptr) { return; }

  const auto* provider = Model::get_config()->get_spatial_settings().get_distance_provider();
  const auto* distances = provider == nullptr ? nullptr : provider->grid_table();
  if (distances == nullptr) { return; }

  const double r_g_0 = r_g_0_;
  const double beta_r = beta_r_;
  const double kappa = kappa_;
  movement_weight_ = make_grid_movement_kernel(
      *distances,
      [r_g_0, beta_r, kappa](double distance) {
        return std::pow(distance + r_g_0, -beta_r) * std::exp(-distance / kappa);
      },
      "BarabasiSM movement weights");
}

std::vector<double> Spatial::BarabasiSM::get_v_relative_out_movement_to_destination(
    const int &from_location,
    const int &number_of_locations,
    const std::vector<double> &relative_distance_vector,
    const std::vector<int> &v_number_of_residents_by_location) const {
  (void)v_number_of_residents_by_location;
  std::vector<double> results(number_of_locations, 0.0);

  if (movement_weight_.size() == static_cast<size_t>(number_of_locations)) {
    const auto weights = movement_weight_.values_from(from_location);
    for (int target_location = 0; target_location < number_of_locations; ++target_location) {
      results[target_location] = weights[static_cast<size_t>(target_location)];
    }
    return results;
  }

  // Compatibility fallback for standalone/unit-test construction and old mode.
  if (relative_distance_vector.size() < static_cast<size_t>(number_of_locations)) {
    throw std::runtime_error(fmt::format(
        "BarabasiSM called without a prepared distance LUT or compatibility distance row"));
  }
  for (int target_location = 0; target_location < number_of_locations; ++target_location) {
    const double distance = relative_distance_vector[target_location];
    if (NumberHelpers::is_zero(distance)) { continue; }
    results[target_location] = std::pow(distance + r_g_0_, -beta_r_) * std::exp(-distance / kappa_);
  }
  return results;
}
