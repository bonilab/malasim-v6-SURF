/*
 * MarshallSM.hxx
 *
 * Gravity model for human migration based upon a distance kernel function.
 *
 * Marshall et al., 2018
 */
#ifndef MARSHALLSM_HXX
#define MARSHALLSM_HXX

#include "Spatial/GIS/LocationPairTable.h"
#include "Spatial/SpatialModel.hxx"
#include "Utils/Helpers/NumberHelpers.h"
#include "Utils/TypeDef.h"

namespace Spatial {
class MarshallSM : public SpatialModel {
public:
  // Disallow copy
  MarshallSM(const MarshallSM&) = delete;
  MarshallSM& operator=(const MarshallSM&) = delete;

  // Disallow move
  MarshallSM(MarshallSM&&) = delete;
  MarshallSM& operator=(MarshallSM&&) = delete;

  [[nodiscard]] double get_tau() const { return tau_; }
  void set_tau(const double &value) { tau_ = value; }

  [[nodiscard]] double get_alpha() const { return alpha_; }
  void set_alpha(const double &value) { alpha_ = value; }

  [[nodiscard]] double get_rho() const { return log_rho_; }
  void set_log_rho(const double &value) { log_rho_ = value; }

  double tau_;
  double alpha_;
  double log_rho_;
  int number_of_locations_;

  // Borrowed, owned by SpatialSettings and outlives this object. Previously this
  // was a by-value copy of the whole n*n matrix.
  const LocationPairTable* spatial_distance_{nullptr};

  // The kernel is a pure function of distance, so it shares the compact
  // representation of the distance table instead of being a second n*n array.
  LocationPairTable kernel_;

  // Precompute the kernel function for the movement model
  void prepare_kernel() {
    const double log_rho = log_rho_;
    const double alpha = alpha_;
    kernel_ = spatial_distance_->map(
        [log_rho, alpha](double distance) { return std::pow(1 + (distance / log_rho), (-alpha)); });
  }

  explicit MarshallSM(double tau, double alpha, double log_rho, int number_of_locations,
                      const LocationPairTable* spatial_distance)
      : tau_(tau),
        alpha_(alpha),
        log_rho_(log_rho),
        number_of_locations_(number_of_locations),
        spatial_distance_(spatial_distance) {}

  ~MarshallSM() override = default;

  void prepare() override { prepare_kernel(); }

  [[nodiscard]] DoubleVector get_v_relative_out_movement_to_destination(
      const int &from_location, const int &number_of_locations,
      const DoubleVector &relative_distance_vector,
      const IntVector &v_number_of_residents_by_location) const override {
    // Note the population size
    auto population = v_number_of_residents_by_location[from_location];

    // Prepare the vector for results
    std::vector<double> results(number_of_locations, 0.0);

    for (auto destination = 0; destination < number_of_locations;
         destination++) {
      // Continue if there is nothing to do
      if (NumberHelpers::is_zero(relative_distance_vector[destination])) {
        continue;
      }

      // Calculate the proportional probability
      double probability = std::pow(population, tau_) * kernel_.at(from_location, destination);
      results[destination] = probability;
    }

    // Done, return the results
    return results;
  }
};
}  // namespace Spatial

#endif
