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
  void set_tau(const double& value) { tau_ = value; }

  [[nodiscard]] double get_alpha() const { return alpha_; }
  void set_alpha(const double& value) {
    alpha_ = value;
    kernel_ = LocationPairTable{};
  }

  [[nodiscard]] double get_rho() const { return log_rho_; }
  void set_log_rho(const double& value) {
    log_rho_ = value;
    kernel_ = LocationPairTable{};
  }

  double tau_;
  double alpha_;
  double log_rho_;
  int number_of_locations_;

  // Borrowed, owned by SpatialSettings and outlives this object.
  const LocationPairTable* spatial_distance_{nullptr};

  // The kernel is a pure function of distance, so it shares the compact
  // representation of the distance table instead of being a second n*n array.
  LocationPairTable kernel_;

  void prepare_kernel() {
    if (spatial_distance_ == nullptr) {
      throw std::runtime_error(fmt::format("{} called without spatial distances", __FUNCTION__));
    }

    const double log_rho = log_rho_;
    const double alpha = alpha_;
    // Preserves the previous behavior that excluded zero-distance pairs,
    // including distinct colocated locations in the dense fallback, and verifies
    // that no real pair maps onto the sentinel.
    kernel_ = spatial_distance_->map_with_zero_sentinel(
        [log_rho, alpha](double distance) { return std::pow(1.0 + (distance / log_rho), -alpha); },
        "MarshallSM kernel");
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
      const int& from_location, const int& number_of_locations,
      const IntVector& v_number_of_residents_by_location) const override {
    if (kernel_.empty()) {
      throw std::runtime_error(fmt::format("{} called without kernel prepared", __FUNCTION__));
    }

    const double source_population_power =
        std::pow(v_number_of_residents_by_location[from_location], tau_);
    const auto kernel_row = kernel_.row_view(from_location);

    DoubleVector results(number_of_locations, 0.0);
    for (int destination = 0; destination < number_of_locations; ++destination) {
      const double kernel_value = kernel_row[static_cast<size_t>(destination)];
      if (NumberHelpers::is_zero(kernel_value)) { continue; }
      results[destination] = source_population_power * kernel_value;
    }
    return results;
  }
};
}  // namespace Spatial

#endif
