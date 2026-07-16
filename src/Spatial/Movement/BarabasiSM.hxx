/*
 * BarabasiSM.hxx
 *
 * Movement model based upon the radius of gyration distribution in
 * https://www.nature.com/articles/nature06958
 *
 * REMINDER Verify the correctness of the equation (2023-05-05)
 */
#ifndef SPATIAL_BARABASISM_HXX
#define SPATIAL_BARABASISM_HXX

#include <cmath>

#include "Spatial/GIS/LocationPairTable.h"
#include "Spatial/SpatialModel.hxx"
#include "Utils/Helpers/NumberHelpers.h"

namespace Spatial {
class BarabasiSM : public SpatialModel {
public:
  // Disallow copy
  BarabasiSM(const BarabasiSM&) = delete;
  BarabasiSM& operator=(const BarabasiSM&) = delete;

  // Disallow move
  BarabasiSM(BarabasiSM&&) = delete;
  BarabasiSM& operator=(BarabasiSM&&) = delete;

  [[nodiscard]] double get_r_g_0() const { return r_g_0_; }
  void set_r_g_0(const double& value) {
    r_g_0_ = value;
    movement_weight_ = LocationPairTable{};
  }

  [[nodiscard]] double get_beta_r() const { return beta_r_; }
  void set_beta_r(const double& value) {
    beta_r_ = value;
    movement_weight_ = LocationPairTable{};
  }

  [[nodiscard]] double get_kappa() const { return kappa_; }
  void set_kappa(const double& value) {
    kappa_ = value;
    movement_weight_ = LocationPairTable{};
  }

  explicit BarabasiSM(double r_g_0, double beta_r, double kappa,
                      const LocationPairTable* spatial_distance)
      : r_g_0_(r_g_0),
        beta_r_(beta_r),
        kappa_(kappa),
        spatial_distance_(spatial_distance) {}

  ~BarabasiSM() override = default;

  void prepare() override {
    if (spatial_distance_ == nullptr) {
      throw std::runtime_error(fmt::format("{} called without spatial distances", __FUNCTION__));
    }

    const double r_g_0 = r_g_0_;
    const double beta_r = beta_r_;
    const double kappa = kappa_;
    movement_weight_ = spatial_distance_->map_with_zero_sentinel(
        [r_g_0, beta_r, kappa](double distance) {
          return std::pow(distance + r_g_0, -beta_r) * std::exp(-distance / kappa);
        },
        "BarabasiSM movement weights");
  }

  [[nodiscard]] DoubleVector get_v_relative_out_movement_to_destination(
      const int& from_location, const int& number_of_locations,
      const IntVector& v_number_of_residents_by_location) const override {
    (void)v_number_of_residents_by_location;
    if (movement_weight_.empty()) {
      throw std::runtime_error(fmt::format("{} called without weights prepared", __FUNCTION__));
    }

    DoubleVector results(number_of_locations, 0.0);
    const auto weights = movement_weight_.row_view(from_location);
    for (int target_location = 0; target_location < number_of_locations; ++target_location) {
      results[target_location] = weights[static_cast<size_t>(target_location)];
    }
    return results;
  }

private:
  double r_g_0_;
  double beta_r_;
  double kappa_;

  const LocationPairTable* spatial_distance_{nullptr};
  LocationPairTable movement_weight_;
};
}  // namespace Spatial

#endif
