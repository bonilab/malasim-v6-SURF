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

#ifdef USE_DISTANCE_LUT
#include "Spatial/GIS/LocationPairTable.h"
#endif
#include "Spatial/SpatialModel.hxx"
#include "Utils/Helpers/NumberHelpers.h"

namespace Spatial {
class BarabasiSM : public SpatialModel {
public:
  // Disallow copy
  BarabasiSM(const BarabasiSM &) = delete;
  BarabasiSM &operator=(const BarabasiSM &) = delete;

  // Disallow move
  BarabasiSM(BarabasiSM &&) = delete;
  BarabasiSM &operator=(BarabasiSM &&) = delete;

  [[nodiscard]] double get_r_g_0() const { return r_g_0_; }
  void set_r_g_0(const double &value) {
    r_g_0_ = value;
#ifdef USE_DISTANCE_LUT
    movement_weight_ = LocationPairTable{};
#endif
  }

  [[nodiscard]] double get_beta_r() const { return beta_r_; }
  void set_beta_r(const double &value) {
    beta_r_ = value;
#ifdef USE_DISTANCE_LUT
    movement_weight_ = LocationPairTable{};
#endif
  }

  [[nodiscard]] double get_kappa() const { return kappa_; }
  void set_kappa(const double &value) {
    kappa_ = value;
#ifdef USE_DISTANCE_LUT
    movement_weight_ = LocationPairTable{};
#endif
  }

  explicit BarabasiSM(double r_g_0, double beta_r, double kappa)
      : r_g_0_(r_g_0), beta_r_(beta_r), kappa_(kappa) {}

  ~BarabasiSM() override = default;

  void prepare() override;

  // Public API intentionally remains identical to 500054a in both build modes.
  [[nodiscard]] std::vector<double> get_v_relative_out_movement_to_destination(
      const int &from_location, const int &number_of_locations,
      const std::vector<double> &relative_distance_vector,
      const std::vector<int> &v_number_of_residents_by_location) const override;

private:
  double r_g_0_;
  double beta_r_;
  double kappa_;

#ifdef USE_DISTANCE_LUT
  LocationPairTable movement_weight_;
#endif
};
}  // namespace Spatial

#endif
