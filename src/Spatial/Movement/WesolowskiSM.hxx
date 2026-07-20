/*
 * WesolowskiSM.hxx
 *
 * Movement model based upon gravity model in
 * https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1004267
 */
#ifndef SPATIAL_WESOLOWSKISM_H
#define SPATIAL_WESOLOWSKISM_H

#include <cmath>

#include "Spatial/GIS/GridPairTable.h"
#include "Spatial/SpatialModel.hxx"
#include "Utils/Helpers/NumberHelpers.h"
#include "Utils/TypeDef.h"

namespace Spatial {
class WesolowskiSM : public SpatialModel {
public:
  // Disallow copy
  WesolowskiSM(const WesolowskiSM &) = delete;
  WesolowskiSM &operator=(const WesolowskiSM &) = delete;

  // Disallow move
  WesolowskiSM(WesolowskiSM &&) = delete;
  WesolowskiSM &operator=(WesolowskiSM &&) = delete;

  double kappa_;
  double alpha_;
  double beta_;
  double gamma_;

  [[nodiscard]] double get_kappa() const { return kappa_; }
  void set_kappa(const double &value) { kappa_ = value; }

  [[nodiscard]] double get_alpha() const { return alpha_; }
  void set_alpha(const double &value) { alpha_ = value; }

  [[nodiscard]] double get_beta() const { return beta_; }
  void set_beta(const double &value) { beta_ = value; }

  [[nodiscard]] double get_gamma() const { return gamma_; }
  void set_gamma(const double &value) {
    gamma_ = value;
    distance_power_ = GridPairTable{};
  }

  explicit WesolowskiSM(double kappa, double alpha, double beta, double gamma)
      : kappa_(kappa), alpha_(alpha), beta_(beta), gamma_(gamma) {}

  ~WesolowskiSM() override = default;

  void prepare() override;

  // Public API intentionally remains identical to 500054a in both build modes.
  [[nodiscard]] DoubleVector get_v_relative_out_movement_to_destination(
      const int &from_location,
      const int &number_of_locations,
      const DoubleVector &relative_distance_vector,
      const IntVector &v_number_of_residents_by_location) const override;

private:
  GridPairTable distance_power_;
};
}  // namespace Spatial

#endif
