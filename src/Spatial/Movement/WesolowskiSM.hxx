/*
 * WesolowskiSM.hxx
 *
 * Movement model based upon gravity model in
 * https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1004267
 */
#ifndef SPATIAL_WESOLOWSKISM_H
#define SPATIAL_WESOLOWSKISM_H

#include <cmath>

#include "Spatial/GIS/LocationPairTable.h"
#include "Spatial/SpatialModel.hxx"
#include "Utils/Helpers/NumberHelpers.h"
#include "Utils/TypeDef.h"

namespace Spatial {
class WesolowskiSM : public SpatialModel {
public:
  // Disallow copy
  WesolowskiSM(const WesolowskiSM&) = delete;
  WesolowskiSM& operator=(const WesolowskiSM&) = delete;

  // Disallow move
  WesolowskiSM(WesolowskiSM&&) = delete;
  WesolowskiSM& operator=(WesolowskiSM&&) = delete;

  double kappa_;
  double alpha_;
  double beta_;
  double gamma_;

  [[nodiscard]] double get_kappa() const { return kappa_; }
  void set_kappa(const double& value) { kappa_ = value; }

  [[nodiscard]] double get_alpha() const { return alpha_; }
  void set_alpha(const double& value) { alpha_ = value; }

  [[nodiscard]] double get_beta() const { return beta_; }
  void set_beta(const double& value) { beta_ = value; }

  [[nodiscard]] double get_gamma() const { return gamma_; }
  void set_gamma(const double& value) {
    gamma_ = value;
    distance_power_ = LocationPairTable{};
  }

  explicit WesolowskiSM(double kappa, double alpha, double beta, double gamma,
                        const LocationPairTable* spatial_distance)
      : kappa_(kappa),
        alpha_(alpha),
        beta_(beta),
        gamma_(gamma),
        spatial_distance_(spatial_distance) {}

  ~WesolowskiSM() override = default;

  void prepare() override {
    if (spatial_distance_ == nullptr) {
      throw std::runtime_error(fmt::format("{} called without spatial distances", __FUNCTION__));
    }

    const double gamma = gamma_;
    distance_power_ = spatial_distance_->map_with_zero_sentinel(
        [gamma](double distance) { return std::pow(distance, gamma); },
        "WesolowskiSM distance powers");
  }

  [[nodiscard]] DoubleVector get_v_relative_out_movement_to_destination(
      const int& from_location, const int& number_of_locations,
      const IntVector& v_number_of_residents_by_location) const override {
    if (distance_power_.empty()) {
      throw std::runtime_error(
          fmt::format("{} called without distance powers prepared", __FUNCTION__));
    }

    DoubleVector results(number_of_locations, 0.0);
    const auto distance_power = distance_power_.row_view(from_location);
    const double source_population_power =
        std::pow(v_number_of_residents_by_location[from_location], alpha_);

    for (int target_location = 0; target_location < number_of_locations; ++target_location) {
      const double denominator = distance_power[static_cast<size_t>(target_location)];
      if (NumberHelpers::is_zero(denominator)) { continue; }

      // Gravity model:
      // N_{ij}=kappa * pop_i^alpha * pop_j^beta / d(i,j)^gamma
      results[target_location] =
          kappa_
          * (source_population_power
             * std::pow(v_number_of_residents_by_location[target_location], beta_))
          / denominator;
    }
    return results;
  }

private:
  const LocationPairTable* spatial_distance_{nullptr};

  // d(i,j)^gamma is a pure function of distance. Precomputing it once removes
  // billions of repeated pow() calls while retaining the compact LUT layout.
  LocationPairTable distance_power_;
};
}  // namespace Spatial

#endif
