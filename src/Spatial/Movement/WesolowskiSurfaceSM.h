/*
 * WesolowskiSurfaceSM.h
 *
 * Movement model based upon gravity model in
 * https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1004267
 * with the travel surface penalty from
 * https://www.nature.com/articles/s41598-022-26878-5 applied to the results.
 */
#ifndef SPATIAL_WESOLOWSKISURFACESM_H
#define SPATIAL_WESOLOWSKISURFACESM_H

#include <cmath>

#include "Spatial/GIS/LocationPairTable.h"
#include "Spatial/SpatialModel.hxx"
#include "Utils/Helpers/NumberHelpers.h"
#include "Utils/TypeDef.h"

namespace Spatial {
class WesolowskiSurfaceSM : public SpatialModel {
public:
  // Disallow copy
  WesolowskiSurfaceSM(const WesolowskiSurfaceSM&) = delete;
  WesolowskiSurfaceSM& operator=(const WesolowskiSurfaceSM&) = delete;

  // Disallow move
  WesolowskiSurfaceSM(WesolowskiSurfaceSM&&) = delete;
  WesolowskiSurfaceSM& operator=(WesolowskiSurfaceSM&&) = delete;

  double kappa_;
  double alpha_;
  double beta_;
  double gamma_;
  int number_of_locations_;

  // Travel surface, computed when the prepare method is called.
  std::vector<double> travel;

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

  explicit WesolowskiSurfaceSM(double kappa, double alpha, double beta, double gamma,
                               int number_of_locations,
                               const LocationPairTable* spatial_distance)
      : kappa_(kappa),
        alpha_(alpha),
        beta_(beta),
        gamma_(gamma),
        number_of_locations_(number_of_locations),
        spatial_distance_(spatial_distance) {}

  ~WesolowskiSurfaceSM() override = default;

  void prepare() override;

  [[nodiscard]] DoubleVector get_v_relative_out_movement_to_destination(
      const int& from_location, const int& number_of_locations,
      const IntVector& v_number_of_residents_by_location) const override {
    if (travel.empty()) {
      throw std::runtime_error(
          fmt::format("{} called without travel surface prepared", __FUNCTION__));
    }
    if (distance_power_.empty()) {
      throw std::runtime_error(
          fmt::format("{} called without distance powers prepared", __FUNCTION__));
    }

    DoubleVector results(number_of_locations, 0.0);
    const auto distance_power = distance_power_.row_view(from_location);
    const double source_population_power =
        std::pow(v_number_of_residents_by_location[from_location], alpha_);
    const double source_travel = travel[from_location];

    for (int destination = 0; destination < number_of_locations; ++destination) {
      const double denominator = distance_power[static_cast<size_t>(destination)];
      if (NumberHelpers::is_zero(denominator)) { continue; }

      // Gravity model:
      // N_{ij}=kappa * pop_i^alpha * pop_j^beta / d(i,j)^gamma
      const double probability =
          kappa_
          * (source_population_power
             * std::pow(v_number_of_residents_by_location[destination], beta_))
          / denominator;

      // Travel penalty: Pr(j|i)' = Pr(j|i) / (1 + t_i + t_j)
      results[destination] = probability / (1.0 + source_travel + travel[destination]);
    }
    return results;
  }

private:
  const LocationPairTable* spatial_distance_{nullptr};
  LocationPairTable distance_power_;
};
}  // namespace Spatial

#endif
