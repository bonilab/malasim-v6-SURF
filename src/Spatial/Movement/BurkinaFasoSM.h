/*
 * BurkinaFaso.hxx
 *
 * Tuned movement model for Burkina Faso based upon Marshall et al. (2018),
 * with a penalty applied based upon the travel time to the nearest city.
 * Intradistrict movement in capital is also penalized as well.
 *
 * Marshall et al., 2018
 */
#ifndef BURKINAFASOSM_HXX
#define BURKINAFASOSM_HXX

#include <utility>

#include "Spatial/GIS/LocationPairTable.h"
#include "Spatial/SpatialModel.hxx"
#include "Utils/Helpers/NumberHelpers.h"
#include "Utils/TypeDef.h"

namespace Spatial {
class BurkinaFasoSM : public SpatialModel {
public:
  BurkinaFasoSM(const BurkinaFasoSM&) = delete;
  BurkinaFasoSM operator=(const BurkinaFasoSM&) = delete;
  BurkinaFasoSM(BurkinaFasoSM&&) = delete;
  BurkinaFasoSM& operator=(BurkinaFasoSM&&) = delete;

  [[nodiscard]] double get_tau() const { return tau_; }
  void set_tau(const double& value) { tau_ = value; }

  [[nodiscard]] double get_alpha() const { return alpha_; }
  void set_alpha(const double& value) {
    alpha_ = value;
    kernel_ = LocationPairTable{};
  }

  [[nodiscard]] double get_rho() const { return rho_; }
  void set_rho(const double& value) {
    rho_ = value;
    kernel_ = LocationPairTable{};
  }

  [[nodiscard]] double get_capital() const { return capital_; }
  void set_capital(const double& value) { capital_ = value; }

  [[nodiscard]] double get_penalty() const { return penalty_; }
  void set_penalty(const double& value) { penalty_ = value; }

private:
  double tau_;
  double alpha_;
  double rho_;
  double capital_;
  double penalty_;
  uint64_t number_of_locations_;

  // Borrowed, owned by SpatialSettings and outlives this object.
  const LocationPairTable* spatial_distance_{nullptr};

  std::vector<double> travel_;
  LocationPairTable kernel_;

  // District id per location, resolved once in prepare() instead of via a
  // string-keyed map lookup per destination.
  bool has_district_level_{false};
  std::vector<int> district_by_location_;

  void prepare_kernel();
  void prepare_districts();

public:
  explicit BurkinaFasoSM(double tau, double alpha, double rho, double capital, double penalty,
                         int number_of_locations, const LocationPairTable* spatial_distance)
      : tau_(tau),
        alpha_(alpha),
        rho_(rho),
        capital_(capital),
        penalty_(penalty),
        number_of_locations_(number_of_locations),
        spatial_distance_(spatial_distance) {}

  ~BurkinaFasoSM() override = default;

  void prepare() override;

  [[nodiscard]] DoubleVector get_v_relative_out_movement_to_destination(
      const int& from_location, const int& number_of_locations,
      const IntVector& v_number_of_residents_by_location) const override;
};
}  // namespace Spatial

#endif
