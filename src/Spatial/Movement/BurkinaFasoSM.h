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

#ifdef USE_DISTANCE_LUT
#include "Spatial/GIS/LocationPairTable.h"
#endif
#include "Spatial/SpatialModel.hxx"
#include "Utils/Helpers/NumberHelpers.h"
#include "Utils/TypeDef.h"

namespace Spatial {
class BurkinaFasoSM : public SpatialModel {
public:
  BurkinaFasoSM(const BurkinaFasoSM &) = delete;
  BurkinaFasoSM operator=(const BurkinaFasoSM &) = delete;
  BurkinaFasoSM(BurkinaFasoSM &&) = delete;
  BurkinaFasoSM &operator=(BurkinaFasoSM &&) = delete;

  [[nodiscard]] double get_tau() const { return tau_; }
  void set_tau(const double &value) { tau_ = value; }

  [[nodiscard]] double get_alpha() const { return alpha_; }
  void set_alpha(const double &value) {
    alpha_ = value;
#ifdef USE_DISTANCE_LUT
    kernel_lut_ = LocationPairTable{};
#endif
  }

  [[nodiscard]] double get_rho() const { return rho_; }
  void set_rho(const double &value) {
    rho_ = value;
#ifdef USE_DISTANCE_LUT
    kernel_lut_ = LocationPairTable{};
#endif
  }

  [[nodiscard]] double get_capital() const { return capital_; }
  void set_capital(const double &value) { capital_ = value; }

  [[nodiscard]] double get_penalty() const { return penalty_; }
  void set_penalty(const double &value) { penalty_ = value; }

private:
  double tau_;
  double alpha_;
  double rho_;
  double capital_;
  double penalty_;
  uint64_t number_of_locations_;
  std::vector<std::vector<double>> spatial_distance_matrix_;

  // These variables are computed when prepare() is called.
  std::vector<double> travel_;
#ifndef USE_DISTANCE_LUT
  std::vector<std::vector<double>> kernel_;
#else
  LocationPairTable constructor_distance_lut_;
  LocationPairTable kernel_lut_;
  bool has_district_level_{false};
  std::vector<int> district_by_location_;
#endif

  void prepare_kernel();
#ifdef USE_DISTANCE_LUT
  void prepare_districts();
#endif

public:
  explicit BurkinaFasoSM(double tau, double alpha, double rho, double capital, double penalty,
                         int number_of_locations,
                         std::vector<std::vector<double>> spatial_distance_matrix);

  ~BurkinaFasoSM() override = default;

  void prepare() override;

  // Public API intentionally remains identical to 500054a in both build modes.
  [[nodiscard]] DoubleVector get_v_relative_out_movement_to_destination(
      const int &from_location, const int &number_of_locations,
      const DoubleVector &relative_distance_vector,
      const IntVector &v_number_of_residents_by_location) const override;
};
}  // namespace Spatial

#endif
