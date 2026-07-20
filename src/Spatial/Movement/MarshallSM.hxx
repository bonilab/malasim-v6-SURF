/*
 * MarshallSM.hxx
 *
 * Gravity model for human migration based upon a distance kernel function.
 *
 * Marshall et al., 2018
 */
#ifndef MARSHALLSM_HXX
#define MARSHALLSM_HXX

#include "Spatial/GIS/GridPairTable.h"
#include "Spatial/SpatialModel.hxx"
#include "Utils/Helpers/NumberHelpers.h"
#include "Utils/TypeDef.h"

namespace Spatial {
class MarshallSM : public SpatialModel {
public:
  // Disallow copy
  MarshallSM(const MarshallSM &) = delete;
  MarshallSM &operator=(const MarshallSM &) = delete;

  // Disallow move
  MarshallSM(MarshallSM &&) = delete;
  MarshallSM &operator=(MarshallSM &&) = delete;

  [[nodiscard]] double get_tau() const { return tau_; }
  void set_tau(const double &value) { tau_ = value; }

  [[nodiscard]] double get_alpha() const { return alpha_; }
  void set_alpha(const double &value) {
    alpha_ = value;
    kernel_lut_ = GridPairTable{};
  }

  [[nodiscard]] double get_rho() const { return log_rho_; }
  void set_log_rho(const double &value) {
    log_rho_ = value;
    kernel_lut_ = GridPairTable{};
  }

  // Existing public data members are retained for source/API compatibility.
  double tau_;
  double alpha_;
  double log_rho_;
  int number_of_locations_;
  double** kernel = nullptr;

  explicit MarshallSM(double tau, double alpha, double log_rho, int number_of_locations);

  ~MarshallSM() override;

  // Retained as public because it was public in 500054a.
  void prepare_kernel();
  void prepare() override;

  // Public API intentionally remains identical to 500054a in both build modes.
  [[nodiscard]] DoubleVector get_v_relative_out_movement_to_destination(
      const int &from_location,
      const int &number_of_locations,
      const DoubleVector &relative_distance_vector,
      const IntVector &v_number_of_residents_by_location) const override;

private:
  void release_dense_kernel();

  GridPairTable kernel_lut_;
};
}  // namespace Spatial

#endif
