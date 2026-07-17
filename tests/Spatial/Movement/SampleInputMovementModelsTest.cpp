#include <gtest/gtest.h>

#include <cmath>
#include <memory>
#include <vector>

#include "Simulation/Model.h"
#ifdef USE_DISTANCE_LUT
#include "Spatial/GIS/LocationPairTable.h"
#endif
#include "Spatial/Movement/BarabasiSM.hxx"
#include "Spatial/Movement/BurkinaFasoSM.h"
#include "Spatial/Movement/MarshallSM.hxx"
#include "Spatial/Movement/WesolowskiSM.hxx"
#include "Spatial/Movement/WesolowskiSurfaceSM.h"
#include "Utils/Cli.h"
#include "fixtures/TestFileGenerators.h"

namespace {
constexpr double kBarabasiRg0 = 5.8;
constexpr double kBarabasiBetaR = 1.65;
constexpr double kBarabasiKappa = 350.0;

constexpr double kWesolowskiKappa = 0.01093251;
constexpr double kWesolowskiAlpha = 0.22268982;
constexpr double kWesolowskiBeta = 0.14319618;
constexpr double kWesolowskiGamma = 0.83741484;

constexpr double kMarshallAlpha = 1.27;
constexpr double kMarshallLogRho = 0.54;
constexpr double kMarshallTau = 1.342;

constexpr double kBurkinaAlpha = 1.27;
constexpr double kBurkinaRho = 0.25;
constexpr double kBurkinaTau = 1.342;
constexpr double kBurkinaCapital = 14.0;
constexpr double kBurkinaPenalty = 12.0;

const DoubleVector kDistances = {0.0, 10.0, 20.0};
const IntVector kResidents = {1000, 2000, 3000};
const std::vector<std::vector<double>> kDistanceMatrix = {
    {0.0, 10.0, 20.0}, {10.0, 0.0, 15.0}, {20.0, 15.0, 0.0}};

const DoubleVector &distance_argument() {
#ifdef USE_DISTANCE_LUT
  static const DoubleVector empty_distance_vector;
  return empty_distance_vector;
#else
  return kDistances;
#endif
}
}  // namespace

class SampleInputMovementModelsTest : public ::testing::Test {
protected:
  void SetUp() override {
    test_fixtures::setup_test_environment();
    Model::get_instance()->release();
    utils::Cli::MaSimAppInput cli_input;
    cli_input.input_path = "test_input.yml";
    Model::set_cli_input(cli_input);
    Model::get_instance()->initialize();
#ifdef USE_DISTANCE_LUT
    // Replace the fixture table with deterministic distances used by these tests.
    Model::get_config()->get_spatial_settings().set_spatial_distance_lut(
        LocationPairTable::make_dense(kDistanceMatrix));
#endif
  }

  void TearDown() override {
    Model::get_instance()->release();
    test_fixtures::cleanup_test_files();
  }
};

TEST_F(SampleInputMovementModelsTest, BarabasiUsesSampleInputParameters) {
  Spatial::BarabasiSM model(kBarabasiRg0, kBarabasiBetaR, kBarabasiKappa);
  model.prepare();
  const auto movement = model.get_v_relative_out_movement_to_destination(
      0, static_cast<int>(kDistances.size()), distance_argument(), kResidents);

  EXPECT_DOUBLE_EQ(model.get_r_g_0(), kBarabasiRg0);
  EXPECT_DOUBLE_EQ(model.get_beta_r(), kBarabasiBetaR);
  EXPECT_DOUBLE_EQ(model.get_kappa(), kBarabasiKappa);
  EXPECT_DOUBLE_EQ(movement[0], 0.0);
  EXPECT_NEAR(movement[1],
              std::pow(kDistances[1] + kBarabasiRg0, -kBarabasiBetaR)
                  * std::exp(-kDistances[1] / kBarabasiKappa),
              1e-14);
}

TEST_F(SampleInputMovementModelsTest, WesolowskiUsesSampleInputParameters) {
  Spatial::WesolowskiSM model(kWesolowskiKappa, kWesolowskiAlpha,
                              kWesolowskiBeta, kWesolowskiGamma);
  model.prepare();
  const auto movement = model.get_v_relative_out_movement_to_destination(
      0, static_cast<int>(kDistances.size()), distance_argument(), kResidents);

  EXPECT_DOUBLE_EQ(model.get_kappa(), kWesolowskiKappa);
  EXPECT_DOUBLE_EQ(model.get_alpha(), kWesolowskiAlpha);
  EXPECT_DOUBLE_EQ(model.get_beta(), kWesolowskiBeta);
  EXPECT_DOUBLE_EQ(model.get_gamma(), kWesolowskiGamma);
  EXPECT_DOUBLE_EQ(movement[0], 0.0);
  EXPECT_NEAR(
      movement[1],
      kWesolowskiKappa
          * (std::pow(kResidents[0], kWesolowskiAlpha)
             * std::pow(kResidents[1], kWesolowskiBeta))
          / std::pow(kDistances[1], kWesolowskiGamma),
      1e-14);
}

TEST_F(SampleInputMovementModelsTest, WesolowskiSurfaceUsesSampleInputParameters) {
  Spatial::WesolowskiSurfaceSM model(kWesolowskiKappa, kWesolowskiAlpha,
                                     kWesolowskiBeta, kWesolowskiGamma,
                                     static_cast<int>(kDistances.size()));
  model.prepare();
  model.travel = {0.0, 0.5, 1.0};
  const auto movement = model.get_v_relative_out_movement_to_destination(
      0, static_cast<int>(kDistances.size()), distance_argument(), kResidents);

  const double gravity =
      kWesolowskiKappa
      * (std::pow(kResidents[0], kWesolowskiAlpha)
         * std::pow(kResidents[1], kWesolowskiBeta))
      / std::pow(kDistances[1], kWesolowskiGamma);
  EXPECT_DOUBLE_EQ(movement[0], 0.0);
  EXPECT_NEAR(movement[1], gravity / (1.0 + model.travel[0] + model.travel[1]), 1e-14);
}

TEST_F(SampleInputMovementModelsTest, MarshallUsesSampleInputParameters) {
  Spatial::MarshallSM model(kMarshallTau, kMarshallAlpha, kMarshallLogRho,
                            static_cast<int>(kDistances.size()), kDistanceMatrix);
  model.prepare();
  const auto movement = model.get_v_relative_out_movement_to_destination(
      0, static_cast<int>(kDistances.size()), distance_argument(), kResidents);

  EXPECT_DOUBLE_EQ(model.get_tau(), kMarshallTau);
  EXPECT_DOUBLE_EQ(model.get_alpha(), kMarshallAlpha);
  EXPECT_DOUBLE_EQ(model.get_rho(), kMarshallLogRho);
  EXPECT_DOUBLE_EQ(movement[0], 0.0);
  EXPECT_NEAR(movement[1],
              std::pow(kResidents[0], kMarshallTau)
                  * std::pow(1.0 + (kDistances[1] / kMarshallLogRho), -kMarshallAlpha),
              1e-10);
}

TEST_F(SampleInputMovementModelsTest, BurkinaFasoUsesSampleInputParameters) {
  Spatial::BurkinaFasoSM model(kBurkinaTau, kBurkinaAlpha, kBurkinaRho,
                               kBurkinaCapital, kBurkinaPenalty,
                               static_cast<int>(kDistances.size()), kDistanceMatrix);
  model.prepare();
  const auto movement = model.get_v_relative_out_movement_to_destination(
      0, static_cast<int>(kDistances.size()), distance_argument(), kResidents);

  EXPECT_DOUBLE_EQ(model.get_tau(), kBurkinaTau);
  EXPECT_DOUBLE_EQ(model.get_alpha(), kBurkinaAlpha);
  EXPECT_DOUBLE_EQ(model.get_rho(), kBurkinaRho);
  EXPECT_DOUBLE_EQ(model.get_capital(), kBurkinaCapital);
  EXPECT_DOUBLE_EQ(model.get_penalty(), kBurkinaPenalty);
  EXPECT_DOUBLE_EQ(movement[0], 0.0);
  EXPECT_NEAR(movement[1],
              std::pow(kResidents[0], kBurkinaTau)
                  * std::pow(1.0 + (kDistances[1] / kBurkinaRho), -kBurkinaAlpha),
              1e-10);
  ASSERT_EQ(movement.size(), kDistances.size());
  for (double value : movement) { EXPECT_GE(value, 0.0); }
}
