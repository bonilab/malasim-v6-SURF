#include <gtest/gtest.h>

#include <memory>
#include <vector>

#include "Simulation/Model.h"
#include "Spatial/GIS/LocationPairTable.h"
#include "Spatial/Movement/WesolowskiSM.hxx"
#include "Utils/Cli.h"
#include "Utils/TypeDef.h"
#include "fixtures/TestFileGenerators.h"

class WesolowskiSMTest : public ::testing::Test {
protected:
  void SetUp() override {
    test_fixtures::setup_test_environment();
    Model::get_instance()->release();
    utils::Cli::MaSimAppInput cli_input;
    cli_input.input_path = "test_input.yml";
    Model::set_cli_input(cli_input);
    Model::get_instance()->initialize();

    kappa = 2.0;
    alpha = 0.1;
    beta = 0.2;
    gamma = 0.3;
    build_model({0.0, 10.0, 20.0});
  }

  void TearDown() override {
    model.reset();
    Model::get_instance()->release();
    test_fixtures::cleanup_test_files();
  }

  void build_model(const std::vector<double>& distances_from_zero) {
    model.reset();
    const size_t n = distances_from_zero.size();
    std::vector<std::vector<double>> matrix(n, std::vector<double>(n, 1.0));
    for (size_t i = 0; i < n; ++i) {
      matrix[i][i] = 0.0;
      matrix[0][i] = distances_from_zero[i];
      matrix[i][0] = distances_from_zero[i];
    }
    spatial_distance = LocationPairTable::make_dense(std::move(matrix));
    model =
        std::make_unique<Spatial::WesolowskiSM>(kappa, alpha, beta, gamma, &spatial_distance);
    model->prepare();
  }

  double kappa;
  double alpha;
  double beta;
  double gamma;
  LocationPairTable spatial_distance;
  std::unique_ptr<Spatial::WesolowskiSM> model;
};

TEST_F(WesolowskiSMTest, InitializeCorrectly) {
  EXPECT_DOUBLE_EQ(model->get_kappa(), kappa);
  EXPECT_DOUBLE_EQ(model->get_alpha(), alpha);
  EXPECT_DOUBLE_EQ(model->get_beta(), beta);
  EXPECT_DOUBLE_EQ(model->get_gamma(), gamma);
}

TEST_F(WesolowskiSMTest, CalculateMovementToSameLocation) {
  const std::vector<int> residents_by_location = {1000, 2000, 3000};
  const auto movement =
      model->get_v_relative_out_movement_to_destination(0, 3, residents_by_location);
  EXPECT_DOUBLE_EQ(movement[0], 0.0);
}

TEST_F(WesolowskiSMTest, CalculateMovementPattern) {
  const std::vector<double> distances = {0.0, 10.0, 20.0};
  const std::vector<int> residents_by_location = {1000, 2000, 3000};

  const auto movement =
      model->get_v_relative_out_movement_to_destination(0, 3, residents_by_location);

  const double expected_movement1 =
      kappa
      * (std::pow(residents_by_location[0], alpha) * std::pow(residents_by_location[1], beta))
      / std::pow(distances[1], gamma);
  const double expected_movement2 =
      kappa
      * (std::pow(residents_by_location[0], alpha) * std::pow(residents_by_location[2], beta))
      / std::pow(distances[2], gamma);

  EXPECT_NEAR(movement[1], expected_movement1, 1e-10);
  EXPECT_NEAR(movement[2], expected_movement2, 1e-10);
}

TEST_F(WesolowskiSMTest, VerifyParameterEffects) {
  const std::vector<int> residents_by_location = {1000, 2000, 3000};
  const auto baseline =
      model->get_v_relative_out_movement_to_destination(0, 3, residents_by_location);

  model->set_kappa(kappa * 2.0);
  const auto updated =
      model->get_v_relative_out_movement_to_destination(0, 3, residents_by_location);

  EXPECT_NEAR(updated[1], baseline[1] * 2.0, 1e-10);
  EXPECT_NEAR(updated[2], baseline[2] * 2.0, 1e-10);
  model->set_kappa(kappa);
}

TEST_F(WesolowskiSMTest, DistanceEffectsOnMovement) {
  build_model({0.0, 5.0, 10.0, 20.0});
  const std::vector<int> residents_by_location = {1000, 1000, 1000, 1000};

  const auto movement =
      model->get_v_relative_out_movement_to_destination(0, 4, residents_by_location);

  EXPECT_GT(movement[1], movement[2]);
  EXPECT_GT(movement[2], movement[3]);
}

TEST_F(WesolowskiSMTest, PopulationEffectsOnMovement) {
  build_model({0.0, 10.0, 10.0, 10.0});
  const std::vector<int> residents_by_location = {1000, 1000, 2000, 4000};

  const auto movement =
      model->get_v_relative_out_movement_to_destination(0, 4, residents_by_location);

  EXPECT_GT(movement[2], movement[1]);
  EXPECT_GT(movement[3], movement[2]);
}
