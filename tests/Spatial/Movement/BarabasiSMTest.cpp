#include <gtest/gtest.h>

#include <memory>
#include <vector>

#include "Simulation/Model.h"
#include "Spatial/GIS/LocationPairTable.h"
#include "Spatial/Movement/BarabasiSM.hxx"
#include "Utils/Cli.h"
#include "fixtures/TestFileGenerators.h"

class BarabasiSMTest : public ::testing::Test {
protected:
  void SetUp() override {
    test_fixtures::setup_test_environment();
    Model::get_instance()->release();
    utils::Cli::MaSimAppInput cli_input;
    cli_input.input_path = "test_input.yml";
    Model::set_cli_input(cli_input);
    Model::get_instance()->initialize();

    r_g_0 = 1.0;
    beta_r = 0.5;
    kappa = 3.0;
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
    model = std::make_unique<Spatial::BarabasiSM>(r_g_0, beta_r, kappa, &spatial_distance);
    model->prepare();
  }

  double r_g_0;
  double beta_r;
  double kappa;
  LocationPairTable spatial_distance;
  std::unique_ptr<Spatial::BarabasiSM> model;
};

TEST_F(BarabasiSMTest, InitializeCorrectly) {
  EXPECT_DOUBLE_EQ(model->get_r_g_0(), r_g_0);
  EXPECT_DOUBLE_EQ(model->get_beta_r(), beta_r);
  EXPECT_DOUBLE_EQ(model->get_kappa(), kappa);
}

TEST_F(BarabasiSMTest, CalculateMovementToSameLocation) {
  const std::vector<int> residents_by_location = {1000, 2000, 3000};
  const auto movement =
      model->get_v_relative_out_movement_to_destination(0, 3, residents_by_location);
  EXPECT_DOUBLE_EQ(movement[0], 0.0);
}

TEST_F(BarabasiSMTest, CalculateMovementPattern) {
  const std::vector<double> distances = {0.0, 10.0, 20.0};
  const std::vector<int> residents_by_location = {1000, 2000, 3000};

  const auto movement =
      model->get_v_relative_out_movement_to_destination(0, 3, residents_by_location);

  const double expected_movement1 =
      std::pow(distances[1] + r_g_0, -beta_r) * std::exp(-distances[1] / kappa);
  const double expected_movement2 =
      std::pow(distances[2] + r_g_0, -beta_r) * std::exp(-distances[2] / kappa);

  EXPECT_NEAR(movement[1], expected_movement1, 1e-10);
  EXPECT_NEAR(movement[2], expected_movement2, 1e-10);
}

TEST_F(BarabasiSMTest, VerifyDecreasingPattern) {
  build_model({0.0, 5.0, 10.0, 20.0});
  const std::vector<int> residents_by_location = {1000, 2000, 3000, 4000};

  const auto movement =
      model->get_v_relative_out_movement_to_destination(0, 4, residents_by_location);

  EXPECT_GT(movement[1], movement[2]);
  EXPECT_GT(movement[2], movement[3]);
}
