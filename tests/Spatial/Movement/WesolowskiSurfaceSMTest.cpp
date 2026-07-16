#include <gtest/gtest.h>

#include <memory>

#include "Simulation/Model.h"
#include "Spatial/GIS/LocationPairTable.h"
#include "Spatial/Movement/WesolowskiSurfaceSM.h"
#include "Utils/Cli.h"
#include "Utils/TypeDef.h"
#include "fixtures/TestFileGenerators.h"

class WesolowskiSurfaceSMTest : public ::testing::Test {
protected:
  void SetUp() override {
    test_fixtures::setup_test_environment();

    test_fixtures::create_test_raster_2_locations("test_init_pop.asc", 1000.0);
    test_fixtures::create_test_raster_2_locations("test_beta.asc", 0.5);
    test_fixtures::create_test_raster_2_locations("test_treatment.asc", 0.6);
    test_fixtures::create_test_raster_2_locations("test_ecozone.asc", 1.0);
    test_fixtures::create_test_raster_2_locations("test_travel.asc", 0.2);

    Model::get_instance()->release();
    utils::Cli::MaSimAppInput cli_input;
    cli_input.input_path = "test_input.yml";
    Model::set_cli_input(cli_input);
    Model::get_instance()->initialize();

    kappa = 2.0;
    alpha = 0.1;
    beta = 0.2;
    gamma = 0.3;
    number_of_locations = 2;
    spatial_distance = LocationPairTable::make_dense({{0.0, 10.0}, {10.0, 0.0}});
    model = std::make_unique<Spatial::WesolowskiSurfaceSM>(
        kappa, alpha, beta, gamma, number_of_locations, &spatial_distance);
    model->prepare();
    model->travel = {0.2, 0.5};
  }

  void TearDown() override {
    model.reset();
    Model::get_instance()->release();
    test_fixtures::cleanup_test_files();
  }

  double kappa;
  double alpha;
  double beta;
  double gamma;
  int number_of_locations;
  LocationPairTable spatial_distance;
  std::unique_ptr<Spatial::WesolowskiSurfaceSM> model;
};

TEST_F(WesolowskiSurfaceSMTest, InitializeCorrectly) {
  EXPECT_DOUBLE_EQ(model->get_kappa(), kappa);
  EXPECT_DOUBLE_EQ(model->get_alpha(), alpha);
  EXPECT_DOUBLE_EQ(model->get_beta(), beta);
  EXPECT_DOUBLE_EQ(model->get_gamma(), gamma);
}

TEST_F(WesolowskiSurfaceSMTest, PrepareMethodWorks) {
  model->travel.clear();
  EXPECT_NO_THROW(model->prepare());
  EXPECT_FALSE(model->travel.empty());
  EXPECT_EQ(model->travel.size(), number_of_locations);
}

TEST_F(WesolowskiSurfaceSMTest, CalculateMovementToSameLocation) {
  const std::vector<int> residents_by_location = {1000, 2000};
  const auto movement = model->get_v_relative_out_movement_to_destination(
      0, number_of_locations, residents_by_location);
  EXPECT_DOUBLE_EQ(movement[0], 0.0);
}

TEST_F(WesolowskiSurfaceSMTest, CalculateMovementPattern) {
  const std::vector<int> residents_by_location = {1000, 2000};
  const auto movement = model->get_v_relative_out_movement_to_destination(
      0, number_of_locations, residents_by_location);

  const double probability =
      kappa
      * (std::pow(residents_by_location[0], alpha) * std::pow(residents_by_location[1], beta))
      / std::pow(10.0, gamma);
  const double expected = probability / (1.0 + model->travel[0] + model->travel[1]);

  EXPECT_NEAR(movement[1], expected, 1e-10);
}

TEST_F(WesolowskiSurfaceSMTest, TravelSurfacePenalty) {
  const std::vector<int> residents_by_location = {1000, 2000};
  const auto original = model->get_v_relative_out_movement_to_destination(
      0, number_of_locations, residents_by_location);

  const auto original_travel = model->travel;
  model->travel = {0.4, 1.0};
  const auto penalized = model->get_v_relative_out_movement_to_destination(
      0, number_of_locations, residents_by_location);

  EXPECT_LT(penalized[1], original[1]);
  model->travel = original_travel;
}

TEST_F(WesolowskiSurfaceSMTest, EmptyTravelSurfaceThrows) {
  const std::vector<int> residents_by_location = {1000, 2000};
  model->travel.clear();

#pragma clang diagnostic ignored "-Wunused-result"
  EXPECT_THROW(model->get_v_relative_out_movement_to_destination(
                   0, number_of_locations, residents_by_location),
               std::runtime_error);
#pragma clang diagnostic push
}
