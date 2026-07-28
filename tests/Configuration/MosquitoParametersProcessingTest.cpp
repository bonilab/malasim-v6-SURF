#include <gtest/gtest.h>

#include <vector>

#include "Configuration/MosquitoParameters.h"
#include "Configuration/Config.h"
#include "Simulation/Model.h"

namespace {

MosquitoParameters make_location_parameters(const std::vector<double>& ifr,
                                            const std::vector<int>& size) {
  MosquitoParameters::LocationBased location_based;
  location_based.set_interrupted_feeding_rate(ifr);
  location_based.set_prmc_size(size);

  MosquitoParameters::MosquitoConfig config;
  config.set_mode("location_based");
  config.set_location_based(location_based);

  MosquitoParameters parameters;
  parameters.set_mosquito_config(config);
  return parameters;
}

}  // namespace

TEST(MosquitoParametersProcessingTest, DistributesSingleLocationValues) {
  auto parameters = make_location_parameters({0.19}, {100});
  std::vector<Spatial::Location> locations(3);

  ASSERT_NO_THROW(parameters.process_config_using_locations(locations));
  for (const auto& location : locations) {
    EXPECT_DOUBLE_EQ(location.mosquito_ifr, 0.19);
    EXPECT_EQ(location.mosquito_size, 100);
  }
}

TEST(MosquitoParametersProcessingTest, AppliesPerLocationValues) {
  auto parameters = make_location_parameters({0.1, 0.2}, {80, 120});
  std::vector<Spatial::Location> locations(2);

  ASSERT_NO_THROW(parameters.process_config_using_locations(locations));
  EXPECT_DOUBLE_EQ(locations[0].mosquito_ifr, 0.1);
  EXPECT_EQ(locations[0].mosquito_size, 80);
  EXPECT_DOUBLE_EQ(locations[1].mosquito_ifr, 0.2);
  EXPECT_EQ(locations[1].mosquito_size, 120);
}

TEST(MosquitoParametersProcessingTest, RejectsMismatchedParameterArrays) {
  auto parameters = make_location_parameters({0.1}, {80, 120});
  std::vector<Spatial::Location> locations(2);

  EXPECT_THROW(parameters.process_config_using_locations(locations), std::invalid_argument);
}

TEST(MosquitoParametersProcessingTest, RejectsArraysThatDoNotMatchLocationCount) {
  auto parameters = make_location_parameters({0.1, 0.2, 0.3}, {80, 120, 160});
  std::vector<Spatial::Location> locations(2);

  EXPECT_THROW(parameters.process_config_using_locations(locations), std::invalid_argument);
}

TEST(MosquitoParametersProcessingTest, RejectsUnsupportedConfigurationMode) {
  MosquitoParameters::MosquitoConfig config;
  MosquitoParameters parameters;
  std::vector<Spatial::Location> locations(1);
  parameters.set_mosquito_config(config);

  EXPECT_NO_THROW(parameters.process_config_using_locations(locations));
  EXPECT_THROW(config.set_mode("unsupported"), std::invalid_argument);
}

TEST(MosquitoParametersProcessingTest, RejectsGridModeWithoutLoadedRaster) {
  auto* model = Model::get_instance();
  model->initialize();
  auto config = std::make_unique<Config>();
  config->get_spatial_settings().set_spatial_data(
      std::make_unique<SpatialData>(&config->get_spatial_settings()));
  model->set_config(std::move(config));

  MosquitoParameters::MosquitoConfig mosquito_config;
  mosquito_config.set_mode("grid_based");
  MosquitoParameters parameters;
  parameters.set_mosquito_config(mosquito_config);
  std::vector<Spatial::Location> locations(1);

  EXPECT_THROW(parameters.process_config_using_locations(locations), std::invalid_argument);
  model->release();
}
