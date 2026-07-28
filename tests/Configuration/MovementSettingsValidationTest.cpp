#include <gtest/gtest.h>

#include "Configuration/MovementSettings.h"

TEST(MovementSettingsValidationTest, RejectsInvalidScalarSetterValues) {
  MovementSettings::BarabasiSM barabasi;
  EXPECT_THROW(barabasi.set_kappa(0), std::invalid_argument);

  MovementSettings::CirculationInfo circulation;
  EXPECT_THROW(circulation.set_max_relative_moving_value(-1), std::invalid_argument);
  EXPECT_THROW(circulation.set_number_of_moving_levels(1), std::runtime_error);
  EXPECT_THROW(circulation.set_circulation_percent(-0.1), std::invalid_argument);
  EXPECT_THROW(circulation.set_circulation_percent(100.1), std::invalid_argument);
}

TEST(MovementSettingsValidationTest, CoversRemainingMovementParameterSetters) {
  MovementSettings::WesolowskiSurfaceSM surface;
  surface.set_kappa(1.0);
  surface.set_alpha(2.0);
  surface.set_beta(3.0);
  surface.set_gamma(4.0);
  EXPECT_DOUBLE_EQ(surface.get_gamma(), 4.0);

  MovementSettings::MarshallSM marshall;
  marshall.set_tau(1.0);
  marshall.set_alpha(2.0);
  marshall.set_log_rho(3.0);
  EXPECT_DOUBLE_EQ(marshall.get_log_rho(), 3.0);

  MovementSettings::BurkinaFasoSM burkina;
  burkina.set_capital(14.0);
  burkina.set_penalty(12.0);
  EXPECT_DOUBLE_EQ(burkina.get_penalty(), 12.0);
}

TEST(MovementSettingsValidationTest, RejectsMissingNestedYamlModels) {
  MovementSettings::BarabasiSM barabasi;
  EXPECT_THROW(YAML::convert<MovementSettings::BarabasiSM>::decode(YAML::Node{}, barabasi),
               std::runtime_error);
  MovementSettings::WesolowskiSM wesolowski;
  EXPECT_THROW(YAML::convert<MovementSettings::WesolowskiSM>::decode(YAML::Node{}, wesolowski),
               std::runtime_error);
  MovementSettings::WesolowskiSurfaceSM surface;
  EXPECT_THROW(
      YAML::convert<MovementSettings::WesolowskiSurfaceSM>::decode(YAML::Node{}, surface),
      std::runtime_error);
  MovementSettings::MarshallSM marshall;
  EXPECT_THROW(YAML::convert<MovementSettings::MarshallSM>::decode(YAML::Node{}, marshall),
               std::runtime_error);
  MovementSettings::BurkinaFasoSM burkina;
  EXPECT_THROW(YAML::convert<MovementSettings::BurkinaFasoSM>::decode(YAML::Node{}, burkina),
               std::runtime_error);
  MovementSettings::MovingLevelDistributionGamma gamma;
  EXPECT_THROW(
      YAML::convert<MovementSettings::MovingLevelDistributionGamma>::decode(YAML::Node{}, gamma),
      std::runtime_error);
  MovementSettings::LengthOfStay stay;
  EXPECT_THROW(YAML::convert<MovementSettings::LengthOfStay>::decode(YAML::Node{}, stay),
               std::runtime_error);
  MovementSettings::MovingLevelDistributionExponential exponential;
  EXPECT_THROW(
      YAML::convert<MovementSettings::MovingLevelDistributionExponential>::decode(YAML::Node{}, exponential),
      std::runtime_error);
  MovementSettings::MovingLevelDistribution distribution;
  EXPECT_THROW(
      YAML::convert<MovementSettings::MovingLevelDistribution>::decode(YAML::Node{}, distribution),
      std::runtime_error);
  MovementSettings::SpatialModelSettings spatial_model;
  EXPECT_THROW(
      YAML::convert<MovementSettings::SpatialModelSettings>::decode(YAML::Node{}, spatial_model),
      std::runtime_error);
  MovementSettings::CirculationInfo circulation;
  EXPECT_THROW(
      YAML::convert<MovementSettings::CirculationInfo>::decode(YAML::Node{}, circulation),
      std::runtime_error);
  MovementSettings settings;
  EXPECT_THROW(YAML::convert<MovementSettings>::decode(YAML::Node{}, settings), std::runtime_error);
}
