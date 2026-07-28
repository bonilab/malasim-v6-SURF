#include <gtest/gtest.h>

#include <cmath>
#include <fstream>
#include <memory>

#include "Environment/SeasonalEquation.h"
#include "Environment/SeasonalRainfall.h"
#include "Simulation/Model.h"
#include "Configuration/Config.h"

TEST(SeasonalEquationTest, BuildsOnePeriodPerLocationAndClampsNegativeAmplitude) {
  SeasonalEquation equation;
  equation.set_raster_base({0.4});
  equation.set_raster_A({1.0});
  equation.set_raster_B({1.0});
  equation.set_raster_phi({182});

  equation.build(2);

  EXPECT_EQ(equation.get_base().size(), 2);
  EXPECT_DOUBLE_EQ(equation.get_seasonal_factor(date::sys_days{date::year{2024} / 1 / 1}, 0),
                   0.4);
  EXPECT_DOUBLE_EQ(equation.get_seasonal_factor(date::sys_days{date::year{2024} / 1 / 1}, 1),
                   0.4);
}

TEST(SeasonalEquationTest, UpdatesLocationsMatchingTheSourceReference) {
  SeasonalEquation equation;
  equation.set_base({0.2, 0.2, 0.9});
  equation.set_A({0.3, 0.3, 0.8});
  equation.set_B({1.0, 1.0, 2.0});
  equation.set_phi({10, 10, 20});
  equation.set_reference_base({0.2, 0.2, 0.9});
  equation.set_reference_A({0.3, 0.3, 0.8});
  equation.set_reference_B({1.0, 1.0, 2.0});
  equation.set_reference_phi({10, 10, 20});

  equation.update_seasonality(0, 2);

  EXPECT_DOUBLE_EQ(equation.get_base()[0], 0.9);
  EXPECT_DOUBLE_EQ(equation.get_A()[1], 0.8);
  EXPECT_DOUBLE_EQ(equation.get_B()[2], 2.0);
}

TEST(SeasonalEquationTest, RasterBuildFailsWithoutEcoclimaticRaster) {
  Model::set_config(std::make_unique<Config>());
  SeasonalEquation equation;
  equation.set_raster(true);
  equation.set_raster_base({0.2});
  equation.set_raster_A({0.3});
  equation.set_raster_B({1.0});
  equation.set_raster_phi({10});

  EXPECT_THROW(equation.build(1), std::invalid_argument);
  Model::get_instance()->release();
}

class SeasonalRainfallTest : public ::testing::Test {
protected:
  void TearDown() override { std::remove(filename_.c_str()); }
  const std::string filename_ = "seasonal_rainfall_test.txt";
};

TEST_F(SeasonalRainfallTest, ReadsAndReturnsDailyFactors) {
  std::ofstream output(filename_);
  output << "0.1 0.2 0.3";
  output.close();

  SeasonalRainfall rainfall;
  rainfall.set_filename(filename_);
  rainfall.set_period(3);

  EXPECT_NO_THROW(rainfall.build());
  EXPECT_DOUBLE_EQ(rainfall.get_seasonal_factor(date::sys_days{date::year{2023} / 1 / 1}, 0), 0.1);
  EXPECT_DOUBLE_EQ(rainfall.get_seasonal_factor(date::sys_days{date::year{2023} / 1 / 3}, 0), 0.3);
}

TEST_F(SeasonalRainfallTest, RejectsInvalidFileDataAndPeriod) {
  {
    std::ofstream output(filename_);
    output << "0.1 1.1";
  }
  SeasonalRainfall rainfall;
  rainfall.set_filename(filename_);
  rainfall.set_period(3);
  EXPECT_THROW(rainfall.build(), std::runtime_error);

  {
    std::ofstream output(filename_);
    output << "0.1";
  }
  EXPECT_THROW(rainfall.build(), std::invalid_argument);
}

TEST_F(SeasonalRainfallTest, RebuildingReplacesPreviouslyReadFactors) {
  {
    std::ofstream output(filename_);
    output << "0.1 0.2";
  }
  SeasonalRainfall rainfall;
  rainfall.set_filename(filename_);
  rainfall.set_period(2);
  ASSERT_NO_THROW(rainfall.build());

  {
    std::ofstream output(filename_);
    output << "0.7 0.8";
  }
  ASSERT_NO_THROW(rainfall.build());
  EXPECT_DOUBLE_EQ(rainfall.get_seasonal_factor(date::sys_days{date::year{2023} / 1 / 1}, 0), 0.7);
}

TEST_F(SeasonalRainfallTest, RejectsMissingAndEmptyFiles) {
  SeasonalRainfall rainfall;
  rainfall.set_filename("missing-seasonal-rainfall.txt");
  rainfall.set_period(1);
  EXPECT_THROW(rainfall.build(), std::runtime_error);

  const std::string empty_filename = "empty-seasonal-rainfall.txt";
  { std::ofstream output(empty_filename); }
  rainfall.set_filename(empty_filename);
  EXPECT_THROW(rainfall.build(), std::runtime_error);
  std::remove(empty_filename.c_str());
}
