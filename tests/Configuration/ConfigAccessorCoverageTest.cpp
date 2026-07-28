#include <gtest/gtest.h>

#include "Configuration/Config.h"

TEST(ConfigAccessorCoverageTest, ExposesDrugAndRaptSettings) {
  Config config;
  EXPECT_NE(&config.get_drug_parameters(), nullptr);
  EXPECT_NE(&config.get_rapt_settings(), nullptr);
}

