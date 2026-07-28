#include <gtest/gtest.h>

#include "Utils/Helpers/NumberHelpers.h"

TEST(NumberHelpersTest, EvaluatesFloatingPointPredicates) {
  EXPECT_TRUE(NumberHelpers::is_equal(1.0, 1.0));
  EXPECT_FALSE(NumberHelpers::is_equal(1.0, 2.0));
  EXPECT_TRUE(NumberHelpers::is_not_equal(1.0, 2.0));
  EXPECT_TRUE(NumberHelpers::is_zero(0.0));
  EXPECT_FALSE(NumberHelpers::is_zero(1.0));
  EXPECT_TRUE(NumberHelpers::is_not_zero(1.0));
  EXPECT_FALSE(NumberHelpers::is_not_zero(0.0));
  EXPECT_TRUE(NumberHelpers::is_enot_qual(1.0, 2.0));
}

TEST(NumberHelpersTest, ConvertsNumbersAndGeneratesSeed) {
  EXPECT_EQ(NumberHelpers::number_to_string(42), "42");
  EXPECT_EQ(NumberHelpers::single_digit_number_to_char(7), '7');
  EXPECT_EQ(NumberHelpers::char_to_single_digit_number('7'), 7);
    EXPECT_NO_THROW({
      const auto seed = NumberHelpers::good_seed(1, 2, 3);
      (void)seed;
    });
}
