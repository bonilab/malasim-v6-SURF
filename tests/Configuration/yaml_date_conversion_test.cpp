#include <date/date.h>
#include <gtest/gtest.h>
#include <yaml-cpp/yaml.h>

#include <stdexcept>

#include "Configuration/ModelSettings.h"

class YamlDateConversionTest : public ::testing::Test {
protected:
  date::year_month_day valid_date;
  YAML::Node valid_node;
  YAML::Node invalid_format_node;
  YAML::Node non_scalar_node;

  void SetUp() override {
    valid_date =
        date::year_month_day{date::year{2023}, date::month{10}, date::day{2}};
    valid_node = YAML::Node("2023/10/02");
    invalid_format_node = YAML::Node("invalid date format (should be yyyy/mm/dd)");
    non_scalar_node = YAML::Node(YAML::NodeType::Sequence);  // Non-scalar node
  }

  void TearDown() override {
    // Cleanup if necessary
  }
};

TEST_F(YamlDateConversionTest, EncodeValidDate) {
  YAML::Node node = YAML::convert<date::year_month_day>::encode(valid_date);

  // Check if the encoded date is correct
  EXPECT_EQ(node.as<std::string>(), "2023/10/02");
}

TEST_F(YamlDateConversionTest, DecodeValidDate) {
  date::year_month_day date;

  ASSERT_NO_THROW(
      { YAML::convert<date::year_month_day>::decode(valid_node, date); });

  EXPECT_EQ(date, valid_date);
}

TEST_F(YamlDateConversionTest, DecodeInvalidDateFormatThrows) {
  date::year_month_day date;

  // Expect an exception due to invalid date format
  EXPECT_THROW(
      {
        YAML::convert<date::year_month_day>::decode(invalid_format_node, date);
      },
      std::runtime_error);
}

TEST_F(YamlDateConversionTest, DecodeNonScalarThrows) {
  date::year_month_day date;

  // Expect an exception due to non-scalar node type
  EXPECT_THROW(
      { YAML::convert<date::year_month_day>::decode(non_scalar_node, date); },
      std::runtime_error);
}


