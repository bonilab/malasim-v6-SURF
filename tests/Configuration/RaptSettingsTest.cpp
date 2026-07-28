#include <gtest/gtest.h>
#include <yaml-cpp/yaml.h>

#include "Utils/YamlFile.h"
#include "date/date.h"

#include "Configuration/RaptSettings.h"

// RaptSettings is a plain config struct: getters/setters, a start-date
// computation, and YAML encode/decode. No Model dependency.
class RaptSettingsTest : public ::testing::Test {
protected:
  RaptSettings settings_;
};

TEST_F(RaptSettingsTest, DefaultsAreSensible) {
  EXPECT_TRUE(settings_.get_is_enabled());
  EXPECT_EQ(settings_.get_period(), 12);
  EXPECT_EQ(settings_.get_therapy_id(), 1);
  EXPECT_DOUBLE_EQ(settings_.get_compliance(), 0.7);
  EXPECT_EQ(settings_.get_age_start(), 18);
  EXPECT_EQ(settings_.get_start_date_as_day(), 0);
}

TEST_F(RaptSettingsTest, SettersRoundTrip) {
  settings_.set_is_enabled(false);
  settings_.set_period(6);
  settings_.set_therapy_id(3);
  settings_.set_compliance(0.95);
  settings_.set_age_start(5);
  const date::year_month_day input_date{date::year{2015}, date::month{5}, date::day{5}};
  settings_.set_start_date(input_date);

  EXPECT_FALSE(settings_.get_is_enabled());
  EXPECT_EQ(settings_.get_period(), 6);
  EXPECT_EQ(settings_.get_therapy_id(), 3);
  EXPECT_DOUBLE_EQ(settings_.get_compliance(), 0.95);
  EXPECT_EQ(settings_.get_age_start(), 5);
  EXPECT_EQ(settings_.get_start_date(), input_date);
}

TEST_F(RaptSettingsTest, ProcessConfigComputesPositiveOffsetWhenStartAfterBeginning) {
  settings_.set_start_date(date::year_month_day{date::year{2010}, date::month{1}, date::day{1}});
  const date::year_month_day input_beginning{date::year{2000}, date::month{1}, date::day{1}};

  settings_.process_config_with_starting_date(input_beginning);

  EXPECT_GT(settings_.get_start_date_as_day(), 0);
}

TEST_F(RaptSettingsTest, ProcessConfigComputesNegativeOffsetWhenStartBeforeBeginning) {
  settings_.set_start_date(date::year_month_day{date::year{2000}, date::month{1}, date::day{1}});
  const date::year_month_day input_beginning{date::year{2010}, date::month{1}, date::day{1}};

  settings_.process_config_with_starting_date(input_beginning);

  EXPECT_LT(settings_.get_start_date_as_day(), 0);
}

TEST_F(RaptSettingsTest, ProcessConfigZeroOffsetWhenSameDate) {
  const date::year_month_day input_beginning{date::year{2000}, date::month{1}, date::day{1}};
  settings_.set_start_date(input_beginning);

  settings_.process_config_with_starting_date(input_beginning);

  EXPECT_EQ(settings_.get_start_date_as_day(), 0);
}

TEST_F(RaptSettingsTest, DecodeAllFieldsWhenEnabled) {
  YAML::Node node;
  node["enabled"] = true;
  node["period"] = 24;
  node["therapy_id"] = 2;
  node["compliance"] = 0.8;
  node["age_start"] = 10;
  node["start_date"] = "2012/06/01";

  RaptSettings decoded;
  ASSERT_NO_THROW(YAML::convert<RaptSettings>::decode(node, decoded));
  EXPECT_TRUE(decoded.get_is_enabled());
  EXPECT_EQ(decoded.get_period(), 24);
  EXPECT_EQ(decoded.get_therapy_id(), 2);
  EXPECT_DOUBLE_EQ(decoded.get_compliance(), 0.8);
  EXPECT_EQ(decoded.get_age_start(), 10);
  EXPECT_EQ(decoded.get_start_date(),
            date::year_month_day(date::year{2012}, date::month{6}, date::day{1}));
}

TEST_F(RaptSettingsTest, DecodeDisabledSkipsOptionalFields) {
  YAML::Node node;
  node["enabled"] = false;

  RaptSettings decoded;
  ASSERT_NO_THROW(YAML::convert<RaptSettings>::decode(node, decoded));
  EXPECT_FALSE(decoded.get_is_enabled());
  // Defaults preserved when disabled and no fields supplied.
  EXPECT_EQ(decoded.get_period(), 12);
}

TEST_F(RaptSettingsTest, DecodeMissingEnabledThrows) {
  YAML::Node node;
  RaptSettings decoded;
  EXPECT_THROW(YAML::convert<RaptSettings>::decode(node, decoded), std::runtime_error);
}

TEST_F(RaptSettingsTest, DecodeMissingPeriodThrowsWhenEnabled) {
  YAML::Node node;
  node["enabled"] = true;
  node["therapy_id"] = 1;
  node["compliance"] = 0.7;
  node["age_start"] = 18;
  node["start_date"] = "2010/01/01";
  RaptSettings decoded;
  EXPECT_THROW(YAML::convert<RaptSettings>::decode(node, decoded), std::runtime_error);
}

TEST_F(RaptSettingsTest, DecodeMissingTherapyIdThrowsWhenEnabled) {
  YAML::Node node;
  node["enabled"] = true;
  node["period"] = 12;
  node["compliance"] = 0.7;
  node["age_start"] = 18;
  node["start_date"] = "2010/01/01";
  RaptSettings decoded;
  EXPECT_THROW(YAML::convert<RaptSettings>::decode(node, decoded), std::runtime_error);
}

TEST_F(RaptSettingsTest, DecodeMissingComplianceThrowsWhenEnabled) {
  YAML::Node node;
  node["enabled"] = true;
  node["period"] = 12;
  node["therapy_id"] = 1;
  node["age_start"] = 18;
  node["start_date"] = "2010/01/01";
  RaptSettings decoded;
  EXPECT_THROW(YAML::convert<RaptSettings>::decode(node, decoded), std::runtime_error);
}

TEST_F(RaptSettingsTest, DecodeMissingAgeStartThrowsWhenEnabled) {
  YAML::Node node;
  node["enabled"] = true;
  node["period"] = 12;
  node["therapy_id"] = 1;
  node["compliance"] = 0.7;
  node["start_date"] = "2010/01/01";
  RaptSettings decoded;
  EXPECT_THROW(YAML::convert<RaptSettings>::decode(node, decoded), std::runtime_error);
}

TEST_F(RaptSettingsTest, DecodeMissingStartDateThrowsWhenEnabled) {
  YAML::Node node;
  node["enabled"] = true;
  node["period"] = 12;
  node["therapy_id"] = 1;
  node["compliance"] = 0.7;
  node["age_start"] = 18;
  RaptSettings decoded;
  EXPECT_THROW(YAML::convert<RaptSettings>::decode(node, decoded), std::runtime_error);
}

TEST_F(RaptSettingsTest, EncodeRoundTripsThroughDecode) {
  settings_.set_is_enabled(true);
  settings_.set_period(36);
  settings_.set_therapy_id(4);
  settings_.set_compliance(0.55);
  settings_.set_age_start(21);
  settings_.set_start_date(date::year_month_day(date::year{2018}, date::month{3}, date::day{9}));

  const YAML::Node encoded = YAML::convert<RaptSettings>::encode(settings_);
  RaptSettings decoded;
  ASSERT_NO_THROW(YAML::convert<RaptSettings>::decode(encoded, decoded));

  EXPECT_EQ(decoded.get_period(), settings_.get_period());
  EXPECT_EQ(decoded.get_therapy_id(), settings_.get_therapy_id());
  EXPECT_DOUBLE_EQ(decoded.get_compliance(), settings_.get_compliance());
  EXPECT_EQ(decoded.get_age_start(), settings_.get_age_start());
  EXPECT_EQ(decoded.get_start_date(), settings_.get_start_date());
}
