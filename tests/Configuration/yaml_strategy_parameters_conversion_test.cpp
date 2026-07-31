#include <gtest/gtest.h>
#include <yaml-cpp/yaml.h>

#include "Configuration/StrategyParameters.h"

// Helper function to compare two std::vector<std::vector<double>>
bool inline compareNestedVectors(const std::vector<std::vector<double>> &vec1,
                                 const std::vector<std::vector<double>> &vec2) {
  if (vec1.size() != vec2.size()) return false;
  for (size_t i = 0; i < vec1.size(); ++i) {
    if (vec1[i] != vec2[i]) return false;  // Compare inner vectors
  }
  return true;
}

class StrategyParametersTest : public ::testing::Test {
protected:
  StrategyParameters strategy_parameters;

  void SetUp() override {
    // Setup StrategyInfo
    StrategyParameters::StrategyInfo strategy_info;
    strategy_info.set_name("SP-AQ-CQ-AL-MFTStrategy");
    strategy_info.set_type("MFT");
    strategy_info.set_therapy_ids({5, 2, 12, 6});
    strategy_info.set_distribution({0.3, 0.3, 0.3, 0.1});

    // Set start_distribution_by_location and peak_distribution_by_location
    strategy_info.set_start_distribution_by_location({{0.1, 0.2}, {0.3, 0.4}});
    strategy_info.set_peak_distribution_by_location({{0.5, 0.6}, {0.7, 0.8}});
    strategy_info.set_public_strategy_id(1);
    strategy_info.set_private_strategy_id(0);
    strategy_info.set_start_public_share(0.2);
    strategy_info.set_peak_public_share(0.8);
    strategy_info.set_start_public_share_by_location({0.2, 0.4});
    strategy_info.set_peak_public_share_by_location({0.8, 0.6});

    // Setup MassDrugAdministration
    StrategyParameters::MassDrugAdministration mda_info;
    mda_info.set_enable(false);
    mda_info.set_mda_therapy_id(8);
    mda_info.set_age_bracket_prob_individual_present_at_mda({10, 40});
    mda_info.set_mean_prob_individual_present_at_mda({0.85, 0.75, 0.85});
    mda_info.set_sd_prob_individual_present_at_mda({0.3, 0.3, 0.3});

    // Set StrategyParameters with a map
    std::map<int, StrategyParameters::StrategyInfo> strategy_db;
    strategy_db[0] = strategy_info;  // Using integer key 0 for the map

    strategy_parameters.set_strategy_db_raw(strategy_db);
    strategy_parameters.set_initial_strategy_id(15);
    strategy_parameters.set_second_line_strategy_id(2);
    strategy_parameters.set_mass_drug_administration(mda_info);
  }
};

// Test encoding functionality
TEST_F(StrategyParametersTest, EncodeStrategyParameters) {
  YAML::Node node = YAML::convert<StrategyParameters>::encode(strategy_parameters);

  EXPECT_EQ(node["initial_strategy_id"].as<int>(), 15);
  EXPECT_EQ(node["second_line_strategy_id"].as<int>(), 2);
  EXPECT_EQ(node["strategy_db"][0]["name"].as<std::string>(), "SP-AQ-CQ-AL-MFTStrategy");
  EXPECT_EQ(node["strategy_db"][0]["therapy_ids"].as<std::vector<int>>(),
            std::vector<int>({5, 2, 12, 6}));

  // Check the new start_distribution_by_location and peak_distribution_by_location fields
  auto start_distribution = node["strategy_db"][0]["start_distribution_by_location"]
                                .as<std::vector<std::vector<double>>>();
  auto peak_distribution = node["strategy_db"][0]["peak_distribution_by_location"]
                               .as<std::vector<std::vector<double>>>();
  EXPECT_TRUE(compareNestedVectors(start_distribution, {{0.1, 0.2}, {0.3, 0.4}}));
  EXPECT_TRUE(compareNestedVectors(peak_distribution, {{0.5, 0.6}, {0.7, 0.8}}));
  EXPECT_EQ(node["strategy_db"][0]["public_strategy_id"].as<int>(), 1);
  EXPECT_EQ(node["strategy_db"][0]["private_strategy_id"].as<int>(), 0);
  EXPECT_DOUBLE_EQ(node["strategy_db"][0]["start_public_share"].as<double>(), 0.2);
  EXPECT_DOUBLE_EQ(node["strategy_db"][0]["peak_public_share"].as<double>(), 0.8);
  EXPECT_EQ(node["strategy_db"][0]["start_public_share_by_location"].as<std::vector<double>>(),
            (std::vector<double>{0.2, 0.4}));

  EXPECT_EQ(node["mass_drug_administration"]["enable"].as<bool>(), false);
  EXPECT_EQ(node["mass_drug_administration"]["age_bracket_prob_individual_present_at_mda"]
                .as<std::vector<int>>(),
            std::vector<int>({10, 40}));
  EXPECT_EQ(node["mass_drug_administration"]["mean_prob_individual_present_at_mda"]
                .as<std::vector<double>>(),
            std::vector<double>({0.85, 0.75, 0.85}));
}

// Test decoding functionality
TEST_F(StrategyParametersTest, DecodeStrategyParameters) {
  YAML::Node node;
  node["initial_strategy_id"] = 15;
  node["second_line_strategy_id"] = 2;

  node["strategy_db"]["0"]["name"] = "SP-AQ-CQ-AL-MFTStrategy";
  node["strategy_db"]["0"]["type"] = "MFT";
  node["strategy_db"]["0"]["therapy_ids"] = std::vector<int>{5, 2, 12, 6};
  node["strategy_db"]["0"]["distribution"] = std::vector<double>{0.3, 0.3, 0.3, 0.1};

  // Set the new start_distribution_by_location and peak_distribution_by_location fields
  node["strategy_db"]["0"]["start_distribution_by_location"] =
      std::vector<std::vector<double>>{{0.1, 0.2}, {0.3, 0.4}};
  node["strategy_db"]["0"]["peak_distribution_by_location"] =
      std::vector<std::vector<double>>{{0.5, 0.6}, {0.7, 0.8}};
  node["strategy_db"]["0"]["public_strategy_id"] = 1;
  node["strategy_db"]["0"]["private_strategy_id"] = 0;
  node["strategy_db"]["0"]["start_public_share"] = 0.2;
  node["strategy_db"]["0"]["peak_public_share"] = 0.8;
  node["strategy_db"]["0"]["start_public_share_by_location"] = std::vector<double>{0.2, 0.4};
  node["strategy_db"]["0"]["peak_public_share_by_location"] = std::vector<double>{0.8, 0.6};

  node["mass_drug_administration"]["enable"] = false;
  node["mass_drug_administration"]["mda_therapy_id"] = 8;
  node["mass_drug_administration"]["age_bracket_prob_individual_present_at_mda"] =
      std::vector<int>{10, 40};
  node["mass_drug_administration"]["mean_prob_individual_present_at_mda"] =
      std::vector<double>{0.85, 0.75, 0.85};
  node["mass_drug_administration"]["sd_prob_individual_present_at_mda"] =
      std::vector<double>{0.3, 0.3, 0.3};

  StrategyParameters decoded_parameters;
  EXPECT_NO_THROW(YAML::convert<StrategyParameters>::decode(node, decoded_parameters));

  EXPECT_EQ(decoded_parameters.get_initial_strategy_id(), 15);
  EXPECT_EQ(decoded_parameters.get_second_line_strategy_id(), 2);
  EXPECT_EQ(decoded_parameters.get_strategy_db_raw().at(0).get_name(), "SP-AQ-CQ-AL-MFTStrategy");
  EXPECT_EQ(decoded_parameters.get_strategy_db_raw().at(0).get_therapy_ids(),
            std::vector<int>({5, 2, 12, 6}));

  // Check the decoded values for the new fields
  EXPECT_TRUE(compareNestedVectors(
      decoded_parameters.get_strategy_db_raw().at(0).get_start_distribution_by_location(),
      {{0.1, 0.2}, {0.3, 0.4}}));
  EXPECT_TRUE(compareNestedVectors(
      decoded_parameters.get_strategy_db_raw().at(0).get_peak_distribution_by_location(),
      {{0.5, 0.6}, {0.7, 0.8}}));
  EXPECT_EQ(decoded_parameters.get_strategy_db_raw().at(0).get_public_strategy_id(), 1);
  EXPECT_EQ(decoded_parameters.get_strategy_db_raw().at(0).get_private_strategy_id(), 0);
  EXPECT_DOUBLE_EQ(decoded_parameters.get_strategy_db_raw().at(0).get_start_public_share(), 0.2);
  EXPECT_DOUBLE_EQ(decoded_parameters.get_strategy_db_raw().at(0).get_peak_public_share(), 0.8);
  EXPECT_EQ(decoded_parameters.get_strategy_db_raw().at(0).get_peak_public_share_by_location(),
            (std::vector<double>{0.8, 0.6}));

  EXPECT_EQ(decoded_parameters.get_mda().get_enable(), false);
  EXPECT_EQ(decoded_parameters.get_mda().get_age_bracket_prob_individual_present_at_mda(),
            std::vector<int>({10, 40}));
  EXPECT_EQ(decoded_parameters.get_mda().get_mean_prob_individual_present_at_mda(),
            std::vector<double>({0.85, 0.75, 0.85}));
}

TEST_F(StrategyParametersTest, DecodeStrategyParametersDefaultsSecondLineStrategyToDisabled) {
  YAML::Node node;
  node["initial_strategy_id"] = 15;
  node["strategy_db"]["0"]["name"] = "SP-AQ-CQ-AL-MFTStrategy";
  node["strategy_db"]["0"]["type"] = "MFT";
  node["strategy_db"]["0"]["therapy_ids"] = std::vector<int>{5, 2, 12, 6};
  node["strategy_db"]["0"]["distribution"] = std::vector<double>{0.3, 0.3, 0.3, 0.1};
  node["mass_drug_administration"]["enable"] = false;
  node["mass_drug_administration"]["mda_therapy_id"] = 8;
  node["mass_drug_administration"]["age_bracket_prob_individual_present_at_mda"] =
      std::vector<int>{10, 40};
  node["mass_drug_administration"]["mean_prob_individual_present_at_mda"] =
      std::vector<double>{0.85, 0.75, 0.85};
  node["mass_drug_administration"]["sd_prob_individual_present_at_mda"] =
      std::vector<double>{0.3, 0.3, 0.3};

  StrategyParameters decoded_parameters;
  ASSERT_TRUE(YAML::convert<StrategyParameters>::decode(node, decoded_parameters));
  EXPECT_EQ(decoded_parameters.get_second_line_strategy_id(), -1);
}

// Test for decoding with missing fields
TEST_F(StrategyParametersTest, DecodeStrategyParametersMissingField) {
  YAML::Node node;
  node["initial_strategy_id"] = 15;  // Missing other fields

  StrategyParameters decoded_parameters;
  EXPECT_THROW(YAML::convert<StrategyParameters>::decode(node, decoded_parameters),
               std::runtime_error);
}

// ---------------------------------------------------------------------------
// DistrictMFT definitions round-trip
//
// StrategyInfo previously had no field for `definitions`, so a DistrictMFT
// strategy encoded to name + type only. It went unnoticed because
// StrategyParameters::process_config builds strategies from the raw YAML node
// kept by set_node(), not from strategy_db_raw_, so the simulation itself never
// read the decoded struct.
// ---------------------------------------------------------------------------

namespace {
YAML::Node make_district_mft_node() {
  YAML::Node node;
  node["initial_strategy_id"] = 0;
  node["strategy_db"]["0"]["name"] = "AngolaDistrictMFT-2021";
  node["strategy_db"]["0"]["type"] = "DistrictMFT";
  node["strategy_db"]["0"]["definitions"]["0"]["district_ids"] = std::vector<int>{4, 12, 13};
  node["strategy_db"]["0"]["definitions"]["0"]["therapy_ids"] =
      std::vector<int>{6, 8, 2, 0, 12, 14, 5};
  node["strategy_db"]["0"]["definitions"]["0"]["distribution"] =
      std::vector<double>{0.765, 0.085, 0.00225, 0.06, 0.00075, 0.027, 0.06};
  node["strategy_db"]["0"]["definitions"]["1"]["district_ids"] = std::vector<int>{2, 7};
  node["strategy_db"]["0"]["definitions"]["1"]["therapy_ids"] = std::vector<int>{6, 8};
  node["strategy_db"]["0"]["definitions"]["1"]["distribution"] = std::vector<double>{0.5, 0.5};

  node["mass_drug_administration"]["enable"] = false;
  node["mass_drug_administration"]["mda_therapy_id"] = 8;
  node["mass_drug_administration"]["age_bracket_prob_individual_present_at_mda"] =
      std::vector<int>{10, 40};
  node["mass_drug_administration"]["mean_prob_individual_present_at_mda"] =
      std::vector<double>{0.85, 0.75, 0.85};
  node["mass_drug_administration"]["sd_prob_individual_present_at_mda"] =
      std::vector<double>{0.3, 0.3, 0.3};
  return node;
}
}  // namespace

TEST_F(StrategyParametersTest, DecodeDistrictMftDefinitions) {
  StrategyParameters decoded;
  ASSERT_TRUE(YAML::convert<StrategyParameters>::decode(make_district_mft_node(), decoded));

  const auto &definitions = decoded.get_strategy_db_raw().at(0).get_definitions();
  ASSERT_EQ(definitions.size(), 2);

  EXPECT_EQ(definitions.at(0).district_ids, (std::vector<int>{4, 12, 13}));
  EXPECT_EQ(definitions.at(0).therapy_ids, (std::vector<int>{6, 8, 2, 0, 12, 14, 5}));
  ASSERT_EQ(definitions.at(0).distribution.size(), 7);
  EXPECT_DOUBLE_EQ(definitions.at(0).distribution[0], 0.765);
  EXPECT_DOUBLE_EQ(definitions.at(0).distribution[6], 0.06);

  EXPECT_EQ(definitions.at(1).district_ids, (std::vector<int>{2, 7}));
  EXPECT_EQ(definitions.at(1).therapy_ids, (std::vector<int>{6, 8}));
  EXPECT_EQ(definitions.at(1).distribution, (std::vector<double>{0.5, 0.5}));
}

TEST_F(StrategyParametersTest, EncodeDistrictMftDefinitions) {
  StrategyParameters decoded;
  ASSERT_TRUE(YAML::convert<StrategyParameters>::decode(make_district_mft_node(), decoded));

  const YAML::Node encoded = YAML::convert<StrategyParameters>::encode(decoded);

  // Both containers must be maps keyed by id. Checking only node[0] would pass
  // for a sequence too, which is exactly how the encode defect went unnoticed:
  // an integer key through operator[] on an empty node appends to a sequence
  // and throws the keys away.
  ASSERT_TRUE(encoded["strategy_db"].IsMap()) << "strategy_db must encode as a map, not a sequence";

  const auto definitions = encoded["strategy_db"][0]["definitions"];
  ASSERT_TRUE(definitions);
  ASSERT_TRUE(definitions.IsMap()) << "definitions must encode as a map, not a sequence";
  ASSERT_EQ(definitions.size(), 2);

  // Look the entries up by key rather than by position.
  EXPECT_EQ(definitions[0]["district_ids"].as<std::vector<int>>(), (std::vector<int>{4, 12, 13}));
  EXPECT_EQ(definitions[0]["therapy_ids"].as<std::vector<int>>(),
            (std::vector<int>{6, 8, 2, 0, 12, 14, 5}));
  EXPECT_EQ(definitions[1]["distribution"].as<std::vector<double>>(),
            (std::vector<double>{0.5, 0.5}));

  // The keys themselves must survive, which is the part a sequence loses.
  std::vector<int> keys;
  for (const auto &element : definitions) { keys.push_back(element.first.as<int>()); }
  EXPECT_EQ(keys, (std::vector<int>{0, 1}));
}

TEST_F(StrategyParametersTest, DistrictMftDefinitionsSurviveRoundTrip) {
  // decode -> encode -> decode must preserve every district assignment. This is
  // the property that was broken: the first encode dropped `definitions`
  // entirely and the second decode produced an empty map.
  StrategyParameters first;
  ASSERT_TRUE(YAML::convert<StrategyParameters>::decode(make_district_mft_node(), first));

  const YAML::Node encoded = YAML::convert<StrategyParameters>::encode(first);

  StrategyParameters second;
  ASSERT_TRUE(YAML::convert<StrategyParameters>::decode(encoded, second));

  const auto &before = first.get_strategy_db_raw().at(0).get_definitions();
  const auto &after = second.get_strategy_db_raw().at(0).get_definitions();
  ASSERT_EQ(after.size(), before.size());
  for (const auto &[key, definition] : before) {
    ASSERT_TRUE(after.contains(key)) << "definition " << key << " lost in round trip";
    EXPECT_EQ(after.at(key).district_ids, definition.district_ids);
    EXPECT_EQ(after.at(key).therapy_ids, definition.therapy_ids);
    EXPECT_EQ(after.at(key).distribution, definition.distribution);
  }
}

TEST_F(StrategyParametersTest, NonDistrictStrategiesOmitDefinitionsKey) {
  // The fixture strategy is an MFT with no definitions, so the key must not
  // appear at all rather than encoding as an empty map.
  const YAML::Node encoded = YAML::convert<StrategyParameters>::encode(strategy_parameters);
  EXPECT_FALSE(encoded["strategy_db"][0]["definitions"]);
  EXPECT_TRUE(strategy_parameters.get_strategy_db_raw().at(0).get_definitions().empty());
}

TEST_F(StrategyParametersTest, EncodeStrategyDbAsMapNotSequence) {
  // Regression guard for the encode defect that DistrictMftDefinitionsSurviveRoundTrip
  // exposed. This applies to every strategy, not just DistrictMFT: the fixture
  // strategy is a plain MFT and its key must survive encoding too.
  const YAML::Node encoded = YAML::convert<StrategyParameters>::encode(strategy_parameters);

  ASSERT_TRUE(encoded["strategy_db"].IsMap());
  ASSERT_FALSE(encoded["strategy_db"].IsSequence());

  std::vector<int> keys;
  for (const auto &element : encoded["strategy_db"]) { keys.push_back(element.first.as<int>()); }
  EXPECT_EQ(keys, (std::vector<int>{0}));
}

TEST_F(StrategyParametersTest, EncodedStrategyParametersAreDecodable) {
  // The narrowest statement of the bug: encode() produced something decode()
  // could not read back, so the two halves of the conversion disagreed.
  const YAML::Node encoded = YAML::convert<StrategyParameters>::encode(strategy_parameters);

  StrategyParameters round_tripped;
  ASSERT_NO_THROW(YAML::convert<StrategyParameters>::decode(encoded, round_tripped));
  EXPECT_EQ(round_tripped.get_initial_strategy_id(),
            strategy_parameters.get_initial_strategy_id());
  ASSERT_TRUE(round_tripped.get_strategy_db_raw().contains(0));
  EXPECT_EQ(round_tripped.get_strategy_db_raw().at(0).get_name(),
            strategy_parameters.get_strategy_db_raw().at(0).get_name());
  EXPECT_EQ(round_tripped.get_strategy_db_raw().at(0).get_therapy_ids(),
            strategy_parameters.get_strategy_db_raw().at(0).get_therapy_ids());
}
