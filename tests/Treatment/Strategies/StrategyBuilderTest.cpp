#include <gtest/gtest.h>
#include <memory>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>
#include <yaml-cpp/yaml.h>

#include "Treatment/Strategies/StrategyBuilder.h"
#include "Treatment/Strategies/IStrategy.h"
#include "Treatment/Strategies/SFTStrategy.h"
#include "Treatment/Strategies/MFTStrategy.h"
#include "Treatment/Strategies/CyclingStrategy.h"
#include "Treatment/Strategies/AdaptiveCyclingStrategy.h"
#include "Treatment/Strategies/MFTRebalancingStrategy.h"
#include "Treatment/Strategies/MFTMultiLocationStrategy.h"
#include "Treatment/Strategies/NestedMFTStrategy.h"
#include "Treatment/Strategies/NestedMFTMultiLocationStrategy.h"
#include "Treatment/Strategies/NovelDrugIntroductionStrategy.h"
#include "Treatment/Strategies/MFTAgeBasedStrategy.h"
#include "Treatment/Strategies/PublicPrivateStrategy.h"
#include "Treatment/Strategies/PublicPrivateMultiLocationStrategy.h"
#include "Treatment/Strategies/DistrictMftStrategy.h"
#include "Treatment/Therapies/Therapy.h"
#include "Treatment/Therapies/TherapyBuilder.h"
#include "Simulation/Model.h"
#include "Spatial/GIS/SpatialData.h"
#include "Utils/Cli.h"
#include "fixtures/TestFileGenerators.h"

class StrategyBuilderTest : public ::testing::Test {
protected:
  void SetUp() override {
    test_fixtures::setup_test_environment();
    Model::get_instance()->release();
    utils::Cli::MaSimAppInput cli_input;
    cli_input.input_path = "test_input.yml";
    Model::set_cli_input(cli_input);
    Model::get_instance()->initialize();
    
    // Use therapies from the therapy database
    therapies.clear();
    // Need at least 2 therapies
    for (int i = 0; i < 2 && i < Model::get_therapy_db().size(); i++) {
      therapies.push_back(Model::get_therapy_db()[i].get());
    }
    
    // Skip test if not enough therapies
    if (therapies.size() < 2) {
      GTEST_SKIP() << "Not enough therapies in the database to run this test";
    }
  }

  void TearDown() override {
    therapies.clear();
    test_fixtures::cleanup_test_files();
  }

  std::vector<Therapy*> therapies;
  
  // Helper to create a YAML node for strategy testing
  YAML::Node create_sft_strategy_node() {
    YAML::Node node;
    node["name"] = "Test SFT Strategy";
    node["type"] = "SFT";
    
    YAML::Node therapy_ids;
    therapy_ids.push_back(therapies[0]->get_id());
    node["therapy_ids"] = therapy_ids;
    
    return node;
  }
  
  YAML::Node create_mft_strategy_node() {
    YAML::Node node;
    node["name"] = "Test MFT Strategy";
    node["type"] = "MFT";
    
    YAML::Node therapy_ids;
    therapy_ids.push_back(therapies[0]->get_id());
    therapy_ids.push_back(therapies[1]->get_id());
    node["therapy_ids"] = therapy_ids;
    
    YAML::Node distributions;
    distributions.push_back(0.3);
    distributions.push_back(0.7);
    node["distributions"] = distributions;
    
    return node;
  }
  
  YAML::Node create_cycling_strategy_node() {
    YAML::Node node;
    node["name"] = "Test Cycling Strategy";
    node["type"] = "Cycling";
    
    YAML::Node therapy_ids;
    therapy_ids.push_back(therapies[0]->get_id());
    therapy_ids.push_back(therapies[1]->get_id());
    node["therapy_ids"] = therapy_ids;
    
    node["cycling_time"] = 30;
    
    return node;
  }
  
  YAML::Node create_adaptive_cycling_strategy_node() {
    YAML::Node node;
    node["name"] = "Test Adaptive Cycling Strategy";
    node["type"] = "AdaptiveCycling";
    
    YAML::Node therapy_ids;
    therapy_ids.push_back(therapies[0]->get_id());
    therapy_ids.push_back(therapies[1]->get_id());
    node["therapy_ids"] = therapy_ids;
    
    node["trigger_value"] = 0.15;
    node["delay_until_actual_trigger"] = 30;
    node["turn_off_days"] = 90;
    
    return node;
  }

  // The fixture raster produced by test_fixtures::create_test_district_raster
  // defines exactly three districts, with 1-based ids 1, 2 and 3. The
  // DistrictMFT builder demands that every district receives an assignment, so
  // these helpers exist to make the coverage explicit in each test.
  static constexpr int k_district_count = 3;

  YAML::Node create_district_mft_node(const std::string &definitions_yaml) {
    YAML::Node node;
    node["name"] = "Test District MFT";
    node["type"] = "DistrictMFT";
    node["definitions"] = YAML::Load(definitions_yaml);
    return node;
  }

  // A definition set covering all three districts. The first distribution is
  // taken verbatim from a production Angola configuration: its decimal values
  // sum to 1.0 but accumulate in binary floating point to 0.99999999..., which
  // is precisely what the previous static_cast<int>(sum) != 1 check rejected.
  YAML::Node create_full_coverage_district_mft_node() {
    return create_district_mft_node(R"(
      0:
        district_ids: [1, 2]
        therapy_ids: [6, 8, 2, 0, 12, 14, 5]
        distribution: [0.765, 0.085, 0.00225, 0.06, 0.00075, 0.027, 0.06]
      1:
        district_ids: [3]
        therapy_ids: [6, 8]
        distribution: [0.5, 0.5]
    )");
  }
};

TEST_F(StrategyBuilderTest, BuildSFTStrategy) {
  YAML::Node node = create_sft_strategy_node();
  
  // Build the strategy
  std::unique_ptr<IStrategy> strategy = StrategyBuilder::build(node, 1);
  
  // Check that we got the right type
  EXPECT_EQ(strategy->get_type(), IStrategy::StrategyType::SFT);
  
  // Check properties
  EXPECT_EQ(strategy->id, 1);
  EXPECT_EQ(strategy->name, "Test SFT Strategy");
  
  // Cast to SFTStrategy and check therapy list
  auto sft_strategy = dynamic_cast<SFTStrategy*>(strategy.get());
  ASSERT_NE(sft_strategy, nullptr);
  ASSERT_EQ(sft_strategy->get_therapy_list().size(), 1);
  EXPECT_EQ(sft_strategy->get_therapy_list()[0]->get_id(), therapies[0]->get_id());
}

TEST_F(StrategyBuilderTest, BuildMFTStrategy) {
  YAML::Node node = create_mft_strategy_node();
  
  // Build the strategy
  std::unique_ptr<IStrategy> strategy = StrategyBuilder::build(node, 2);
  
  // Check that we got the right type
  EXPECT_EQ(strategy->get_type(), IStrategy::StrategyType::MFT);
  
  // Check properties
  EXPECT_EQ(strategy->id, 2);
  EXPECT_EQ(strategy->name, "Test MFT Strategy");
  
  // Cast to MFTStrategy and check therapy list and distributions
  auto mft_strategy = dynamic_cast<MFTStrategy*>(strategy.get());
  ASSERT_NE(mft_strategy, nullptr);
  
  ASSERT_EQ(mft_strategy->therapy_list.size(), 2);
  EXPECT_EQ(mft_strategy->therapy_list[0]->get_id(), therapies[0]->get_id());
  EXPECT_EQ(mft_strategy->therapy_list[1]->get_id(), therapies[1]->get_id());
  
  ASSERT_EQ(mft_strategy->distribution.size(), 2);
  EXPECT_DOUBLE_EQ(mft_strategy->distribution[0], 0.3);
  EXPECT_DOUBLE_EQ(mft_strategy->distribution[1], 0.7);
}

TEST_F(StrategyBuilderTest, BuildCyclingStrategy) {
  YAML::Node node = create_cycling_strategy_node();
  
  // Build the strategy
  std::unique_ptr<IStrategy> strategy = StrategyBuilder::build(node, 3);
  
  // Check that we got the right type
  EXPECT_EQ(strategy->get_type(), IStrategy::StrategyType::Cycling);
  
  // Check properties
  EXPECT_EQ(strategy->id, 3);
  EXPECT_EQ(strategy->name, "Test Cycling Strategy");
  
  // Cast to CyclingStrategy and check therapy list and cycling time
  auto cycling_strategy = dynamic_cast<CyclingStrategy*>(strategy.get());
  ASSERT_NE(cycling_strategy, nullptr);
  
  ASSERT_EQ(cycling_strategy->therapy_list.size(), 2);
  EXPECT_EQ(cycling_strategy->therapy_list[0]->get_id(), therapies[0]->get_id());
  EXPECT_EQ(cycling_strategy->therapy_list[1]->get_id(), therapies[1]->get_id());
  
  EXPECT_EQ(cycling_strategy->cycling_time, 30);
}

TEST_F(StrategyBuilderTest, BuildAdaptiveCyclingStrategy) {
  YAML::Node node = create_adaptive_cycling_strategy_node();
  
  // Build the strategy
  std::unique_ptr<IStrategy> strategy = StrategyBuilder::build(node, 4);
  
  // Check that we got the right type
  EXPECT_EQ(strategy->get_type(), IStrategy::StrategyType::AdaptiveCycling);
  
  // Check properties
  EXPECT_EQ(strategy->id, 4);
  EXPECT_EQ(strategy->name, "Test Adaptive Cycling Strategy");
  
  // Cast to AdaptiveCyclingStrategy and check therapy list and parameters
  auto adaptive_cycling_strategy = dynamic_cast<AdaptiveCyclingStrategy*>(strategy.get());
  ASSERT_NE(adaptive_cycling_strategy, nullptr);
  
  ASSERT_EQ(adaptive_cycling_strategy->therapy_list.size(), 2);
  EXPECT_EQ(adaptive_cycling_strategy->therapy_list[0]->get_id(), therapies[0]->get_id());
  EXPECT_EQ(adaptive_cycling_strategy->therapy_list[1]->get_id(), therapies[1]->get_id());
  
  EXPECT_DOUBLE_EQ(adaptive_cycling_strategy->trigger_value, 0.15);
  EXPECT_EQ(adaptive_cycling_strategy->delay_until_actual_trigger, 30);
  EXPECT_EQ(adaptive_cycling_strategy->turn_off_days, 90);
}

TEST_F(StrategyBuilderTest, AddTherapies) {
  YAML::Node node;
  YAML::Node therapy_ids;
  therapy_ids.push_back(therapies[0]->get_id());
  therapy_ids.push_back(therapies[1]->get_id());
  node["therapy_ids"] = therapy_ids;
  
  // Create a strategy to add therapies to
  auto strategy = std::make_unique<SFTStrategy>();
  
  // Add therapies from the node
  StrategyBuilder::add_therapies(node, strategy.get());
  
  // Check therapies were added correctly
  ASSERT_EQ(strategy->get_therapy_list().size(), 2);
  EXPECT_EQ(strategy->get_therapy_list()[0]->get_id(), therapies[0]->get_id());
  EXPECT_EQ(strategy->get_therapy_list()[1]->get_id(), therapies[1]->get_id());
}

TEST_F(StrategyBuilderTest, AddDistributions) {
  YAML::Node node;
  YAML::Node distributions;
  distributions.push_back(0.2);
  distributions.push_back(0.3);
  distributions.push_back(0.5);
  node["distributions"] = distributions;
  
  // Create vector to add distributions to
  std::vector<double> dist_vector;
  
  // Add distributions from the node - must pass the sequence node, not the parent node
  StrategyBuilder::add_distributions(node["distributions"], dist_vector);
  
  // Check distributions were added correctly
  ASSERT_EQ(dist_vector.size(), 3);
  EXPECT_DOUBLE_EQ(dist_vector[0], 0.2);
  EXPECT_DOUBLE_EQ(dist_vector[1], 0.3);
  EXPECT_DOUBLE_EQ(dist_vector[2], 0.5);
}

TEST_F(StrategyBuilderTest, InvalidStrategyType) {
  YAML::Node node;
  node["name"] = "Invalid Strategy";
  node["type"] = "InvalidType";

  // An unrecognised type name is a configuration error, so build() reports it
  // as std::invalid_argument, consistent with every other field validation in
  // StrategyBuilder.
  //
  // This assertion previously passed for the wrong reason. build() looked the
  // name up with std::map::operator[], which default-inserts a zero-valued
  // entry, so "InvalidType" resolved to StrategyType::SFT (0) and the SFT
  // builder then died on the missing therapy_ids node with a YAML::InvalidNode
  // - which happens to derive from std::runtime_error. The test was green while
  // the type dispatch was silently broken.
  try {
    StrategyBuilder::build(node, 5);
    FAIL() << "Expected std::invalid_argument for an unrecognised strategy type";
  } catch (const std::invalid_argument &ex) {
    // The message must name the offending type, otherwise a typo in a large
    // strategy_db is very hard to locate.
    EXPECT_NE(std::string(ex.what()).find("InvalidType"), std::string::npos)
        << "exception message should name the unrecognised type: " << ex.what();
  }
}

TEST_F(StrategyBuilderTest, StrategyTypeMapCoversEveryBuilderCase) {
  // Regression guard for the omission that motivated the fix: DistrictMFT and
  // MFTAgeBased had builders and enum values but no entry here, so no
  // configuration could ever reach them.
  const std::vector<std::pair<std::string, IStrategy::StrategyType>> expected = {
      {"SFT", IStrategy::StrategyType::SFT},
      {"Cycling", IStrategy::StrategyType::Cycling},
      {"AdaptiveCycling", IStrategy::StrategyType::AdaptiveCycling},
      {"MFT", IStrategy::StrategyType::MFT},
      {"MFTRebalancing", IStrategy::StrategyType::MFTRebalancing},
      {"NestedMFT", IStrategy::StrategyType::NestedMFT},
      {"MFTMultiLocation", IStrategy::StrategyType::MFTMultiLocation},
      {"NestedMFTMultiLocation", IStrategy::StrategyType::NestedMFTMultiLocation},
      {"NovelDrugIntroduction", IStrategy::StrategyType::NovelDrugIntroduction},
      {"DistrictMFT", IStrategy::StrategyType::DistrictMft},
      {"MFTAgeBased", IStrategy::StrategyType::MFTAgeBased},
      {"PublicPrivate", IStrategy::StrategyType::PublicPrivate},
      {"PublicPrivateMultiLocation", IStrategy::StrategyType::PublicPrivateMultiLocation},
  };

  for (const auto &[name, type] : expected) {
    const auto entry = IStrategy::strategy_type_map.find(name);
    ASSERT_NE(entry, IStrategy::strategy_type_map.end())
        << "strategy_type_map is missing the type name: " << name;
    EXPECT_EQ(entry->second, type) << "wrong enum value registered for: " << name;
  }

  // If a new StrategyType is added, it needs a name here and a case in
  // StrategyBuilder::build, otherwise it is unreachable from configuration.
  EXPECT_EQ(IStrategy::strategy_type_map.size(), expected.size())
      << "strategy_type_map has entries not covered by this test";
}

TEST_F(StrategyBuilderTest, MissingTherapies) {
  YAML::Node node;
  node["name"] = "Missing Therapies";
  node["type"] = "SFT";
  
  // No therapy_ids specified
  EXPECT_THROW(StrategyBuilder::build(node, 6), std::runtime_error);
}

TEST_F(StrategyBuilderTest, BuildsAdditionalStrategyVariants) {
  auto rebalancing = create_mft_strategy_node();
  rebalancing["type"] = "MFTRebalancing";
  rebalancing["distribution"] = rebalancing["distributions"];
  rebalancing["update_duration_after_rebalancing"] = 30;
  rebalancing["delay_until_actual_trigger"] = 7;
  ASSERT_NE(StrategyBuilder::build(rebalancing, 5), nullptr);

  auto age_based = create_mft_strategy_node();
  age_based["type"] = "MFTAgeBased";
  age_based["age_boundaries"] = YAML::Load("[18]");
  ASSERT_NE(StrategyBuilder::build(age_based, 6), nullptr);

  auto nested = create_mft_strategy_node();
  nested["type"] = "NestedMFT";
  nested["start_distribution"] = YAML::Load("[0.5, 0.5]");
  nested["peak_distribution"] = YAML::Load("[0.2, 0.8]");
  nested["peak_after"] = 30;
  nested["strategy_ids"] = YAML::Load("[]");
  ASSERT_NE(StrategyBuilder::build(nested, 7), nullptr);

  auto multi_location = create_mft_strategy_node();
  multi_location["type"] = "MFTMultiLocation";
  multi_location["start_distribution_by_location"] = YAML::Load("[[0.5, 0.5]]");
  multi_location["peak_distribution_by_location"] = YAML::Load("[[0.2, 0.8]]");
  multi_location["peak_after"] = 30;
  ASSERT_NE(StrategyBuilder::build(multi_location, 8), nullptr);

  auto nested_multi_location = multi_location;
  nested_multi_location["type"] = "NestedMFTMultiLocation";
  nested_multi_location["strategy_ids"] = YAML::Load("[]");
  ASSERT_NE(StrategyBuilder::build(nested_multi_location, 9), nullptr);

  auto novel = YAML::Load(R"(
    name: novel
    type: NovelDrugIntroduction
    start_distribution: [0.5, 0.5]
    peak_distribution: [0.2, 0.8]
    peak_after: 30
    strategy_ids: []
    newly_introduced_strategy_id: 1
    tf_threshold: 0.2
    replacement_fraction: 0.5
    replacement_duration: 10
  )");
  ASSERT_NE(StrategyBuilder::build(novel, 10), nullptr);
}

TEST_F(StrategyBuilderTest, ValidatesPublicPrivateStrategyFields) {
  auto public_private = YAML::Load(R"(
    name: public-private
    type: PublicPrivate
    public_strategy_id: 0
    private_strategy_id: 1
    start_public_share: 0.3
    peak_public_share: 0.7
    peak_after: 30
  )");

  if (Model::get_strategy_db().size() > 1) {
    auto strategy = StrategyBuilder::build(public_private, 2);
    ASSERT_NE(strategy, nullptr);
    auto* result = dynamic_cast<PublicPrivateStrategy*>(strategy.get());
    ASSERT_NE(result, nullptr);
    EXPECT_DOUBLE_EQ(result->start_public_share, 0.3);
    EXPECT_DOUBLE_EQ(result->peak_public_share, 0.7);
  }

  public_private["start_public_share"] = 1.1;
  EXPECT_THROW(StrategyBuilder::build(public_private, 2), std::invalid_argument);
}

TEST_F(StrategyBuilderTest, ValidatesPublicPrivateMultiLocationFields) {
  YAML::Node node = YAML::Load(R"(
    name: public-private-multi
    type: PublicPrivateMultiLocation
    public_strategy_id: 0
    private_strategy_id: 1
    peak_after: 30
  )");
  YAML::Node start_shares;
  YAML::Node peak_shares;
  for (int location = 0; location < Model::get_config()->number_of_locations(); ++location) {
    start_shares.push_back(0.3);
    peak_shares.push_back(0.7);
  }
  node["start_public_share_by_location"] = start_shares;
  node["peak_public_share_by_location"] = peak_shares;

  auto strategy = StrategyBuilder::build(node, 2);
  ASSERT_NE(dynamic_cast<PublicPrivateMultiLocationStrategy*>(strategy.get()), nullptr);

  node["start_public_share_by_location"][0] = 1.1;
  EXPECT_THROW(StrategyBuilder::build(node, 2), std::invalid_argument);
  node["start_public_share_by_location"][0] = 0.3;
  node.remove("peak_after");
  EXPECT_THROW(StrategyBuilder::build(node, 2), std::invalid_argument);
  node["peak_after"] = 30;
  node["public_strategy_id"] = 2;
  node["private_strategy_id"] = 2;
  EXPECT_THROW(StrategyBuilder::build(node, 3), std::invalid_argument);

  node["public_strategy_id"] = 0;
  node["private_strategy_id"] = 1;
  node["peak_after"] = -1;
  EXPECT_THROW(StrategyBuilder::build(node, 2), std::invalid_argument);
  node["peak_after"] = 30;
  node["start_public_share_by_location"] = YAML::Load("[0.3, 0.3]");
  EXPECT_THROW(StrategyBuilder::build(node, 2), std::invalid_argument);
  node.remove("start_public_share_by_location");
  EXPECT_THROW(StrategyBuilder::build(node, 2), std::invalid_argument);

  node["start_public_share_by_location"] = start_shares;
  node.remove("private_strategy_id");
  EXPECT_THROW(StrategyBuilder::build(node, 2), std::invalid_argument);
  node["private_strategy_id"] = 9999;
  EXPECT_THROW(StrategyBuilder::build(node, 2), std::invalid_argument);
}

// ---------------------------------------------------------------------------
// DistrictMFT
//
// These exercise the DistrictMFT path through StrategyBuilder::build, which had
// no YAML-level coverage at all. DistrictMftStrategyTest constructs the object
// directly and so never touched the parsing, validation, or type dispatch that
// a real configuration goes through.
// ---------------------------------------------------------------------------

TEST_F(StrategyBuilderTest, FixtureProvidesThreeDistricts) {
  // The DistrictMFT tests below hard-code district ids 1..3. If the fixture
  // raster ever changes, fail here with an obvious reason rather than leaving
  // the coverage assertions to fail cryptically.
  const auto* boundary = Model::get_spatial_data()->get_boundary("district");
  ASSERT_NE(boundary, nullptr) << "the fixture must define a 'district' admin level";
  EXPECT_EQ(boundary->min_unit_id, 1);
  EXPECT_EQ(boundary->max_unit_id, k_district_count);
  EXPECT_EQ(boundary->unit_count, k_district_count);
}

TEST_F(StrategyBuilderTest, BuildDistrictMftStrategy) {
  YAML::Node node = create_full_coverage_district_mft_node();

  std::unique_ptr<IStrategy> strategy = StrategyBuilder::build(node, 11);

  // Before DistrictMFT was registered in strategy_type_map this returned an
  // SFTStrategy, or threw, depending on whether a stray therapy_ids key existed.
  ASSERT_NE(strategy, nullptr);
  EXPECT_EQ(strategy->get_type(), IStrategy::StrategyType::DistrictMft);
  EXPECT_NE(dynamic_cast<DistrictMftStrategy*>(strategy.get()), nullptr);
  EXPECT_EQ(strategy->id, 11);
  EXPECT_EQ(strategy->name, "Test District MFT");
}

TEST_F(StrategyBuilderTest, DistrictMftAcceptsDistributionThatSumsToOneInDecimal) {
  // Guards the truncation bug directly: static_cast<int>(sum) != 1 rejected
  // this distribution even though the values sum to 1.0 exactly in decimal.
  YAML::Node node = create_district_mft_node(R"(
    0:
      district_ids: [1, 2, 3]
      therapy_ids: [6, 8, 2, 0, 12, 14, 5]
      distribution: [0.765, 0.085, 0.0045, 0.063, 0.015, 0.006, 0.0615]
  )");

  EXPECT_NO_THROW({
    auto strategy = StrategyBuilder::build(node, 12);
    EXPECT_EQ(strategy->get_type(), IStrategy::StrategyType::DistrictMft);
  });
}

TEST_F(StrategyBuilderTest, DistrictMftRejectsDistributionSumMismatch) {
  // A distribution that is genuinely wrong must still be caught. The tolerance
  // added for floating-point accumulation is 1e-6, far below this error.
  YAML::Node node = create_district_mft_node(R"(
    0:
      district_ids: [1, 2, 3]
      therapy_ids: [6, 8]
      distribution: [0.5, 0.4]
  )");

  EXPECT_THROW(StrategyBuilder::build(node, 13), std::invalid_argument);
}

TEST_F(StrategyBuilderTest, DistrictMftRejectsIncompleteDistrictCoverage) {
  // District 3 is left unassigned. Partial coverage must be rejected at load
  // time, since get_therapy would otherwise hit a null MFT mid-run.
  YAML::Node node = create_district_mft_node(R"(
    0:
      district_ids: [1, 2]
      therapy_ids: [6, 8]
      distribution: [0.5, 0.5]
  )");

  EXPECT_THROW(StrategyBuilder::build(node, 14), std::invalid_argument);
}

TEST_F(StrategyBuilderTest, DistrictMftRejectsDuplicateDistrictAssignment) {
  YAML::Node node = create_district_mft_node(R"(
    0:
      district_ids: [1, 2, 3]
      therapy_ids: [6, 8]
      distribution: [0.5, 0.5]
    1:
      district_ids: [2]
      therapy_ids: [6, 8]
      distribution: [0.5, 0.5]
  )");

  EXPECT_THROW(StrategyBuilder::build(node, 15), std::invalid_argument);
}

TEST_F(StrategyBuilderTest, DistrictMftRejectsOutOfRangeDistrictId) {
  // The fixture has districts 1..3, so 4 is out of range.
  YAML::Node node = create_district_mft_node(R"(
    0:
      district_ids: [1, 2, 3, 4]
      therapy_ids: [6, 8]
      distribution: [0.5, 0.5]
  )");

  EXPECT_THROW(StrategyBuilder::build(node, 16), std::invalid_argument);
}

TEST_F(StrategyBuilderTest, DistrictMftRejectsTherapyIdEqualToDatabaseSize) {
  // Off-by-one guard. The bound check was `id > therapy_db.size()`, which let
  // an id exactly equal to the size through and then indexed out of bounds in
  // get_therapy.
  const auto therapy_count = static_cast<int>(Model::get_therapy_db().size());
  ASSERT_GT(therapy_count, 0);

  YAML::Node node = create_district_mft_node(
      "0:\n"
      "  district_ids: [1, 2, 3]\n"
      "  therapy_ids: [" + std::to_string(therapy_count) + ", 8]\n"
      "  distribution: [0.5, 0.5]\n");

  EXPECT_THROW(StrategyBuilder::build(node, 17), std::invalid_argument);
}

TEST_F(StrategyBuilderTest, DistrictMftRejectsMismatchedTherapyAndDistributionSizes) {
  YAML::Node node = create_district_mft_node(R"(
    0:
      district_ids: [1, 2, 3]
      therapy_ids: [6, 8, 2]
      distribution: [0.5, 0.5]
  )");

  EXPECT_THROW(StrategyBuilder::build(node, 18), std::invalid_argument);
}

TEST_F(StrategyBuilderTest, DistrictMftAcceptsZeroWeightedTherapy) {
  // A zero share is a legitimate way to keep an arm in the therapy list while
  // switching it off for a sweep, so it must not be rejected outright.
  YAML::Node node = create_district_mft_node(R"(
    0:
      district_ids: [1, 2, 3]
      therapy_ids: [6, 8, 2]
      distribution: [0.5, 0.5, 0.0]
  )");

  EXPECT_NO_THROW(StrategyBuilder::build(node, 19));
}

TEST_F(StrategyBuilderTest, DistrictMftStoresSharesAsDouble) {
  // Reverting percentages to float reintroduces a residual of roughly 1.5e-8
  // below 1.0, which get_therapy can fall into.
  static_assert(
      std::is_same_v<decltype(DistrictMftStrategy::MftStrategy::percentages)::value_type, double>,
      "DistrictMftStrategy::MftStrategy::percentages must hold double, not float");
  SUCCEED();
}
