#include <gtest/gtest.h>
#include <yaml-cpp/yaml.h>

#include <functional>
#include <stdexcept>
#include <string>

#include "Configuration/Config.h"
#include "Simulation/Model.h"
#include "Utils/Cli.h"
#include "fixtures/TestFileGenerators.h"

namespace {

constexpr char TEST_CONFIG_FILE[] = "test_input.yml";

YAML::Node make_candidate(const double p_ci_symp, const double z, const double kappa,
                          const double midpoint, const double p_seek_base,
                          const double mutation_probability_per_locus) {
  YAML::Node candidate;
  candidate["p_ci_symp"] = p_ci_symp;
  candidate["z"] = z;
  candidate["kappa"] = kappa;
  candidate["midpoint"] = midpoint;
  candidate["p_seek_base"] = p_seek_base;
  candidate["mutation_probability_per_locus"] = mutation_probability_per_locus;
  return candidate;
}

void add_candidate_section(YAML::Node &config, const bool random_selection) {
  auto candidates = config["immune_system_paprameter_candidates"];
  candidates["used_in_simulation"] = 2;
  candidates["random_selection"] = random_selection;
  candidates["candidates"][2] = make_candidate(0.25, 2.0, 0.2, 0.3, 0.65, 0.002);
  candidates["candidates"][7] = make_candidate(0.35, 3.0, 0.3, 0.4, 0.75, 0.003);
  candidates["candidates"][15] = make_candidate(0.45, 4.0, 0.4, 0.5, 0.85, 0.004);
}

class ConfigTest : public ::testing::Test {
protected:
  void SetUp() override { Model::get_instance()->release(); }

  void TearDown() override {
    Model::get_instance()->release();
    test_fixtures::cleanup_test_files();
  }

  static void write_config(const std::function<void(YAML::Node &)> &modify = nullptr) {
    test_fixtures::setup_test_environment(TEST_CONFIG_FILE, modify);
  }

  static bool initialize_model() {
    utils::Cli::MaSimAppInput input;
    input.input_path = TEST_CONFIG_FILE;
    Model::set_cli_input(input);
    return Model::get_instance()->initialize();
  }
};

TEST_F(ConfigTest, LoadsCompleteConfiguration) {
  write_config();

  ASSERT_TRUE(initialize_model());

  auto* actual_config = Model::get_config();
  ASSERT_NE(actual_config, nullptr);
  EXPECT_GT(actual_config->number_of_locations(), 0U);
  EXPECT_EQ(actual_config->number_of_locations(), actual_config->location_db().size());
  EXPECT_GT(actual_config->number_of_age_classes(), 0);
  EXPECT_GT(actual_config->number_of_tracking_days(), 0);
  EXPECT_FALSE(actual_config->get_population_events().get_events_raw().empty());
}

TEST_F(ConfigTest, MissingFileReturnsFalse) {
  Config config;
  EXPECT_FALSE(config.load("missing_config.yml"));
}

TEST_F(ConfigTest, RejectsPopulationEventOutsideSimulationWindow) {
  write_config(
      [](YAML::Node &config) { config["population_events"][0]["info"][0]["date"] = "1999/12/31"; });

  EXPECT_THROW(initialize_model(), std::invalid_argument);
}

TEST_F(ConfigTest, AppliesImmuneSystemCandidateOverrides) {
  write_config([](YAML::Node &config) { add_candidate_section(config, false); });

  ASSERT_TRUE(initialize_model());

  auto* actual_config = Model::get_config();
  ASSERT_NE(actual_config, nullptr);
  EXPECT_TRUE(actual_config->has_immune_system_parameter_candidates());
  EXPECT_EQ(actual_config->get_immune_system_parameter_candidates().get_used_in_simulation(), 2);
  EXPECT_DOUBLE_EQ(
      actual_config->get_immune_system_parameters().get_immune_effect_on_progression_to_clinical(),
      2.0);
  EXPECT_DOUBLE_EQ(
      actual_config->get_immune_system_parameters().get_factor_effect_age_mature_immunity(), 0.2);
  EXPECT_DOUBLE_EQ(actual_config->get_immune_system_parameters().get_midpoint(), 0.3);
  EXPECT_DOUBLE_EQ(actual_config->get_epidemiological_parameters()
                       .get_allow_new_coinfection_to_cause_symptoms()
                       .get_probability(),
                   0.25);
  EXPECT_DOUBLE_EQ(actual_config->get_epidemiological_parameters()
                       .get_age_based_probability_of_seeking_treatment()
                       .get_power()
                       .base,
                   0.65);
  EXPECT_DOUBLE_EQ(actual_config->get_genotype_parameters().get_mutation_probability_per_locus(),
                   0.002);
}

TEST_F(ConfigTest, RandomCandidateSelectionIsReproducibleWithConfiguredSeed) {
  write_config([](YAML::Node &config) {
    config["model_settings"]["initial_seed_number"] = 12345;
    add_candidate_section(config, true);
  });

  ASSERT_TRUE(initialize_model());
  const int expected_candidate =
      Model::get_config()->get_immune_system_parameter_candidates().get_used_in_simulation();
  EXPECT_TRUE(expected_candidate == 2 || expected_candidate == 7 || expected_candidate == 15);

  Model::get_instance()->release();
  ASSERT_TRUE(initialize_model());
  const int actual_candidate =
      Model::get_config()->get_immune_system_parameter_candidates().get_used_in_simulation();
  EXPECT_EQ(actual_candidate, expected_candidate);
}

TEST_F(ConfigTest, ReloadClearsOmittedOptionalSections) {
  write_config([](YAML::Node &config) {
    add_candidate_section(config, false);
    config["rapt_settings"]["enabled"] = true;
    config["rapt_settings"]["period"] = 24;
  });
  ASSERT_TRUE(initialize_model());
  ASSERT_TRUE(Model::get_config()->has_immune_system_parameter_candidates());
  ASSERT_EQ(Model::get_config()->get_rapt_settings().get_period(), 24);

  write_config([](YAML::Node &config) {
    config.remove("immune_system_paprameter_candidates");
    config.remove("rapt_settings");
  });

  Model::get_config()->reload();

  EXPECT_FALSE(Model::get_config()->has_immune_system_parameter_candidates());
  EXPECT_TRUE(
      Model::get_config()->get_immune_system_parameter_candidates().get_candidates().empty());
  EXPECT_EQ(Model::get_config()->get_rapt_settings().get_period(), 12);
}

}  // namespace
