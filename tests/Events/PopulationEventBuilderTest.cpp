#include <gtest/gtest.h>

#include "Configuration/Config.h"
#include "Events/Population/AnnualBetaUpdateEvent.hxx"
#include "Events/Population/AnnualCoverageUpdateEvent.hxx"
#include "Events/Population/ChangeCirculationPercentEvent.hxx"
#include "Events/Population/ChangeInterruptedFeedingRateEvent.h"
#include "Events/Population/ChangeMutationMaskEvent.h"
#include "Events/Population/ChangeMutationProbabilityPerLocusEvent.h"
#include "Events/Population/ChangeWithinHostInducedFreeRecombinationEvent.h"
#include "Events/Population/PopulationEventBuilder.h"
#include "Events/Population/SingleRoundMDAEvent.h"
#include "Events/Population/TurnOffMutationEvent.h"
#include "Events/Population/TurnOnMutationEvent.h"
#include "Simulation/Model.h"
#include "Utils/Cli.h"
#include "fixtures/TestFileGenerators.h"

class PopulationEventBuilderTest : public ::testing::Test {
protected:
  void SetUp() override {
    test_fixtures::setup_test_environment();
    utils::Cli::MaSimAppInput cli_input;
    cli_input.input_path = "test_input.yml";
    Model::set_cli_input(cli_input);
    ASSERT_TRUE(Model::get_instance()->initialize());
  }

  void TearDown() override {
    Model::get_instance()->release();
    test_fixtures::cleanup_test_files();
  }

  Config* config() { return Model::get_config(); }
};

TEST_F(PopulationEventBuilderTest, BuildsMutationAndMdaEventsWithDatesAndValues) {
  const auto mutation = PopulationEventBuilder::build_turn_on_mutation_event(
      YAML::Load("- date: 2024/01/11\n  mutation_probability: 0.25"), config());
  ASSERT_EQ(mutation.size(), 1);
  auto* turn_on = dynamic_cast<TurnOnMutationEvent*>(mutation[0].get());
  ASSERT_NE(turn_on, nullptr);
  EXPECT_EQ(turn_on->get_time(),
            (date::sys_days{date::year{2024} / 1 / 11}
             - date::sys_days{config()->get_simulation_timeframe().get_starting_date()})
                .count());
  EXPECT_DOUBLE_EQ(turn_on->mutation_probability, 0.25);

  const auto turn_off = PopulationEventBuilder::build_turn_off_mutation_event(
      YAML::Load("- date: 2024/01/12"), config());
  ASSERT_EQ(turn_off.size(), 1);
  EXPECT_EQ(turn_off[0]->get_time(),
            (date::sys_days{date::year{2024} / 1 / 12}
             - date::sys_days{config()->get_simulation_timeframe().get_starting_date()})
                .count());

  const auto mda = PopulationEventBuilder::build_single_round_mda_event(
      YAML::Load("- date: 2024/01/02\n  fraction_population_targeted: [0.2]\n  days_to_complete_all_treatments: 5"),
      config());
  ASSERT_EQ(mda.size(), 1);
  auto* mda_event = dynamic_cast<SingleRoundMDAEvent*>(mda[0].get());
  ASSERT_NE(mda_event, nullptr);
  EXPECT_EQ(mda_event->get_fraction_population_targeted().size(), config()->number_of_locations());
  EXPECT_DOUBLE_EQ(mda_event->get_fraction_population_targeted()[0], 0.2);
  EXPECT_EQ(mda_event->get_days_to_complete(), 5);
}

TEST_F(PopulationEventBuilderTest, BuildsConfigurationChangeEvents) {
  const auto recombination = PopulationEventBuilder::build_change_within_host_induced_free_recombination_events(
      YAML::Load("- date: 2024/01/03\n  value: true"), config());
  ASSERT_EQ(recombination.size(), 1);
  auto* recombination_event =
      dynamic_cast<ChangeWithinHostInducedFreeRecombinationEvent*>(recombination[0].get());
  ASSERT_NE(recombination_event, nullptr);
  EXPECT_TRUE(recombination_event->value);
  EXPECT_EQ(recombination_event->get_time(),
            (date::sys_days{date::year{2024} / 1 / 3}
             - date::sys_days{config()->get_simulation_timeframe().get_starting_date()})
                .count());

  const auto probability = PopulationEventBuilder::build_change_mutation_probability_per_locus_events(
      YAML::Load("- date: 2024/01/04\n  mutation_probability_per_locus: 0.5"), config());
  ASSERT_EQ(probability.size(), 1);
  EXPECT_DOUBLE_EQ(dynamic_cast<ChangeMutationProbabilityPerLocusEvent*>(probability[0].get())->value,
                   0.5);

  const auto mask = PopulationEventBuilder::build_change_mutation_mask_events(
      YAML::Load("- date: 2024/01/05\n  mutation_mask: [true, false, true]"), config());
  ASSERT_EQ(mask.size(), 1);
  EXPECT_EQ(dynamic_cast<ChangeMutationMaskEvent*>(mask[0].get())->mask,
            std::vector<bool>({true, false, true}));

  const auto ifr = PopulationEventBuilder::build_change_interrupted_feeding_rate_event(
      YAML::Load("- location: 0\n  date: 2024/01/06\n  interrupted_feeding_rate: 0.3"), config());
  ASSERT_EQ(ifr.size(), 1);
  auto* ifr_event = dynamic_cast<ChangeInterruptedFeedingRateEvent*>(ifr[0].get());
  ASSERT_NE(ifr_event, nullptr);
  EXPECT_EQ(ifr_event->location, 0);
  EXPECT_DOUBLE_EQ(ifr_event->ifr, 0.3);
}

TEST_F(PopulationEventBuilderTest, RejectsInvalidAnnualAndCirculationDefinitions) {
  EXPECT_THROW(PopulationEventBuilder::build_annual_beta_update_event(
                   YAML::Load("[]"), config()),
               std::invalid_argument);
  EXPECT_THROW(PopulationEventBuilder::build_annual_coverage_update_event(
                   YAML::Load("- date: 2024/01/01\n  rate: 0.1\n- date: 2024/02/01\n  rate: 0.1"), config()),
               std::invalid_argument);
  EXPECT_THROW(PopulationEventBuilder::build_change_circulation_percent_event(
                   YAML::Load("- date: 2024/01/01\n  circulation_percent: 1.1"), config()),
               std::invalid_argument);
  EXPECT_THROW(PopulationEventBuilder::build_change_circulation_percent_event(
                   YAML::Load("- date: 2024/01/01\n  circulation_percent: -0.1"), config()),
               std::invalid_argument);
}

TEST_F(PopulationEventBuilderTest, BuildsThroughPublicDispatcher) {
  YAML::Node node = YAML::Load(
      "name: change_mutation_probability_per_locus\n"
      "info:\n"
      "  - date: 2024/01/08\n"
      "    mutation_probability_per_locus: 0.2");
  const auto events = PopulationEventBuilder::build(node);
  ASSERT_EQ(events.size(), 1);
  EXPECT_EQ(events[0]->name(), ChangeMutationProbabilityPerLocusEvent::EVENT_NAME);
  EXPECT_EQ(events[0]->get_time(),
            (date::sys_days{date::year{2024} / 1 / 8}
             - date::sys_days{config()->get_simulation_timeframe().get_starting_date()})
                .count());
}
