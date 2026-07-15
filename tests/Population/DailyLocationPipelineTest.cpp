#include <gtest/gtest.h>

#include <algorithm>

#include "Mosquito/Mosquito.h"
#include "Population/Population.h"
#include "Simulation/Model.h"
#include "Utils/Cli.h"
#include "Utils/Constants.h"
#include "fixtures/TestFileGenerators.h"

class DailyLocationPipelineTest : public ::testing::Test {
protected:
  void SetUp() override {
    test_fixtures::create_complete_test_environment();
    utils::Cli::MaSimAppInput cli_input;
    cli_input.input_path = "test_input.yml";
    Model::set_cli_input(cli_input);
    ASSERT_TRUE(Model::get_instance()->initialize());
  }

  void TearDown() override { test_fixtures::cleanup_test_files(); }
};

TEST_F(DailyLocationPipelineTest, NewbornIsIncludedInAlignedSamplingState) {
  auto* population = Model::get_population();
  constexpr int location = 0;
  auto* newborn = population->give_1_birth(location);

  population->update_current_foi();

  const auto &alive_people = population->all_alive_persons_by_location()[location];
  EXPECT_NE(std::find(alive_people.begin(), alive_people.end(), newborn), alive_people.end());
  EXPECT_EQ(alive_people.size(), population->individual_foi_by_location()[location].size());
  EXPECT_EQ(alive_people.size(),
            population->individual_relative_biting_by_location()[location].size());
  EXPECT_EQ(alive_people.size(),
            population->individual_relative_moving_by_location()[location].size());
}

TEST_F(DailyLocationPipelineTest, LocationBirthEventAppendsNewbornSamplingEntries) {
  auto* population = Model::get_population();
  constexpr int location = 0;
  population->update_current_foi();
  const auto initial_size = population->all_alive_persons_by_location()[location].size();
  auto demographic = Model::get_config()->get_population_demographic();
  demographic.set_birth_rate(Constants::DAYS_IN_YEAR);
  Model::get_config()->set_population_demographic(demographic);

  population->perform_birth_event_at_location(location);

  const auto &alive_people = population->all_alive_persons_by_location()[location];
  EXPECT_GT(alive_people.size(), initial_size);
  EXPECT_EQ(alive_people.size(), population->individual_foi_by_location()[location].size());
  EXPECT_EQ(alive_people.size(),
            population->individual_relative_biting_by_location()[location].size());
  EXPECT_EQ(alive_people.size(),
            population->individual_relative_moving_by_location()[location].size());
}

TEST_F(DailyLocationPipelineTest, DeadPersonIsRemovedFromAlignedSamplingStateBeforeDeletion) {
  auto* population = Model::get_population();
  constexpr int location = 0;
  auto* person = population->give_1_birth(location);
  population->update_current_foi();
  const auto expected_size = population->all_alive_persons_by_location()[location].size() - 1;

  person->set_host_state(Person::DEAD);
  population->clear_dead_people_at_location(location);

  const auto &alive_people = population->all_alive_persons_by_location()[location];
  EXPECT_EQ(alive_people.size(), expected_size);
  EXPECT_EQ(alive_people.size(), population->individual_foi_by_location()[location].size());
  EXPECT_EQ(alive_people.size(),
            population->individual_relative_biting_by_location()[location].size());
  EXPECT_EQ(alive_people.size(),
            population->individual_relative_moving_by_location()[location].size());
}

TEST_F(DailyLocationPipelineTest, ZeroFoiLocationDoesNotPreventLaterCohortClearing) {
  auto* config = Model::get_config();
  ASSERT_GE(config->number_of_locations(), 2);
  auto* population = Model::get_population();
  auto* mosquito = Model::get_mosquito();
  auto* genotype = Model::get_genotype_db()->at(0);
  constexpr int tracking_index = 0;

  population->current_force_of_infection_by_location()[0] = 0.0;
  population->current_force_of_infection_by_location()[1] = 0.0;
  mosquito->genotypes_table[tracking_index][0][0] = genotype;
  mosquito->genotypes_table[tracking_index][1][0] = genotype;

  mosquito->infect_new_cohort_in_prmc(config, Model::get_random(), population, tracking_index);

  EXPECT_EQ(mosquito->genotypes_table[tracking_index][0][0], nullptr);
  EXPECT_EQ(mosquito->genotypes_table[tracking_index][1][0], nullptr);
}

TEST_F(DailyLocationPipelineTest, DailyUpdateLeavesEveryLocationSamplingStateAligned) {
  EXPECT_NO_THROW(Model::get_instance()->daily_update());

  auto* population = Model::get_population();
  for (auto location = 0; location < Model::get_config()->number_of_locations(); ++location) {
    const auto &alive_people = population->all_alive_persons_by_location()[location];
    EXPECT_EQ(alive_people.size(), population->individual_foi_by_location()[location].size());
    EXPECT_EQ(alive_people.size(),
              population->individual_relative_biting_by_location()[location].size());
    EXPECT_EQ(alive_people.size(),
              population->individual_relative_moving_by_location()[location].size());
    for (const auto* person : alive_people) { EXPECT_NE(person->get_host_state(), Person::DEAD); }
  }
}
