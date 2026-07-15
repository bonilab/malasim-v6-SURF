#include <cmath>
#include <memory>

#include "Population/ImmuneSystem/ImmuneSystem.h"
#include "Population/ImmuneSystem/ImmuneSystemConstants.h"
#include "Population/Person/Person.h"
#include "Simulation/Model.h"
#include "Utils/Cli.h"
#include "fixtures/TestFileGenerators.h"
#include "gtest/gtest.h"

namespace {

class ImmuneSystemModelTest : public ::testing::Test {
protected:
  void SetUp() override {
    test_fixtures::setup_test_environment("test_input.yml");
    Model::get_instance()->release();
    utils::Cli::MaSimAppInput cli_input;
    cli_input.input_path = "test_input.yml";
    Model::set_cli_input(cli_input);
    ASSERT_TRUE(Model::get_instance()->initialize());

    person_ = std::make_unique<Person>();
    person_->initialize();
    person_->set_age(0);
    person_->set_latest_update_time(0);
    Model::get_scheduler()->set_current_time(0);
  }

  void TearDown() override {
    person_.reset();
    Model::get_instance()->release();
    test_fixtures::cleanup_test_files();
  }

  std::unique_ptr<Person> person_;
};

}  // namespace

TEST(ImmuneSystemTest, ConstructionAndSetters) {
  ImmuneSystem immune;
  EXPECT_EQ(immune.person(), nullptr);
  Person person;
  immune.set_person(&person);
  EXPECT_EQ(immune.person(), &person);
}

TEST(ImmuneSystemTest, UsesNonInfantModeByDefault) {
  ImmuneSystem immune;
  EXPECT_EQ(immune.mode(), ImmuneSystemMode::NON_INFANT);
}

TEST(ImmuneSystemTest, SetLatestImmuneValue) {
  ImmuneSystem immune;
  immune.set_latest_immune_value(0.5);
  EXPECT_DOUBLE_EQ(immune.get_latest_immune_value(), 0.5);
}

TEST(ImmuneSystemTest, IncreaseFlag) {
  ImmuneSystem immune;
  EXPECT_FALSE(immune.increase());
}

TEST(ImmuneSystemTest, SwitchToNonInfantPreservesValue) {
  ImmuneSystem immune;
  immune.initialize_as_infant();
  immune.set_latest_immune_value(0.25);

  immune.switch_to_non_infant();

  EXPECT_EQ(immune.mode(), ImmuneSystemMode::NON_INFANT);
  EXPECT_DOUBLE_EQ(immune.get_latest_immune_value(), 0.25);
}

TEST_F(ImmuneSystemModelTest, InfantModeUsesMaternalImmunityDecay) {
  auto* immune_system = person_->get_immune_system();
  immune_system->initialize_as_infant();
  immune_system->set_latest_immune_value(0.75);
  Model::get_scheduler()->set_current_time(10);

  const auto expected_value = 0.75 * std::exp(-immune::K_INFANT_IMMUNE_DECAY_RATE * 10);

  EXPECT_DOUBLE_EQ(immune_system->get_current_value(), expected_value);
}

TEST_F(ImmuneSystemModelTest, ModeSwitchIsContinuousAfterDailyUpdate) {
  auto* immune_system = person_->get_immune_system();
  immune_system->initialize_as_infant();
  immune_system->set_latest_immune_value(0.75);
  Model::get_scheduler()->set_current_time(10);

  const auto expected_value = 0.75 * std::exp(-immune::K_INFANT_IMMUNE_DECAY_RATE * 10);
  immune_system->update();
  person_->set_latest_update_time(Model::get_scheduler()->current_time());
  immune_system->switch_to_non_infant();

  EXPECT_EQ(immune_system->mode(), ImmuneSystemMode::NON_INFANT);
  EXPECT_DOUBLE_EQ(immune_system->get_latest_immune_value(), expected_value);
  EXPECT_DOUBLE_EQ(immune_system->get_current_value(), expected_value);
}
