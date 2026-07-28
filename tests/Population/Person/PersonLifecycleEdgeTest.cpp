#include <gtest/gtest.h>

#include "Events/BirthdayEvent.h"
#include "Population/Person/Person.h"
#include "Treatment/Therapies/SCTherapy.h"
#include "PersonTestBase.h"

class PersonLifecycleEdgeTest : public PersonTestBase {
protected:
  void SetUp() override {
    PersonTestBase::SetUp();
    mock_scheduler_->set_current_time(5);
  }
};

TEST_F(PersonLifecycleEdgeTest, ComputesPartialTherapyComplianceAcrossProbabilityBranches) {
  auto epidemiology = mock_config_->get_epidemiological_parameters();
  epidemiology.set_p_compliance(0.5);
  epidemiology.set_min_dosing_days(3);
  mock_config_->set_epidemiological_parameters(epidemiology);

  EXPECT_CALL(*mock_random_, random_flat(_, _)).WillOnce(Return(0.9));
  EXPECT_EQ(Person::complied_dosing_days(1), 3);
  EXPECT_CALL(*mock_random_, random_flat(_, _)).WillOnce(Return(0.1));
  EXPECT_EQ(Person::complied_dosing_days(1), 1);

  epidemiology.set_p_compliance(1.0);
  mock_config_->set_epidemiological_parameters(epidemiology);
  EXPECT_EQ(Person::complied_dosing_days(2), 2);
}

TEST_F(PersonLifecycleEdgeTest, RejectsInvalidBirthdayDelayAndSchedulesLifecycleEvents) {
  EXPECT_THROW(person_->schedule_birthday_event(0), std::invalid_argument);
  EXPECT_NO_THROW(person_->schedule_birthday_event(2));
  EXPECT_TRUE(person_->has_birthday_event());

  person_->schedule_rapt_event(2);
  person_->schedule_receive_mda_therapy_event(nullptr, 2);
  person_->schedule_receive_therapy_event(nullptr, nullptr, 2);
  person_->schedule_switch_immune_system_mode_event(2);
  person_->schedule_move_parasite_to_blood(nullptr, 2);
  person_->schedule_mature_gametocyte_event(nullptr);
  person_->schedule_update_by_drug_event(nullptr);
  person_->schedule_move_to_target_location_next_day_event(1);
  person_->schedule_return_to_residence_event(2);

  EXPECT_TRUE(person_->has_update_by_having_drug_event());
  EXPECT_GE(person_->get_events().size(), 8u);
}

TEST_F(PersonLifecycleEdgeTest, ComplianceOneUsesTheFullDosingDayCount) {
  auto epidemiology = mock_config_->get_epidemiological_parameters();
  epidemiology.set_p_compliance(1.0);
  mock_config_->set_epidemiological_parameters(epidemiology);
  EXPECT_CALL(*mock_random_, random_flat(_, _)).Times(0);
  EXPECT_EQ(Person::complied_dosing_days(7), 7);
}

TEST_F(PersonLifecycleEdgeTest, ComputesVariableSingleCompoundTherapyCompletionDay) {
  SCTherapy therapy;
  therapy.dosing_day = {1, 3};
  therapy.calculate_max_dosing_day();
  therapy.set_full_compliance(true);
  EXPECT_EQ(Person::complied_dosing_days(&therapy), 3);

  therapy.set_full_compliance(false);
  therapy.pr_completed_days = {0.25, 0.75};
  EXPECT_CALL(*mock_random_, random_flat(_, _)).WillOnce(Return(0.1));
  EXPECT_EQ(Person::complied_dosing_days(&therapy), 1);
  EXPECT_CALL(*mock_random_, random_flat(_, _)).WillOnce(Return(0.5));
  EXPECT_EQ(Person::complied_dosing_days(&therapy), 2);
}

TEST_F(PersonLifecycleEdgeTest, RejectsPastEvents) {
  auto event = std::make_unique<BirthdayEvent>(person_.get());
  event->set_time(0);
  EXPECT_THROW(person_->schedule_basic_event(std::move(event)), std::invalid_argument);
}
