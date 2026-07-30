#include "Events/EndClinicalEvent.h"
#include "Events/ProgressToClinicalEvent.h"
#include "Events/ReportTreatmentFailureDeathEvent.h"
#include "Events/TestTreatmentFailureEvent.h"
#include "PersonTestBase.h"
#include "gmock/gmock.h"

using ::testing::Eq;

class PersonClinicalTest : public PersonTestBase {
protected:
  static constexpr core::Age K_TEST_AGE = 0;
  static constexpr double K_BASE_MORTALITY = 0.4;

  void SetUp() override {
    PersonTestBase::SetUp();

    // Establish valid, deterministic demographic state without making fixture
    // initialization part of each test's notification expectations.
    person_->set_population(nullptr);
    person_->set_age(K_TEST_AGE);
    person_->set_population(mock_population_);

    const auto age_class_count = mock_config_->number_of_age_classes();
    mock_config_->get_population_demographic().set_mortality_when_treatment_fail_by_age_class(
        std::vector<double>(age_class_count, K_BASE_MORTALITY));

    mock_scheduler_->set_current_time(5);
  }
};

TEST_F(PersonClinicalTest, ScheduleEndClinicalEvent) {
  EXPECT_CALL(*mock_random_, random_normal_int(_, _)).WillOnce(Return(7));

  auto test_parasite = std::make_unique<ClonalParasitePopulation>();
  person_->schedule_end_clinical_event(test_parasite.get());

  EXPECT_TRUE(person_->has_event<EndClinicalEvent>());
  EXPECT_EQ(person_->get_events().size(), 1);

  // cast to EndClinicalEvent
  auto event = dynamic_cast<EndClinicalEvent*>(person_->get_events().begin()->second.get());
  EXPECT_EQ(event->get_time(), 12);
  EXPECT_EQ(event->clinical_caused_parasite(), test_parasite.get());
}

TEST_F(PersonClinicalTest, ScheduleProgressToClinicalEventUnderFive) {
  // Age 3 remains in the age class established by the fixture.
  EXPECT_CALL(*mock_population_, notify_change(Eq(person_.get()), Person::Property::AGE, _, _));
  person_->set_age(3);

  auto test_parasite = std::make_unique<ClonalParasitePopulation>();
  person_->schedule_progress_to_clinical_event(test_parasite.get());
  EXPECT_TRUE(person_->has_event<ProgressToClinicalEvent>());
  EXPECT_EQ(person_->get_events().size(), 1);

  // cast to ProgressToClinicalEvent
  auto event = dynamic_cast<ProgressToClinicalEvent*>(person_->get_events().begin()->second.get());
  EXPECT_EQ(event->get_time(), 10);
  EXPECT_EQ(event->clinical_caused_parasite(), test_parasite.get());
}

TEST_F(PersonClinicalTest, ScheduleProgressToClinicalEventOverFive) {
  // 1 for initial age change,
  EXPECT_CALL(*mock_population_, notify_change(_, _, _, _)).Times(2);

  person_->set_age(6);
  auto test_parasite = std::make_unique<ClonalParasitePopulation>();
  person_->schedule_progress_to_clinical_event(test_parasite.get());
  EXPECT_TRUE(person_->has_event<ProgressToClinicalEvent>());
  EXPECT_EQ(person_->get_events().size(), 1);

  // cast to ProgressToClinicalEvent
  auto event = dynamic_cast<ProgressToClinicalEvent*>(person_->get_events().begin()->second.get());
  EXPECT_EQ(event->get_time(), 12);
  EXPECT_EQ(event->clinical_caused_parasite(), test_parasite.get());
}

TEST_F(PersonClinicalTest, CancelClinicalEvents) {
  EXPECT_CALL(*mock_population_, notify_change(_, _, _, _)).Times(2);

  person_->set_age(25);

  // Schedule two events and keep track of the first one
  auto test_parasite1 = std::make_unique<ClonalParasitePopulation>();
  auto test_parasite2 = std::make_unique<ClonalParasitePopulation>();
  person_->schedule_progress_to_clinical_event(test_parasite1.get());
  person_->schedule_progress_to_clinical_event(test_parasite2.get());
  auto test_event = person_->get_events().begin()->second.get();

  // Cancel all except the first event
  person_->cancel_all_other_progress_to_clinical_events_except(test_event);

  // Count executable events
  int executable_event_count = 0;
  for (const auto &[_, event] : person_->get_events()) {
    if (auto* progress_event = dynamic_cast<ProgressToClinicalEvent*>(event.get())) {
      if (event->is_executable()) { executable_event_count++; }
    }
  }
  EXPECT_EQ(executable_event_count, 1);
}

TEST_F(PersonClinicalTest, ProbabilityProgressToClinical) {
  const double expected_probability = 0.7;
  EXPECT_CALL(*mock_immune_system_, get_clinical_progression_probability())
      .WillOnce(Return(expected_probability));

  EXPECT_DOUBLE_EQ(person_->get_probability_progress_to_clinical(), expected_probability);
}
TEST_F(PersonClinicalTest, DeathProgressionWithoutTreatment) {
  const double input_at_death_threshold = K_BASE_MORTALITY;
  const double input_above_death_threshold = K_BASE_MORTALITY + 0.1;

  EXPECT_CALL(*mock_random_, random_flat(0.0, 1.0)).WillOnce(Return(input_at_death_threshold));
  EXPECT_TRUE(person_->will_progress_to_death_when_receive_no_treatment());

  EXPECT_CALL(*mock_random_, random_flat(0.0, 1.0)).WillOnce(Return(input_above_death_threshold));
  EXPECT_FALSE(person_->will_progress_to_death_when_receive_no_treatment());
}

TEST_F(PersonClinicalTest, DeathProgressionWithTreatment) {
  const double expected_treated_mortality = K_BASE_MORTALITY * 0.1;
  const double input_at_death_threshold = expected_treated_mortality;
  const double input_above_death_threshold = expected_treated_mortality + 0.01;

  EXPECT_CALL(*mock_random_, random_flat(0.0, 1.0)).WillOnce(Return(input_at_death_threshold));
  EXPECT_TRUE(person_->will_progress_to_death_when_recieve_treatment());

  EXPECT_CALL(*mock_random_, random_flat(0.0, 1.0)).WillOnce(Return(input_above_death_threshold));
  EXPECT_FALSE(person_->will_progress_to_death_when_recieve_treatment());
}

TEST_F(PersonClinicalTest, TestTreatmentFailureScheduling) {
  const int testing_day = 28;
  const int therapy_id = 1;
  auto test_parasite = std::make_unique<ClonalParasitePopulation>();

  person_->schedule_test_treatment_failure_event(test_parasite.get(), testing_day, therapy_id);
  EXPECT_FALSE(person_->get_events().empty());

  // cast to TestTreatmentFailureEvent
  auto event =
      dynamic_cast<TestTreatmentFailureEvent*>(person_->get_events().begin()->second.get());
  // current time is 5
  EXPECT_EQ(event->get_time(), testing_day + 5);
  EXPECT_EQ(event->clinical_caused_parasite(), test_parasite.get());
  EXPECT_EQ(event->therapy_id(), therapy_id);
}

TEST_F(PersonClinicalTest, ReportTreatmentFailureDeathScheduling) {
  EXPECT_CALL(*mock_population_, notify_change(_, Person::Property::AGE_CLASS, _, _)).Times(1);
  EXPECT_CALL(*mock_population_, notify_change(_, Person::Property::LOCATION, _, _)).Times(1);

  const int therapy_id = 1;
  const int testing_day = 30;

  person_->set_age_class(1);
  person_->set_location(1);

  person_->schedule_report_treatment_failure_death_event(therapy_id, testing_day);
  EXPECT_FALSE(person_->get_events().empty());

  // cast to ReportTreatmentFailureDeathEvent
  auto event =
      dynamic_cast<ReportTreatmentFailureDeathEvent*>(person_->get_events().begin()->second.get());
  // current time is 5
  EXPECT_EQ(event->get_time(), testing_day + 5);
  EXPECT_EQ(event->therapy_id(), therapy_id);
  EXPECT_EQ(event->age_class(), person_->get_age_class());
  EXPECT_EQ(event->location_id(), person_->get_location());
}
