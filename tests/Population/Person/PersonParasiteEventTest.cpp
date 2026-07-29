#include "Events/MatureGametocyteEvent.h"
#include "Events/MoveParasiteToBloodEvent.h"
#include "Events/UpdateWhenDrugIsPresentEvent.h"
#include "Parasites/Genotype.h"
#include "PersonParasiteTestBase.h"

TEST_F(PersonParasiteTest, ScheduleMoveParasiteToBloodEvent) {
  auto genotype_ptr = std::make_unique<Genotype>("aaabbbccc");
  genotype_ptr->set_genotype_id(1);
  auto genotype_raw = genotype_ptr.get();
  Model::get_genotype_db()->add(std::move(genotype_ptr));

  person_->schedule_move_parasite_to_blood(genotype_raw, 5);

  ASSERT_EQ(person_->get_events().size(), 1);
  auto event = dynamic_cast<MoveParasiteToBloodEvent*>(person_->get_events().begin()->second.get());
  ASSERT_NE(event, nullptr);
  EXPECT_EQ(event->infection_genotype(), genotype_raw);
  EXPECT_EQ(event->get_time(), current_time_ + 5);
}

TEST_F(PersonParasiteTest, ScheduleMatureGametocyteEventOver5) {
  EXPECT_CALL(*mock_population_, notify_change(_, Person::Property::AGE, _, _)).Times(1);
  EXPECT_CALL(*mock_population_, notify_change(_, Person::Property::AGE_CLASS, _, _)).Times(1);

  auto genotype_ptr = std::make_unique<Genotype>("aaabbbccc");
  genotype_ptr->set_genotype_id(1);
  auto genotype_raw = genotype_ptr.get();
  Model::get_genotype_db()->add(std::move(genotype_ptr));
  auto parasite = person_->add_new_parasite_to_blood(genotype_raw);
  person_->set_age(10);

  person_->schedule_mature_gametocyte_event(parasite);

  ASSERT_EQ(person_->get_events().size(), 1);
  auto event = dynamic_cast<MatureGametocyteEvent*>(person_->get_events().begin()->second.get());
  ASSERT_NE(event, nullptr);
  EXPECT_EQ(event->blood_parasite(), parasite);
  EXPECT_EQ(event->get_time(), current_time_
                                   + Model::get_config()
                                         ->get_epidemiological_parameters()
                                         .get_days_mature_gametocyte_over_five());
}

TEST_F(PersonParasiteTest, ScheduleMatureGametocyteEventUnder5) {
  EXPECT_CALL(*mock_population_, notify_change(_, Person::Property::AGE, _, _));
  EXPECT_CALL(*mock_population_, notify_change(_, Person::Property::AGE_CLASS, _, _));

  auto genotype_ptr = std::make_unique<Genotype>("aaabbbccc");
  genotype_ptr->set_genotype_id(1);
  auto genotype_raw = genotype_ptr.get();
  Model::get_genotype_db()->add(std::move(genotype_ptr));
  auto parasite = person_->add_new_parasite_to_blood(genotype_raw);
  person_->set_age(4);

  person_->schedule_mature_gametocyte_event(parasite);

  ASSERT_EQ(person_->get_events().size(), 1);
  auto event = dynamic_cast<MatureGametocyteEvent*>(person_->get_events().begin()->second.get());
  ASSERT_NE(event, nullptr);
  EXPECT_EQ(event->blood_parasite(), parasite);
  EXPECT_EQ(event->get_time(), current_time_
                                   + Model::get_config()
                                         ->get_epidemiological_parameters()
                                         .get_days_mature_gametocyte_under_five());
}

TEST_F(PersonParasiteTest, ScheduleUpdateByDrugEvent) {
  auto genotype_ptr = std::make_unique<Genotype>("aaabbbccc");
  genotype_ptr->set_genotype_id(1);
  auto genotype_raw = genotype_ptr.get();
  Model::get_genotype_db()->add(std::move(genotype_ptr));
  auto parasite = person_->add_new_parasite_to_blood(genotype_raw);

  person_->schedule_update_by_drug_event(parasite);

  ASSERT_EQ(person_->get_events().size(), 1);
  auto event =
      dynamic_cast<UpdateWhenDrugIsPresentEvent*>(person_->get_events().begin()->second.get());
  ASSERT_NE(event, nullptr);
  EXPECT_EQ(event->clinical_caused_parasite(), parasite);
  EXPECT_EQ(event->get_time(), current_time_ + 1);
}

TEST_F(PersonParasiteTest, InfectionFromInfectiousBite) {
  EXPECT_CALL(*mock_immune_system_, get_current_value()).WillOnce(::testing::Return(0.5));

  const double prob = person_->p_infection_from_an_infectious_bite();
  const auto expected_prob = (1 - 0.5) / 8.333 + 0.04;
  EXPECT_DOUBLE_EQ(prob, expected_prob);
}
