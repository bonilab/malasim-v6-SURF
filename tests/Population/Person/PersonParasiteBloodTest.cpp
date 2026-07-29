#include "Events/MoveParasiteToBloodEvent.h"
#include "Parasites/Genotype.h"
#include "PersonParasiteTestBase.h"

TEST_F(PersonParasiteTest, AddNewParasiteToBlood) {
  auto genotype_ptr = std::make_unique<Genotype>("aaabbbccc");
  genotype_ptr->set_genotype_id(1);

  auto parasite = person_->add_new_parasite_to_blood(genotype_ptr.get());
  EXPECT_NE(parasite, nullptr);
}

TEST_F(PersonParasiteTest, InfectedByWhenLiverParasiteTypeIsNotSet) {
  EXPECT_CALL(*mock_population_, notify_change(_, Person::Property::HOST_STATE, _, _)).Times(1);

  auto genotype_ptr = std::make_unique<Genotype>("aaabbbccc");
  genotype_ptr->set_genotype_id(1);
  auto genotype_raw = genotype_ptr.get();
  Model::get_genotype_db()->add(std::move(genotype_ptr));

  person_->infected_by(1);

  EXPECT_EQ(person_->get_host_state(), Person::HostStates::EXPOSED);
  EXPECT_EQ(person_->liver_parasite_type(), genotype_raw);
  EXPECT_EQ(person_->get_events().size(), 1);

  auto event = dynamic_cast<MoveParasiteToBloodEvent*>(person_->get_events().begin()->second.get());
  ASSERT_NE(event, nullptr);
  EXPECT_EQ(event->infection_genotype(), genotype_raw);
  EXPECT_EQ(event->get_time(), current_time_ + 7);
}

TEST_F(PersonParasiteTest, InfectedByWhenLiverParasiteTypeIsSet) {
  EXPECT_CALL(*mock_population_, notify_change(_, Person::Property::HOST_STATE, _, _)).Times(1);

  auto genotype_ptr = std::make_unique<Genotype>("aaabbbccc");
  genotype_ptr->set_genotype_id(1);
  auto genotype_raw = genotype_ptr.get();
  Model::get_genotype_db()->add(std::move(genotype_ptr));

  person_->infected_by(1);
  EXPECT_EQ(person_->get_host_state(), Person::HostStates::EXPOSED);
  EXPECT_EQ(person_->liver_parasite_type(), genotype_raw);

  auto new_genotype_ptr = std::make_unique<Genotype>("abcabcabc");
  new_genotype_ptr->set_genotype_id(2);
  Model::get_genotype_db()->add(std::move(new_genotype_ptr));

  person_->infected_by(2);

  EXPECT_EQ(person_->get_host_state(), Person::HostStates::EXPOSED);
  EXPECT_EQ(person_->liver_parasite_type(), genotype_raw);
  EXPECT_EQ(person_->get_events().size(), 1);

  auto event = dynamic_cast<MoveParasiteToBloodEvent*>(person_->get_events().begin()->second.get());
  ASSERT_NE(event, nullptr);
  EXPECT_EQ(event->infection_genotype(), genotype_raw);
  EXPECT_EQ(event->get_time(), current_time_ + 7);
}

TEST_F(PersonParasiteTest, HasDetectableParasiteWhenNoParasiteInBlood) {
  EXPECT_FALSE(person_->has_detectable_parasite());
}

TEST_F(PersonParasiteTest, HasDetectableParasiteWhenParasiteInBlood) {
  auto genotype_ptr = std::make_unique<Genotype>("aaabbbccc");
  genotype_ptr->set_genotype_id(1);
  auto genotype_raw = genotype_ptr.get();
  Model::get_genotype_db()->add(std::move(genotype_ptr));

  auto parasite = person_->add_new_parasite_to_blood(genotype_raw);
  parasite->set_last_update_log10_parasite_density(2.0);
  EXPECT_TRUE(person_->has_detectable_parasite());
}

TEST_F(PersonParasiteTest, HasParasiteInBloodsButUnderDetectionThreshold) {
  auto genotype_ptr = std::make_unique<Genotype>("aaabbbccc");
  genotype_ptr->set_genotype_id(1);
  auto genotype_raw = genotype_ptr.get();
  Model::get_genotype_db()->add(std::move(genotype_ptr));

  auto parasite = person_->add_new_parasite_to_blood(genotype_raw);
  parasite->set_last_update_log10_parasite_density(-2.0);
  EXPECT_FALSE(person_->has_detectable_parasite());
}

TEST_F(PersonParasiteTest, IsGametocytaemicWhenParasiteInBloodAndCanProduceGametocytes) {
  EXPECT_FALSE(person_->is_gametocytaemic());

  auto genotype_ptr = std::make_unique<Genotype>("aaabbbccc");
  genotype_ptr->set_genotype_id(1);
  auto genotype_raw = genotype_ptr.get();
  Model::get_genotype_db()->add(std::move(genotype_ptr));

  auto parasite = person_->add_new_parasite_to_blood(genotype_raw);
  parasite->set_gametocyte_level(1.0);
  EXPECT_TRUE(person_->is_gametocytaemic());
}

TEST_F(PersonParasiteTest, NotGametocytaemicWhenNoParasiteInBlood) {
  EXPECT_FALSE(person_->is_gametocytaemic());
}

TEST_F(PersonParasiteTest, NotGametocytaemicWhenParasiteInBloodButNoGametocytes) {
  auto genotype_ptr = std::make_unique<Genotype>("aaabbbccc");
  genotype_ptr->set_genotype_id(1);
  auto genotype_raw = genotype_ptr.get();
  Model::get_genotype_db()->add(std::move(genotype_ptr));

  auto parasite = person_->add_new_parasite_to_blood(genotype_raw);
  parasite->set_gametocyte_level(0.0);
  EXPECT_FALSE(person_->is_gametocytaemic());
}
