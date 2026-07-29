#include "Events/MoveParasiteToBloodEvent.h"
#include "Parasites/Genotype.h"
#include "PersonParasiteTestBase.h"

TEST_F(PersonParasiteTest, ChangeStateWhenNoParasiteInBloodWhenNoParasiteInLiver) {
  EXPECT_CALL(*mock_population_, notify_change(_, Person::Property::HOST_STATE, _, _)).Times(2);

  person_->set_host_state(Person::HostStates::CLINICAL);
  person_->change_state_when_no_parasite_in_blood();
  EXPECT_EQ(person_->get_host_state(), Person::HostStates::SUSCEPTIBLE);
}

TEST_F(PersonParasiteTest, ChangeStateWhenNoParasiteInBloodWhenParasiteInLiver) {
  EXPECT_CALL(*mock_population_, notify_change(_, Person::Property::HOST_STATE, _, _)).Times(2);

  person_->set_host_state(Person::HostStates::CLINICAL);
  auto genotype_ptr = std::make_unique<Genotype>("aaabbbccc");
  genotype_ptr->set_genotype_id(1);
  auto genotype_raw = genotype_ptr.get();
  Model::get_genotype_db()->add(std::move(genotype_ptr));
  person_->set_liver_parasite_type(genotype_raw);

  person_->change_state_when_no_parasite_in_blood();
  EXPECT_EQ(person_->get_host_state(), Person::HostStates::EXPOSED);
}

TEST_F(PersonParasiteTest, RandomlyChooseParasiteFrom2Infections) {
  EXPECT_CALL(*mock_random_, random_uniform(2)).WillOnce(::testing::Return(1));
  EXPECT_CALL(*mock_population_, notify_change(_, Person::Property::HOST_STATE, _, _)).Times(1);

  auto genotype1_ptr = std::make_unique<Genotype>("aaabbbccc");
  genotype1_ptr->set_genotype_id(1);
  auto genotype2_ptr = std::make_unique<Genotype>("aaabbbccc");
  genotype2_ptr->set_genotype_id(2);
  auto genotype2_raw = genotype2_ptr.get();
  Model::get_genotype_db()->add(std::move(genotype1_ptr));
  Model::get_genotype_db()->add(std::move(genotype2_ptr));

  person_->get_today_infections().push_back(1);
  person_->get_today_infections().push_back(2);
  person_->randomly_choose_parasite();

  EXPECT_EQ(person_->get_host_state(), Person::HostStates::EXPOSED);
  EXPECT_EQ(person_->liver_parasite_type(), genotype2_raw);
  EXPECT_TRUE(person_->get_today_infections().empty());
  ASSERT_EQ(person_->get_events().size(), 1);
  auto event = dynamic_cast<MoveParasiteToBloodEvent*>(person_->get_events().begin()->second.get());
  ASSERT_NE(event, nullptr);
  EXPECT_EQ(event->infection_genotype(), genotype2_raw);
  EXPECT_EQ(event->get_time(), current_time_ + 7);
}

TEST_F(PersonParasiteTest, RandomlyChooseParasiteFrom1Infection) {
  EXPECT_CALL(*mock_random_, random_uniform(1)).Times(0);
  EXPECT_CALL(*mock_population_, notify_change(_, Person::Property::HOST_STATE, _, _)).Times(1);

  auto genotype_ptr = std::make_unique<Genotype>("aaabbbccc");
  genotype_ptr->set_genotype_id(1);
  auto genotype_raw = genotype_ptr.get();
  Model::get_genotype_db()->add(std::move(genotype_ptr));
  person_->get_today_infections().push_back(1);

  person_->randomly_choose_parasite();

  EXPECT_EQ(person_->get_host_state(), Person::HostStates::EXPOSED);
  EXPECT_EQ(person_->liver_parasite_type(), genotype_raw);
  EXPECT_TRUE(person_->get_today_infections().empty());
  ASSERT_EQ(person_->get_events().size(), 1);
  auto event = dynamic_cast<MoveParasiteToBloodEvent*>(person_->get_events().begin()->second.get());
  ASSERT_NE(event, nullptr);
  EXPECT_EQ(event->infection_genotype(), genotype_raw);
  EXPECT_EQ(event->get_time(), current_time_ + 7);
}

TEST_F(PersonParasiteTest, RandomlyChooseParasiteFrom0Infections) {
  EXPECT_CALL(*mock_random_, random_uniform(0)).Times(0);
  EXPECT_CALL(*mock_population_, notify_change(_, Person::Property::HOST_STATE, _, _)).Times(0);

  person_->randomly_choose_parasite();

  EXPECT_EQ(person_->get_host_state(), Person::HostStates::SUSCEPTIBLE);
  EXPECT_TRUE(person_->get_today_infections().empty());
  EXPECT_TRUE(person_->get_events().empty());
}

TEST_F(PersonParasiteTest, RandomlyChooseParasiteFrom4Infections) {
  EXPECT_CALL(*mock_random_, random_uniform(4)).WillOnce(::testing::Return(2));
  EXPECT_CALL(*mock_population_, notify_change(_, Person::Property::HOST_STATE, _, _)).Times(1);

  auto genotype1_ptr = std::make_unique<Genotype>("aaabbbccc");
  genotype1_ptr->set_genotype_id(1);
  auto genotype2_ptr = std::make_unique<Genotype>("aaabbbccc");
  genotype2_ptr->set_genotype_id(2);
  auto genotype3_ptr = std::make_unique<Genotype>("aaabbbccc");
  genotype3_ptr->set_genotype_id(3);
  auto genotype4_ptr = std::make_unique<Genotype>("aaabbbccc");
  genotype4_ptr->set_genotype_id(4);
  auto genotype3_raw = genotype3_ptr.get();
  Model::get_genotype_db()->add(std::move(genotype1_ptr));
  Model::get_genotype_db()->add(std::move(genotype2_ptr));
  Model::get_genotype_db()->add(std::move(genotype3_ptr));
  Model::get_genotype_db()->add(std::move(genotype4_ptr));

  person_->get_today_infections() = {1, 2, 3, 4};
  person_->randomly_choose_parasite();

  EXPECT_EQ(person_->get_host_state(), Person::HostStates::EXPOSED);
  EXPECT_EQ(person_->liver_parasite_type(), genotype3_raw);
  EXPECT_TRUE(person_->get_today_infections().empty());
  ASSERT_EQ(person_->get_events().size(), 1);
  auto event = dynamic_cast<MoveParasiteToBloodEvent*>(person_->get_events().begin()->second.get());
  ASSERT_NE(event, nullptr);
  EXPECT_EQ(event->infection_genotype(), genotype3_raw);
  EXPECT_EQ(event->get_time(), current_time_ + 7);
}
