#include "PersonTestBase.h"

#include "Core/types.h"
#include "Utils/Constants.h"

class PersonBehaviorTest : public PersonTestBase {};

TEST_F(PersonBehaviorTest, AgeDependentBitingFactorCoversInfantChildAndAdultBands) {
  mock_scheduler_->set_current_time(0);

  person_->set_age(0);
  person_->set_birthday(-10);
  EXPECT_NEAR(person_->get_age_dependent_biting_factor(), 6.5 / 61.5, 1e-12);
  person_->set_birthday(-100);
  EXPECT_NEAR(person_->get_age_dependent_biting_factor(), 8.0 / 61.5, 1e-12);
  person_->set_birthday(-200);
  EXPECT_NEAR(person_->get_age_dependent_biting_factor(), 9.0 / 61.5, 1e-12);
  person_->set_birthday(-400);
  EXPECT_NEAR(person_->get_age_dependent_biting_factor(), 9.5 / 61.5, 1e-12);

  person_->set_age(1);
  EXPECT_NEAR(person_->get_age_dependent_biting_factor(), 11.0 / 61.5, 1e-12);
  person_->set_age(2);
  EXPECT_NEAR(person_->get_age_dependent_biting_factor(), 13.5 / 61.5, 1e-12);
  person_->set_age(3);
  EXPECT_NEAR(person_->get_age_dependent_biting_factor(), 15.5 / 61.5, 1e-12);
  person_->set_age(10);
  EXPECT_NEAR(person_->get_age_dependent_biting_factor(), (17.5 + 6.0 * 2.75) / 61.5, 1e-12);
  person_->set_age(20);
  EXPECT_DOUBLE_EQ(person_->get_age_dependent_biting_factor(), 1.0);
}

TEST_F(PersonBehaviorTest, ProbabilityPresentAtMdaUsesAgeBracket) {
  auto mda = mock_config_->get_strategy_parameters().get_mda();
  mda.set_age_bracket_prob_individual_present_at_mda({5, 15, 100});
  mock_config_->get_strategy_parameters().set_mass_drug_administration(mda);
  person_->set_prob_present_at_mda_by_age({0.2, 0.5, 0.9});

  person_->set_age(4);
  EXPECT_DOUBLE_EQ(person_->prob_present_at_mda(), 0.2);
  person_->set_age(10);
  EXPECT_DOUBLE_EQ(person_->prob_present_at_mda(), 0.5);
  person_->set_age(25);
  EXPECT_DOUBLE_EQ(person_->prob_present_at_mda(), 0.9);
}
