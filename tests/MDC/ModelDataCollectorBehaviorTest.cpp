#include <gtest/gtest.h>

#include "MDC/ModelDataCollector.h"
#include "Parasites/Genotype.h"
#include "Simulation/Model.h"
#include "Utils/Cli.h"
#include "fixtures/TestFileGenerators.h"

class ModelDataCollectorBehaviorTest : public ::testing::Test {
 protected:
  void SetUp() override {
    test_fixtures::setup_test_environment();
    utils::Cli::MaSimAppInput cli_input;
    cli_input.input_path = "test_input.yml";
    Model::set_cli_input(cli_input);
    ASSERT_TRUE(Model::get_instance()->initialize());
    collector_ = Model::get_mdc();
    collector_->initialize();
    collector_->set_recording(true);
  }

  void TearDown() override {
    Model::get_instance()->release();
    test_fixtures::cleanup_test_files();
  }

  ModelDataCollector* collector_{nullptr};
};

TEST_F(ModelDataCollectorBehaviorTest, RecordsTreatmentFailureDeathInfectionAndRecrudescence) {
  constexpr core::LocationId location = 0;
  collector_->record_1_recrudescence_treatment(location, 85, 1, 0);
  collector_->record_1_non_treated_case(location, 85, 1);
  collector_->record_1_malaria_death(location, 85, true);
  collector_->record_1_malaria_death(location, 25, false);
  collector_->record_1_infection(location);
  collector_->record_1_treatment_failure_by_therapy(location, 1, 0);
  collector_->record_1_treatment_success_by_therapy(location, 1, 0);
  collector_->update_utl_vector();

  EXPECT_EQ(collector_->monthly_number_of_recrudescence_treatment_by_location()[location], 1);
  EXPECT_EQ(collector_->monthly_nontreatment_by_location()[location], 1);
  EXPECT_EQ(collector_->monthly_number_of_new_infections_by_location()[location], 1);
  EXPECT_EQ(collector_->monthly_treatment_failure_by_location()[location], 1);
  EXPECT_EQ(collector_->monthly_treatment_success_by_location()[location], 1);
  EXPECT_NO_THROW(collector_->monthly_update());
}

TEST_F(ModelDataCollectorBehaviorTest, ComputesYearlyAndPartialYearEir) {
  constexpr core::LocationId location = 0;
  auto* scheduler = Model::get_scheduler();
  scheduler->set_current_time(0);
  collector_->yearly_update();

  collector_->collect_number_of_bites(location, 365);
  collector_->update_person_days_by_years(location, 365);
  scheduler->set_current_time(365);
  EXPECT_NO_THROW(collector_->yearly_update());
  EXPECT_NO_THROW(collector_->calculate_eir());
  EXPECT_NO_THROW(collector_->update_after_run());
}

TEST_F(ModelDataCollectorBehaviorTest, HandlesBloodSlideRangesAndMutationByDrug) {
  constexpr core::LocationId location = 0;
  EXPECT_DOUBLE_EQ(collector_->get_blood_slide_prevalence(location, 0, 5), 0.0);
  EXPECT_DOUBLE_EQ(collector_->get_blood_slide_prevalence(location, 0, 20), 0.0);
  EXPECT_DOUBLE_EQ(collector_->get_blood_slide_prevalence(location, 10, 20), 0.0);

  auto from = Model::get_genotype_db()->at(0);
  auto to = Model::get_genotype_db()->at(0);
  EXPECT_NO_THROW(collector_->record_1_mutation_by_drug(location, from, to, 0));
}
