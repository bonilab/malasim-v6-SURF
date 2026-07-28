#include <gtest/gtest.h>

#include "Configuration/Config.h"
#include "Events/Population/AnnualBetaUpdateEvent.hxx"
#include "Events/Population/AnnualCoverageUpdateEvent.hxx"
#include "Events/Population/ChangeCirculationPercentEvent.hxx"
#include "Events/Population/ChangeInterruptedFeedingRateEvent.h"
#include "Events/Population/ChangeTreatmentStrategyEvent.h"
#include "Events/Population/ChangeMutationMaskEvent.h"
#include "Events/Population/ChangeMutationProbabilityPerLocusEvent.h"
#include "Events/Population/ChangeTreatmentCoverageEvent.h"
#include "Events/Population/ChangeWithinHostInducedFreeRecombinationEvent.h"
#include "Events/Population/PopulationEventBuilder.h"
#include "Events/Population/SingleRoundMDAEvent.h"
#include "Events/Population/TurnOffMutationEvent.h"
#include "Events/Population/TurnOnMutationEvent.h"
#include "Events/Population/ChangeCirculationPercentEvent.hxx"
#include "Events/Population/DistrictImportationDailyEvent.h"
#include "Events/Population/ImportationPeriodicallyRandomEvent.h"
#include "Events/Population/ImportationEvent.h"
#include "Events/Population/ImportationPeriodicallyEvent.h"
#include "Events/Population/IntroduceAmodiaquineMutantEvent.h"
#include "Events/Population/IntroduceLumefantrineMutantEvent.h"
#include "Events/Population/Introduce580YMutantEvent.h"
#include "Events/Population/IntroduceTripleMutantToDPMEvent.h"
#include "Events/Population/IntroducePlas2CopyParasiteEvent.h"
#include "Events/Population/IntroduceMutantEvent.hxx"
#include "Events/Population/IntroduceMutantRasterEvent.hxx"
#include "Events/Population/ModifyNestedMFTEvent.h"
#include "Events/Population/UpdateBetaRasterEvent.hxx"
#include "Events/Population/RotateStrategyEvent.h"
#include "Events/Population/IntroduceParasitesPeriodicallyEventV2.h"
#include "Parasites/Genotype.h"
#include "Simulation/Model.h"
#include "Utils/Cli.h"
#include "fixtures/TestFileGenerators.h"

#include <fstream>

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

  EXPECT_THROW(PopulationEventBuilder::build_importation_periodically_random_event(
                   YAML::Load("- date: 2024/01/02\n  genotype_id: 0\n  count: 1\n  log_parasite_density: 3.0"), config()),
               std::invalid_argument);
  EXPECT_THROW(PopulationEventBuilder::build_importation_periodically_random_event(
                   YAML::Load("- date: 2024/01/01\n  genotype_id: -1\n  count: 1\n  log_parasite_density: 3.0"), config()),
               std::invalid_argument);
  EXPECT_THROW(PopulationEventBuilder::build_importation_periodically_random_event(
                   YAML::Load("- date: 2024/01/01\n  genotype_id: 9999\n  count: 1\n  log_parasite_density: 3.0"), config()),
               std::invalid_argument);
  EXPECT_THROW(PopulationEventBuilder::build_importation_periodically_random_event(
                   YAML::Load("- date: 2024/01/01\n  genotype_id: 0\n  count: 0\n  log_parasite_density: 3.0"), config()),
               std::invalid_argument);
  EXPECT_THROW(PopulationEventBuilder::build_importation_periodically_random_event(
                   YAML::Load("- date: 2024/01/01\n  genotype_id: 0\n  count: 1\n  log_parasite_density: 0.0"), config()),
               std::invalid_argument);
}

TEST_F(PopulationEventBuilderTest, BuildsAdditionalPopulationEventTypes) {
  const auto annual_beta = PopulationEventBuilder::build_annual_beta_update_event(
      YAML::Load("- date: 2024/01/01\n  rate: 0.1"), config());
  ASSERT_EQ(annual_beta.size(), 1);
  EXPECT_EQ(annual_beta[0]->name(), AnnualBetaUpdateEvent::EVENT_NAME);

  const auto annual_coverage = PopulationEventBuilder::build_annual_coverage_update_event(
      YAML::Load("- date: 2024/01/01\n  rate: 0.1"), config());
  ASSERT_EQ(annual_coverage.size(), 1);
  EXPECT_EQ(annual_coverage[0]->name(), AnnualCoverageUpdateEvent::EVENT_NAME);

  const auto circulation = PopulationEventBuilder::build_change_circulation_percent_event(
      YAML::Load("- date: 2024/01/01\n  circulation_percent: 0.4"), config());
  ASSERT_EQ(circulation.size(), 1);
  EXPECT_EQ(circulation[0]->name(), ChangeCirculationPercentEvent::EVENT_NAME);

  const std::string raster_name = "builder_beta.asc";
  test_fixtures::create_test_raster_file(raster_name, 1, 1, 0.5);
  const auto raster = PopulationEventBuilder::build_update_beta_raster_event(
      YAML::Load("- date: 2024/01/01\n  beta_raster: builder_beta.asc"), config());
  ASSERT_EQ(raster.size(), 1);
  EXPECT_EQ(raster[0]->name(), UpdateBetaRasterEvent::EVENT_NAME);
  std::remove(raster_name.c_str());

  const auto district = PopulationEventBuilder::build_import_district_mutant_daily_events(
      YAML::Load("- district: 1\n  daily_rate: 0.2\n  start_date: 2024/01/01\n  alleles:\n    - chromosome: 5\n      locus: 86\n      allele: Y"),
      config());
  ASSERT_EQ(district.size(), 1);
  EXPECT_EQ(district[0]->name(), DistrictImportationDailyEvent::EVENT_NAME);

  const auto random_import = PopulationEventBuilder::build_importation_periodically_random_event(
      YAML::Load("- date: 2024/01/01\n  genotype_id: 0\n  count: 2\n  log_parasite_density: 3.0"),
      config());
  ASSERT_EQ(random_import.size(), 1);
  EXPECT_EQ(random_import[0]->name(), ImportationPeriodicallyRandomEvent::EVENT_NAME);
}

TEST_F(PopulationEventBuilderTest, BuildsMutantEventDefinitionsAndSkipsInvalidLocations) {
  const auto alleles = "    - chromosome: 5\n      locus: 86\n      allele: Y\n";
  const auto amodiaquine = PopulationEventBuilder::build_introduce_amodiaquine_mutant_parasite_events(
      YAML::Load(std::string("- location: 0\n  date: 2024/01/01\n  fraction: 0.2\n  alleles:\n") + alleles), config());
  ASSERT_EQ(amodiaquine.size(), 1);
  EXPECT_EQ(amodiaquine[0]->name(), IntroduceAmodiaquineMutantEvent::EVENT_NAME);

  const auto lumefantrine = PopulationEventBuilder::build_introduce_lumefantrine_mutant_parasite_events(
      YAML::Load(std::string("- location: 0\n  date: 2024/01/01\n  fraction: 0.2\n  alleles:\n") + alleles), config());
  ASSERT_EQ(lumefantrine.size(), 1);
  EXPECT_EQ(lumefantrine[0]->name(), IntroduceLumefantrineMutantEvent::EVENT_NAME);

  const auto mutant_580y = PopulationEventBuilder::build_introduce_580Y_mutant_events(
      YAML::Load(std::string("- location: 0\n  date: 2024/01/01\n  fraction: 0.2\n  alleles:\n") + alleles), config());
  ASSERT_EQ(mutant_580y.size(), 1);
  EXPECT_EQ(mutant_580y[0]->name(), Introduce580YMutantEvent::EVENT_NAME);

  const auto triple = PopulationEventBuilder::build_introduce_triple_mutant_to_dpm_parasite_events(
      YAML::Load(std::string("- location: 0\n  date: 2024/01/01\n  fraction: 0.2\n  alleles:\n") + alleles), config());
  ASSERT_EQ(triple.size(), 1);
  EXPECT_EQ(triple[0]->name(), IntroduceTripleMutantToDPMEvent::EVENT_NAME);

  const auto skipped = PopulationEventBuilder::build_introduce_580Y_mutant_events(
      YAML::Load("- location: 999\n  date: 2024/01/01\n  fraction: 0.2\n  alleles: []"), config());
  EXPECT_TRUE(skipped.empty());
}

TEST_F(PopulationEventBuilderTest, BuildsImportationAndPeriodicDefinitions) {
  const std::string genotype = "||||YF1||TTHFIMG,x||||||FNCMYRIPRPCRA|1";
  const auto importation = PopulationEventBuilder::build_introduce_parasite_events(
      YAML::Load("- location: 0\n  parasite_info:\n    - genotype_aa_sequence: \"" + genotype + "\""
                   + "\n      number_of_cases: 3\n      date: 2024/01/01"),
      config());
  ASSERT_EQ(importation.size(), 1);
  auto* importation_event = dynamic_cast<ImportationEvent*>(importation[0].get());
  ASSERT_NE(importation_event, nullptr);
  EXPECT_EQ(importation_event->get_location(), 0);
  EXPECT_EQ(importation_event->get_number_of_cases(), 3);

  const auto periodic = PopulationEventBuilder::build_introduce_parasites_periodically_events(
      YAML::Load("- location: 0\n  parasite_info:\n    - genotype_aa_sequence: \"" + genotype + "\""
                   + "\n      number_of_cases: 2\n      duration: 30\n      start_date: 2024/01/01"),
      config());
  ASSERT_EQ(periodic.size(), 1);
  auto* periodic_event = dynamic_cast<ImportationPeriodicallyEvent*>(periodic[0].get());
  ASSERT_NE(periodic_event, nullptr);
  EXPECT_EQ(periodic_event->get_location(), 0);
  EXPECT_EQ(periodic_event->get_duration(), 30);
  EXPECT_EQ(periodic_event->get_number_of_cases(), 2);

  const auto plas2 = PopulationEventBuilder::build_introduce_plas2_parasite_events(
      YAML::Load("- location: 0\n  date: 2024/01/01\n  fraction: 0.4"), config());
  ASSERT_EQ(plas2.size(), 1);
  EXPECT_EQ(plas2[0]->name(), IntroducePlas2CopyParasiteEvent::EVENT_NAME);

  const auto nested = PopulationEventBuilder::build_modify_nested_mft_strategy_event(
      YAML::Load("- date: 2024/01/01\n  strategy_id: 0"), config());
  ASSERT_EQ(nested.size(), 1);
  EXPECT_EQ(nested[0]->name(), ModifyNestedMFTEvent::EVENT_NAME);
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

TEST_F(PopulationEventBuilderTest, ExecutesSimpleMutationAndMosquitoConfigurationEvents) {
  auto turn_on = std::make_unique<TurnOnMutationEvent>(0, 0.25);
  turn_on->set_executable(true);
  turn_on->execute();
  EXPECT_DOUBLE_EQ(config()->get_genotype_parameters().get_mutation_probability_per_locus(), 0.25);

  auto turn_off = std::make_unique<TurnOffMutationEvent>(0);
  turn_off->set_executable(true);
  turn_off->execute();
  EXPECT_DOUBLE_EQ(config()->get_genotype_parameters().get_mutation_probability_per_locus(), 0.0);

  auto recombination = std::make_unique<ChangeWithinHostInducedFreeRecombinationEvent>(false, 0);
  recombination->set_executable(true);
  recombination->execute();
  EXPECT_FALSE(config()->get_mosquito_parameters().get_within_host_induced_free_recombination());

  auto mask = std::make_unique<ChangeMutationMaskEvent>(std::vector<bool>{true, false}, 0);
  mask->set_executable(true);
  mask->execute();
  EXPECT_EQ(config()->get_genotype_parameters().get_mutation_mask(), std::vector<bool>({true, false}));

  auto ifr = std::make_unique<ChangeInterruptedFeedingRateEvent>(0, 0.42, 0);
  ifr->set_executable(true);
  ifr->execute();
  EXPECT_DOUBLE_EQ(config()->location_db()[0].mosquito_ifr, 0.42);

  auto probability = std::make_unique<ChangeMutationProbabilityPerLocusEvent>(0.33, 0);
  probability->set_executable(true);
  probability->execute();
  EXPECT_DOUBLE_EQ(config()->get_genotype_parameters().get_mutation_probability_per_locus(), 0.33);
}

TEST_F(PopulationEventBuilderTest, ExecutesAnnualAndCirculationUpdates) {
  config()->location_db()[0].beta = 1.0;
  auto annual_beta = std::make_unique<AnnualBetaUpdateEvent>(0.5, 0);
  annual_beta->set_executable(true);
  annual_beta->execute();
  EXPECT_DOUBLE_EQ(config()->location_db()[0].beta, 1.5);

  auto* treatment_coverage = Model::get_treatment_coverage();
  ASSERT_NE(treatment_coverage, nullptr);
  treatment_coverage->p_treatment_under_5[0] = 0.2;
  treatment_coverage->p_treatment_over_5[0] = 0.4;
  auto annual_coverage = std::make_unique<AnnualCoverageUpdateEvent>(0.5, 0);
  annual_coverage->set_executable(true);
  annual_coverage->execute();
  EXPECT_DOUBLE_EQ(treatment_coverage->p_treatment_under_5[0], 0.6);
  EXPECT_DOUBLE_EQ(treatment_coverage->p_treatment_over_5[0], 0.7);

  auto circulation_info = config()->get_movement_settings().get_circulation_info();
  circulation_info.set_circulation_percent(0.1);
  config()->get_movement_settings().set_circulation_info(circulation_info);
  auto circulation = std::make_unique<ChangeCirculationPercentEvent>(0.8, 0);
  circulation->set_executable(true);
  circulation->execute();
  EXPECT_NEAR(config()->get_movement_settings().get_circulation_info().get_circulation_percent(),
              0.8, 1e-6);
}

TEST_F(PopulationEventBuilderTest, BuildsMutantEventsFromUnitsAndRasters) {
  const auto unit_events = PopulationEventBuilder::build_introduce_mutant_event(
      YAML::Load(R"(
        - day: 2024/01/10
          unit_id: 0
          fraction: 0.25
          alleles:
            - chromosome: 1
              locus: 2
              allele: A
      )"), config(), "district");
  ASSERT_EQ(unit_events.size(), 1);
  EXPECT_EQ(unit_events[0]->name(), IntroduceMutantEvent::EVENT_NAME);
  unit_events[0]->set_executable(true);
  EXPECT_NO_THROW(unit_events[0]->execute());

  {
    std::ofstream raster("mutant_event.asc");
    raster << "ncols 8\n"
           << "nrows 1\n"
           << "xllcorner 0\n"
           << "yllcorner 0\n"
           << "cellsize 1\n"
           << "NODATA_value -9999\n"
           << "1 0 1 0 1 0 1 0\n";
  }
  const auto raster_events = PopulationEventBuilder::build_introduce_mutant_raster_event(
      YAML::Load(R"(
        - date: 2024/01/11
          raster: mutant_event.asc
          fraction: 0.5
          alleles:
            - chromosome: 1
              locus: 3
              allele: T
      )"), config());
  ASSERT_EQ(raster_events.size(), 1);
  EXPECT_EQ(raster_events[0]->name(), IntroduceMutantRasterEvent::EVENT_NAME);
  raster_events[0]->set_executable(true);
  EXPECT_NO_THROW(raster_events[0]->execute());
}

TEST_F(PopulationEventBuilderTest, ExecutesZeroCaseImportationEventsAndReschedules) {
  auto importation = std::make_unique<ImportationEvent>(0, 0, 0, 0);
  importation->set_executable(true);
  EXPECT_NO_THROW(importation->execute());

  auto periodic = std::make_unique<ImportationPeriodicallyEvent>(0, 1, 0, 0, 0);
  periodic->set_executable(true);
  EXPECT_NO_THROW(periodic->execute());

  auto random_periodic = std::make_unique<ImportationPeriodicallyRandomEvent>(0, 0, 0, 1.0);
  random_periodic->set_executable(true);
  EXPECT_NO_THROW(random_periodic->execute());
}

TEST_F(PopulationEventBuilderTest, ExecutesStrategyAndDistrictEventsWithNoCases) {
  auto strategy_change = std::make_unique<ChangeTreatmentStrategyEvent>(0, 0);
  strategy_change->set_executable(true);
  EXPECT_NO_THROW(strategy_change->execute());
  ASSERT_NE(Model::get_treatment_strategy(), nullptr);
  EXPECT_EQ(Model::get_treatment_strategy()->id, 0);

  auto district = std::make_unique<DistrictImportationDailyEvent>(0, 0.0, 0,
                                                                   std::vector<std::tuple<int, int, char>>{});
  district->set_executable(true);
  EXPECT_NO_THROW(district->execute());
}

TEST_F(PopulationEventBuilderTest, DispatchesRemainingSupportedEventDefinitions) {
  const auto genotype_sequence = Model::get_genotype_db()->at(0)->get_aa_sequence();

  YAML::Node importation;
  importation["name"] = ImportationEvent::EVENT_NAME;
  importation["info"][0]["location"] = 0;
  importation["info"][0]["parasite_info"][0]["genotype_aa_sequence"] = genotype_sequence;
  importation["info"][0]["parasite_info"][0]["number_of_cases"] = 0;
  importation["info"][0]["parasite_info"][0]["date"] = "2024/01/01";
  EXPECT_EQ(PopulationEventBuilder::build(importation).size(), 1);

  YAML::Node periodic;
  periodic["name"] = ImportationPeriodicallyEvent::EVENT_NAME;
  periodic["info"][0]["location"] = 0;
  periodic["info"][0]["parasite_info"][0]["genotype_aa_sequence"] = genotype_sequence;
  periodic["info"][0]["parasite_info"][0]["duration"] = 1;
  periodic["info"][0]["parasite_info"][0]["number_of_cases"] = 0;
  periodic["info"][0]["parasite_info"][0]["start_date"] = "2024/01/01";
  EXPECT_EQ(PopulationEventBuilder::build(periodic).size(), 1);

  YAML::Node random_periodic;
  random_periodic["name"] = ImportationPeriodicallyRandomEvent::EVENT_NAME;
  random_periodic["info"][0]["date"] = "2024/02/01";
  random_periodic["info"][0]["genotype_id"] = 0;
  random_periodic["info"][0]["count"] = 1;
  random_periodic["info"][0]["log_parasite_density"] = 1.0;
  EXPECT_EQ(PopulationEventBuilder::build(random_periodic).size(), 1);

  YAML::Node annual_beta;
  annual_beta["name"] = AnnualBetaUpdateEvent::EVENT_NAME;
  annual_beta["info"][0]["date"] = "2024/01/01";
  annual_beta["info"][0]["rate"] = 0.1;
  EXPECT_EQ(PopulationEventBuilder::build(annual_beta).size(), 1);

  YAML::Node mask;
  mask["name"] = ChangeMutationMaskEvent::EVENT_NAME;
  mask["info"][0]["date"] = "2024/01/01";
  mask["info"][0]["mutation_mask"] = YAML::Load("[true, false]");
  EXPECT_EQ(PopulationEventBuilder::build(mask).size(), 1);

  YAML::Node unknown;
  unknown["name"] = "unsupported_event";
  EXPECT_TRUE(PopulationEventBuilder::build(unknown).empty());
}

TEST_F(PopulationEventBuilderTest, DispatchesConfigurationAndRasterEventDefinitions) {
  YAML::Node coverage;
  coverage["name"] = ChangeTreatmentCoverageEvent::EVENT_NAME;
  coverage["info"][0] = YAML::Load(R"(
type: SteadyTCM
date: 2000/1/1
p_treatment_under_5_by_location: [0.1]
p_treatment_over_5_by_location: [0.2]
)");
  EXPECT_EQ(PopulationEventBuilder::build(coverage).size(), 1);

  YAML::Node strategy;
  strategy["name"] = ChangeTreatmentStrategyEvent::EVENT_NAME;
  strategy["info"][0]["date"] = "2000/1/1";
  strategy["info"][0]["strategy_id"] = 0;
  EXPECT_EQ(PopulationEventBuilder::build(strategy).size(), 1);

  YAML::Node single_mda;
  single_mda["name"] = SingleRoundMDAEvent::EVENT_NAME;
  single_mda["info"][0]["date"] = "2000/1/1";
  single_mda["info"][0]["fraction_population_targeted"] = YAML::Load("[0.0]");
  single_mda["info"][0]["days_to_complete_all_treatments"] = 1;
  EXPECT_EQ(PopulationEventBuilder::build(single_mda).size(), 1);

  YAML::Node nested;
  nested["name"] = ModifyNestedMFTEvent::EVENT_NAME;
  nested["info"][0]["date"] = "2000/1/1";
  nested["info"][0]["strategy_id"] = 0;
  EXPECT_EQ(PopulationEventBuilder::build(nested).size(), 1);

  YAML::Node raster;
  raster["name"] = UpdateBetaRasterEvent::EVENT_NAME;
  raster["info"][0]["date"] = "2000/1/1";
  raster["info"][0]["beta_raster"] = "builder_beta_dispatch.asc";
  test_fixtures::create_test_raster_file("builder_beta_dispatch.asc", 1, 1, 0.5);
  EXPECT_EQ(PopulationEventBuilder::build(raster).size(), 1);
  std::remove("builder_beta_dispatch.asc");

  YAML::Node rotate;
  rotate["name"] = RotateStrategyEvent::EVENT_NAME;
  rotate["info"][0]["date"] = "2000/1/1";
  rotate["info"][0]["years"] = 1;
  rotate["info"][0]["first_strategy_id"] = 0;
  rotate["info"][0]["second_strategy_id"] = 0;
  EXPECT_EQ(PopulationEventBuilder::build(rotate).size(), 1);

  YAML::Node district;
  district["name"] = DistrictImportationDailyEvent::EVENT_NAME;
  district["info"][0]["district"] = 1;
  district["info"][0]["daily_rate"] = 0.1;
  district["info"][0]["start_date"] = "2000/1/1";
  district["info"][0]["alleles"] = YAML::Load("[]");
  EXPECT_EQ(PopulationEventBuilder::build(district).size(), 1);

  YAML::Node mutant;
  mutant["name"] = IntroduceMutantEvent::EVENT_NAME;
  mutant["admin_level"] = "district";
  mutant["info"] = YAML::Load(R"(
- day: 2000/1/1
  unit_id: 0
  fraction: 0.1
  alleles:
    - chromosome: 1
      locus: 2
      allele: A
)");
  EXPECT_EQ(PopulationEventBuilder::build(mutant).size(), 1);
}
