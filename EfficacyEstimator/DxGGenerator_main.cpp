/*
 * DxGGenerator_main.cpp
 *
 * Alternative executable built against the malaria simulation that can be used
 * to calculate the drug efficacies for the various genotypes.
 *
 * Includes recurrence_test mode ported from TMS.
 */

#include <math.h>
#include <cstdlib>
#include <iomanip>
#include <iostream>

#include "Configuration/Config.h"
#include "Utils/Random.h"
#include "Events/ProgressToClinicalEvent.h"
#include "IndividualsFileReporter.h"
#include "MDC/ModelDataCollector.h"
#include "Simulation/Model.h"
#include "PkPdReporter.h"
#include "Population/ImmuneSystem/ImmuneSystem.h"
#include "Population/Person/Person.h"
#include "Population/Population.h"
#include "Utils/Index/PersonIndexAll.h"
#include "Treatment/Strategies/IStrategy.h"
#include "Treatment/Strategies/SFTStrategy.h"
#include "Treatment/Therapies/SCTherapy.h"
#include "Treatment/Therapies/Therapy.h"
#include "Parasites/Genotype.h"

using namespace std;

bool validate_config_for_ee(utils::Cli::DxGAppInput& input);
double getEfficacyForTherapy(std::string g_str, Model* p_model,
                             utils::Cli::DxGAppInput& input, int therapy_id);
double getEfficacyForTherapyCRT(Model* p_model,
                                utils::Cli::DxGAppInput& input, int therapy_id);
double getEfficacyForTherapyRecurrenceTest(std::string g_str, Model* p_model,
                                           utils::Cli::DxGAppInput& input,
                                           int therapy_id);

inline double round(double val) {
  if (val < 0) return ceil(val - 0.5);
  return floor(val + 0.5);
}

int main(int argc, char** argv) {
  utils::Cli::get_instance().parse(argc, argv);

  // Sync parsed input path into Model's cli_input_ before initialization.
  utils::Cli::MaSimAppInput cli_input{};
  cli_input.input_path = utils::Cli::get_instance().get_input_path();
  Model::set_cli_input(std::move(cli_input));

  Model::get_instance()->initialize();

  auto p_model = Model::get_instance();
  auto input = utils::Cli::get_instance().get_dxg_app_input();

  if (input.as_iiv != -1) {
    for (auto& sd : p_model->get_drug_db()
                        ->at(0)
                        ->age_group_specific_drug_concentration_sd()) {
      sd = utils::Cli::get_instance().get_dxg_app_input().as_iiv;
    }
  }

  if (input.as_ec50 != -1) {
    //TODO: Check this
    p_model->get_drug_db()->at(0)->base_EC50 = input.as_ec50;
  }

  std::cout << std::setprecision(5);
  int max_therapy_id{0}, min_therapy_id{0};

  if (input.therapy_list.empty()) {
    if (input.therapies.empty()) {
      min_therapy_id = 0;
      max_therapy_id = 0;
    } else if (input.therapies.size() == 1) {
      min_therapy_id = input.therapies[0];
      max_therapy_id = input.therapies[0];
    } else if (input.therapies.size() == 2) {
      min_therapy_id = input.therapies[0];
      max_therapy_id = input.therapies[1];
    }
  }

  // ==================== RECURRENCE TEST MODE (from TMS) ====================
  if (input.is_recurrence_test) {
    if (input.genotypes.empty()) {
      std::cout << "List of genotypes is empty for recurrence test"
                << std::endl;
      exit(0);
    }

    // Print header
    std::cout << "ID,Genotype";
    if (input.therapy_list.empty()) {
      for (auto therapy_id = min_therapy_id; therapy_id <= max_therapy_id;
           therapy_id++) {
        std::cout << "," << *p_model->get_therapy_db()[therapy_id];
      }
    } else {
      for (auto therapy_id : input.therapy_list) {
        std::cout << "," << *p_model->get_therapy_db()[therapy_id];
      }
    }
    std::cout << std::endl;

    // Run each genotype x therapy combination
    for (int g_index = 0; g_index < input.genotypes.size(); g_index++) {
      std::stringstream ss;
      ss << g_index << "," << input.genotypes[g_index];

      if (input.therapy_list.empty()) {
        for (auto therapy_id = min_therapy_id; therapy_id <= max_therapy_id;
             therapy_id++) {
          double efficacy = getEfficacyForTherapyRecurrenceTest(
              input.genotypes[g_index], p_model, input, therapy_id);
          ss << "," << efficacy;
        }
      } else {
        for (auto therapy_id : input.therapy_list) {
          double efficacy = getEfficacyForTherapyRecurrenceTest(
              input.genotypes[g_index], p_model, input, therapy_id);
          ss << "," << efficacy;
        }
      }
      std::cout << ss.str() << std::endl;
    }
  }
  // ==================== CRT CALIBRATION MODE ====================
  else if (input.is_crt_calibration) {
    if (input.genotypes.empty()) {
      std::cout << "List of population genotypes is empty" << std::endl;
      exit(0);
    }
    std::stringstream ss;
    if (input.therapy_list.empty()) {
      for (auto therapy_id = min_therapy_id; therapy_id <= max_therapy_id;
           therapy_id++) {
        std::cout << *p_model->get_therapy_db()[therapy_id] << "\t";
      }
    } else {
      for (auto therapy_id : input.therapy_list) {
        std::cout << *p_model->get_therapy_db()[therapy_id] << "\t";
      }
    }
    std::cout << std::endl;
    if (input.therapy_list.empty()) {
      for (auto therapy_id = min_therapy_id; therapy_id <= max_therapy_id;
           therapy_id++) {
        double efficacy = getEfficacyForTherapyCRT(p_model, input, therapy_id);
        ss << efficacy
           << (therapy_id == max_therapy_id ? "" : "\t");
      }
    } else {
      for (int t_index = 0; t_index < input.therapy_list.size(); t_index++) {
        double efficacy = getEfficacyForTherapyCRT(
            p_model, input, input.therapy_list[t_index]);
        ss << efficacy
           << (input.therapy_list[t_index] == input.therapy_list.size() - 1
                   ? ""
                   : "\t");
      }
    }
    std::cout << ss.str() << std::endl;
  }
  // ==================== EE CALIBRATION MODE ====================
  else if (input.is_ee_calibration) {
    if (!validate_config_for_ee(input)) {
      std::cout << "Parameters for Efficacy Estimator are not correct"
                << std::endl;
      exit(0);
    } else if (input.genotypes.size() > 1) {
      std::cout << "Only 1 genotype is accepted using Efficacy Estimator"
                << std::endl;
      exit(0);
    } else {
      if (input.population_size != p_model->get_population()->size()) {
        p_model->get_config()->location_db()[0].population_size =
            input.population_size;
        p_model->set_population(std::make_unique<Population>());
        p_model->get_population()->initialize();
      }

      auto start_drug_id = input.is_art ? 0 : 1;
      for (int i = 0; i < input.number_of_drugs_in_combination; i++) {
        auto* dt = p_model->get_drug_db()->at(i + start_drug_id).get();
        dt->set_name(fmt::format("D{}", i));
        dt->set_drug_half_life(input.half_life[i]);
        dt->set_maximum_parasite_killing_rate(input.k_max[i]);
        dt->set_n(input.slope[i]);
        dt->set_k(4);
        for (double& mda : dt->age_specific_drug_absorption()) {
          mda = input.mean_drug_absorption[i];
        }
      }

      auto* scTherapy =
          static_cast<SCTherapy*>(p_model->get_therapy_db()[0].get());
      scTherapy->drug_ids.clear();
      scTherapy->dosing_day.clear();
      for (int i = 0; i < input.number_of_drugs_in_combination; i++) {
        scTherapy->drug_ids.push_back(i + start_drug_id);
        scTherapy->dosing_day.push_back(input.dosing_days[i]);
      }

      p_model->get_reporters().clear();
      p_model->add_reporter(std::make_unique<PkPdReporter>(&input));

      p_model->get_genotype_db()->clear();
      std::vector<Genotype*> genotype_inputs;
      for (auto genotype_str : input.genotypes) {
        genotype_inputs.push_back(
            p_model->get_genotype_db()->get_genotype(genotype_str));
      }

      auto* genotype = p_model->get_genotype_db()->get_genotype(
          genotype_inputs.front()->get_aa_sequence());

      for (auto& person :
           p_model->get_population()->all_persons()->v_person()) {
        auto density = p_model->get_config()
                           ->get_parasite_parameters()
                           .get_parasite_density_levels()
                           .get_log_parasite_density_from_liver();
        auto* blood_parasite = person->add_new_parasite_to_blood(genotype);

        person->get_immune_system()->set_increase(true);
        person->set_host_state(Person::EXPOSED);

        blood_parasite->set_gametocyte_level(
            p_model->get_config()
                ->get_epidemiological_parameters()
                .get_gametocyte_level_full());
        blood_parasite->set_last_update_log10_parasite_density(density);

        const int days_to_clinical =
            (person->get_age() <= 5)
                ? Model::get_config()
                      ->get_epidemiological_parameters()
                      .get_days_to_clinical_under_five()
                : Model::get_config()
                      ->get_epidemiological_parameters()
                      .get_days_to_clinical_over_five();
        auto event =
            std::make_unique<ProgressToClinicalEvent>(person.get());
        event->set_time(person->calculate_future_time(days_to_clinical));
        event->set_clinical_caused_parasite(blood_parasite);
        person->schedule_basic_event(std::move(event));
      }

      p_model->run();
      const auto result =
          1 - p_model->get_mdc()->blood_slide_prevalence_by_location()[0];
    }
  }
  // ==================== DEFAULT MODE ====================
  else {
    std::cout << "ID\tGenotype\t";
    if (input.therapy_list.empty()) {
      for (auto therapy_id = min_therapy_id; therapy_id <= max_therapy_id;
           therapy_id++) {
        std::cout << *p_model->get_therapy_db()[therapy_id] << "\t";
      }
    } else {
      for (auto therapy_id : input.therapy_list) {
        std::cout << *p_model->get_therapy_db()[therapy_id] << "\t";
      }
    }
    std::cout << std::endl;
    for (int g_index = 0; g_index < input.genotypes.size(); g_index++) {
      std::stringstream ss;
      if (input.is_old_format) {
        ss << g_index << "\t"
           << p_model->get_mosquito()->get_old_genotype_string2(
                  input.genotypes[g_index])
           << "\t";
      } else {
        ss << g_index << "\t"
           << p_model->get_mosquito()->get_old_genotype_string(
                  input.genotypes[g_index])
           << "\t";
      }
      if (input.therapy_list.empty()) {
        for (auto therapy_id = min_therapy_id; therapy_id <= max_therapy_id;
             therapy_id++) {
          double efficacy = getEfficacyForTherapy(
              input.genotypes[g_index], p_model, input, therapy_id);
          ss << efficacy << "\t";
        }
      } else {
        for (int t_index = 0; t_index < input.therapy_list.size();
             t_index++) {
          double efficacy = getEfficacyForTherapy(
              input.genotypes[g_index], p_model, input,
              input.therapy_list[t_index]);
          ss << efficacy << "\t";
        }
      }
      std::cout << ss.str() << std::endl;
    }
  }

  return 0;
}

// ==================== RECURRENCE TEST (TMS logic in malasim API) ====================
double getEfficacyForTherapyRecurrenceTest(std::string g_str, Model* p_model,
                                           utils::Cli::DxGAppInput& input,
                                           int therapy_id) {
  Therapy* mainTherapy = p_model->get_therapy_db()[therapy_id].get();

  // Verify and set the strategy
  auto* strategy =
      dynamic_cast<SFTStrategy*>(p_model->get_treatment_strategy());
  if (strategy == nullptr) {
    std::cerr
        << "The recurrence_test mode can only be used with a SFTStrategy!"
        << std::endl;
    exit(EXIT_FAILURE);
  }
  strategy->get_therapy_list().clear();
  strategy->add_therapy(mainTherapy);

  // Reset reporters and add PkPdReporter with CSV per-person output
  p_model->get_reporters().clear();

  auto pkpd_reporter = std::make_unique<PkPdReporter>(&input);
  auto* pkpd_ptr = pkpd_reporter.get();
  p_model->add_reporter(std::move(pkpd_reporter));

  // Initialize with prefix for per-genotype/therapy output file
  std::string path_prefix = fmt::format("{}_{}", therapy_id, g_str);
  pkpd_ptr->initialize(0, path_prefix);

  // Infect all persons with the specified genotype
  auto* genotype = p_model->get_genotype_db()->get_genotype(g_str);

  for (auto& person :
       p_model->get_population()->all_persons()->v_person()) {
    person->get_all_clonal_parasite_populations()->clear();
    auto density = Model::get_config()
                       ->get_parasite_parameters()
                       .get_parasite_density_levels()
                       .get_log_parasite_density_from_liver();
    auto* blood_parasite = person->add_new_parasite_to_blood(genotype);

    person->get_immune_system()->set_increase(true);
    person->set_host_state(Person::EXPOSED);

    blood_parasite->set_gametocyte_level(
        Model::get_config()
            ->get_epidemiological_parameters()
            .get_gametocyte_level_full());
    blood_parasite->set_last_update_log10_parasite_density(density);

    auto event =
        std::make_unique<ProgressToClinicalEvent>(person.get());
    event->set_time(0);
    event->set_clinical_caused_parasite(blood_parasite);
    person->schedule_basic_event(std::move(event));
  }

  p_model->run();

  // TMS uses treatment failure rate instead of blood slide prevalence
  const auto result =
      1 - (1.0 *
           p_model->get_mdc()->monthly_treatment_failure_by_location()[0] /
           p_model->get_population()->size());

  // Reset simulation state for next run
  p_model->set_population(std::make_unique<Population>());
  p_model->set_scheduler(std::make_unique<Scheduler>());

  p_model->get_scheduler()->initialize(
      Model::get_config()->get_simulation_timeframe().get_starting_date(),
      Model::get_config()->get_simulation_timeframe().get_ending_date());
  p_model->get_population()->initialize();

  return result;
}

// ==================== EXISTING: default efficacy ====================
double getEfficacyForTherapy(std::string g_str, Model* p_model,
                             utils::Cli::DxGAppInput& input, int therapy_id) {
  Therapy* mainTherapy = p_model->get_therapy_db()[therapy_id].get();
  dynamic_cast<SFTStrategy*>(p_model->get_treatment_strategy())
      ->get_therapy_list()
      .clear();
  dynamic_cast<SFTStrategy*>(p_model->get_treatment_strategy())
      ->add_therapy(mainTherapy);

  p_model->get_reporters().clear();
  if (!input.output_file.empty()) {
    p_model->add_reporter(std::make_unique<PkPdReporter>(&input));
  } else {
    p_model->add_reporter(std::make_unique<PkPdReporter>());
  }

  for (auto& person :
       p_model->get_population()->all_persons()->v_person()) {
    auto density = Model::get_config()
                       ->get_parasite_parameters()
                       .get_parasite_density_levels()
                       .get_log_parasite_density_from_liver();
    auto* genotype = p_model->get_genotype_db()->get_genotype(g_str);
    auto* blood_parasite = person->add_new_parasite_to_blood(genotype);

    person->get_immune_system()->set_increase(true);
    person->set_host_state(Person::EXPOSED);

    blood_parasite->set_gametocyte_level(
        Model::get_config()
            ->get_epidemiological_parameters()
            .get_gametocyte_level_full());
    blood_parasite->set_last_update_log10_parasite_density(density);

    const int days_to_clinical =
        (person->get_age() <= 5)
            ? Model::get_config()
                  ->get_epidemiological_parameters()
                  .get_days_to_clinical_under_five()
            : Model::get_config()
                  ->get_epidemiological_parameters()
                  .get_days_to_clinical_over_five();
    auto event =
        std::make_unique<ProgressToClinicalEvent>(person.get());
    event->set_time(person->calculate_future_time(days_to_clinical));
    event->set_clinical_caused_parasite(blood_parasite);
    person->schedule_basic_event(std::move(event));
  }

  p_model->run();
  const auto result =
      1 - p_model->get_mdc()->blood_slide_prevalence_by_location()[0];

  p_model->set_population(std::make_unique<Population>());
  p_model->set_scheduler(std::make_unique<Scheduler>());
  p_model->get_scheduler()->initialize(
      Model::get_config()->get_simulation_timeframe().get_starting_date(),
      Model::get_config()->get_simulation_timeframe().get_ending_date());
  p_model->get_population()->initialize();

  return result;
}

// ==================== EXISTING: CRT calibration ====================
double getEfficacyForTherapyCRT(Model* p_model,
                                utils::Cli::DxGAppInput& input,
                                int therapy_id) {
  Therapy* mainTherapy = p_model->get_therapy_db()[therapy_id].get();
  dynamic_cast<SFTStrategy*>(p_model->get_treatment_coverage())
      ->get_therapy_list()
      .clear();
  dynamic_cast<SFTStrategy*>(p_model->get_treatment_coverage())
      ->add_therapy(mainTherapy);

  p_model->get_reporters().clear();
  if (!input.output_file.empty()) {
    p_model->add_reporter(std::make_unique<PkPdReporter>(&input));
  } else {
    p_model->add_reporter(std::make_unique<PkPdReporter>());
  }

  for (auto& person :
       p_model->get_population()->all_persons()->v_person()) {
    std::string g_str = "";
    int infect_prob = Model::get_random()->random_uniform(1, 104);
    if (infect_prob < 74) {
      g_str = input.genotypes[2];
    } else if (infect_prob < 91) {
      g_str = input.genotypes[1];
    } else {
      g_str = input.genotypes[0];
    }
    auto* genotype = p_model->get_genotype_db()->get_genotype(g_str);
    auto* blood_parasite = person->add_new_parasite_to_blood(genotype);
    auto density = Model::get_config()
                       ->get_parasite_parameters()
                       .get_parasite_density_levels()
                       .get_log_parasite_density_from_liver();

    person->get_immune_system()->set_increase(true);
    person->set_host_state(Person::EXPOSED);

    blood_parasite->set_gametocyte_level(
        Model::get_config()
            ->get_epidemiological_parameters()
            .get_gametocyte_level_full());
    blood_parasite->set_last_update_log10_parasite_density(density);

    const int days_to_clinical =
        (person->get_age() <= 5)
            ? Model::get_config()
                  ->get_epidemiological_parameters()
                  .get_days_to_clinical_under_five()
            : Model::get_config()
                  ->get_epidemiological_parameters()
                  .get_days_to_clinical_over_five();
    auto event =
        std::make_unique<ProgressToClinicalEvent>(person.get());
    event->set_time(person->calculate_future_time(days_to_clinical));
    event->set_clinical_caused_parasite(blood_parasite);
    person->schedule_basic_event(std::move(event));
  }

  p_model->run();
  const auto result =
      1 - p_model->get_mdc()->blood_slide_prevalence_by_location()[0];

  p_model->set_population(std::make_unique<Population>());
  p_model->set_scheduler(std::make_unique<Scheduler>());
  p_model->get_scheduler()->initialize(
      Model::get_config()->get_simulation_timeframe().get_starting_date(),
      Model::get_config()->get_simulation_timeframe().get_ending_date());
  p_model->get_population()->initialize();

  return result;
}

// ==================== EXISTING: EE validation ====================
bool validate_config_for_ee(utils::Cli::DxGAppInput& input) {
  input.number_of_drugs_in_combination = input.half_life.size();

  if (input.number_of_drugs_in_combination > 5) {
    std::cerr << "Error: Number of drugs in combination should not greater "
                 "than 5"
              << std::endl;
    return false;
  }
  if (input.k_max.size() != input.number_of_drugs_in_combination ||
      input.EC50.size() != input.number_of_drugs_in_combination ||
      input.slope.size() != input.number_of_drugs_in_combination ||
      input.dosing_days.size() != input.number_of_drugs_in_combination) {
    std::cerr << "Error: Wrong number of drugs in combination" << std::endl;
    return false;
  }
  for (auto k : input.k_max) {
    if (k >= 1 || k < 0) {
      std::cerr << "Error: k_max should be in range of (0,1]" << std::endl;
      return false;
    }
  }
  for (auto ec50 : input.EC50) {
    if (ec50 < 0) {
      std::cerr << "Error: EC50 should be greater than 0." << std::endl;
      return false;
    }
  }
  for (auto n : input.slope) {
    if (n < 0) {
      std::cerr << "Error: n should greater than 0." << std::endl;
      return false;
    }
  }
  for (auto dosing : input.dosing_days) {
    if (dosing < 0) {
      std::cerr << "Error: dosing should greater than 0." << std::endl;
      return false;
    }
  }
  if (input.mean_drug_absorption.empty()) {
    for (int i = 0; i < input.number_of_drugs_in_combination; ++i) {
      input.mean_drug_absorption.push_back(1.0);
    }
  } else if (input.mean_drug_absorption.size() !=
             input.number_of_drugs_in_combination) {
    std::cerr << "Error: Wrong number of drugs in combination" << std::endl;
    return false;
  }
  return true;
}