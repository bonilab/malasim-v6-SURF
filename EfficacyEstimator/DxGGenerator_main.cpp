/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/*
 * File:   main.cpp
 * Author: Merlin
 *
 * Created on January 12, 2017, 4:31 PM
 */
// TODO: make it works, input for genotype will be string or a file with a list of string

#include <spdlog/spdlog.h>

#include <cmath>  // log10
#include <cstdlib>
#include <iomanip>
#include <iostream>

#include "Configuration/Config.h"
#include "Events/ProgressToClinicalEvent.h"
#include "MDC/ModelDataCollector.h"
#include "Parasites/Genotype.h"
#include "PkPdReporter.h"
#include "Population/ImmuneSystem/ImmuneSystem.h"
#include "Population/Person/Person.h"
#include "Population/Population.h"
#include "Simulation/Model.h"
#include "Treatment/Strategies/IStrategy.h"
#include "Treatment/Strategies/SFTStrategy.h"
#include "Treatment/Therapies//SCTherapy.h"
#include "Treatment/Therapies/Therapy.h"
#include "Utils/Index/PersonIndexAll.h"
#include "Utils/Random.h"

bool validate_config_for_ee(utils::DxGAppInput &input);
double get_efficacy_for_therapy(std::string g_str, Model* p_model, utils::DxGAppInput &input,
                                int therapy_id);
double get_efficacy_for_therapy_crt(Model* p_model, utils::DxGAppInput &input, int therapy_id);

// efficacy_map efficacies;

inline double round(double val) {
  if (val < 0) return ceil(val - 0.5);
  return floor(val + 0.5);
}

/*
 *
 */
int main(int argc, char** argv) {
  auto input = utils::Cli::parse_dxg_args(argc, argv);
  Model::get_instance()->initialize();

  auto* p_model = Model::get_instance();

  if (input.as_iiv != -1) {
    for (auto &sd : Model::get_drug_db()->at(0)->age_group_specific_drug_concentration_sd()) {
      sd = input.as_iiv;
    }
  }

  if (input.as_ec50 != -1) {
    // TODO: fix it
    //    p_model->get_config()->EC50_power_n_table()[0][0] = pow(as_ec50,
    //    p_model->get_config()->drug_db()->at(0)->n());
  }
  std::cout << std::setprecision(5);
  int max_therapy_id{0};
  int min_therapy_id{0};

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

  // TODO: Genotype should be imported  from input files

  if (input.is_crt_calibration) {
    if (input.genotypes.empty()) {
      std::cout << "List of population genotypes is empty" << '\n';
      exit(0);
    }
    std::stringstream ss;
    if (input.therapy_list.empty()) {
      for (auto therapy_id = min_therapy_id; therapy_id <= max_therapy_id; therapy_id++) {
        std::cout << *Model::get_therapy_db()[therapy_id] << "\t";
      }
    } else {
      for (auto therapy_id : input.therapy_list) {
        std::cout << *Model::get_therapy_db()[therapy_id] << "\t";
      }
    }
    std::cout << '\n';
    if (input.therapy_list.empty()) {
      for (auto therapy_id = min_therapy_id; therapy_id <= max_therapy_id; therapy_id++) {
        double efficacy = get_efficacy_for_therapy_crt(p_model, input, therapy_id);
        ss << efficacy << (therapy_id == max_therapy_id ? "" : "\t");
      }
    } else {
      for (int t_index = 0; t_index < input.therapy_list.size(); t_index++) {
        double efficacy = get_efficacy_for_therapy_crt(p_model, input, input.therapy_list[t_index]);
        ss << efficacy
           << (input.therapy_list[t_index] == input.therapy_list.size() - 1 ? "" : "\t");
      }
    }
    std::cout << ss.str() << '\n';
  } else if (input.is_ee_calibration) {
    if (!validate_config_for_ee(input)) {
      std::cout << "Parameters for Efficacy Estimator are not correct" << '\n';
      exit(0);
    } else if (input.genotypes.size() > 1) {
      std::cout << "Only 1 genotype is accepted using Efficacy Estimator" << '\n';
      exit(0);
    } else {
      // ==== override population size ======
      if (input.population_size != Model::get_population()->size()) {
        Model::get_config()->location_db()[0].population_size = input.population_size;
        Model::set_population(std::make_unique<Population>());
        Model::get_population()->initialize();
      }
      // ==== override drug type info ========
      auto start_drug_id = input.is_art ? 0 : 1;
      for (int i = 0; i < input.number_of_drugs_in_combination; i++) {
        auto* dt = Model::get_drug_db()->at(i + start_drug_id).get();
        dt->set_name(fmt::format("D{}", i));
        dt->set_drug_half_life(input.half_life[i]);
        dt->set_maximum_parasite_killing_rate(input.k_max[i]);
        dt->set_n(input.slope[i]);
        // TODO: add app arguments later
        //    dt->set_p_mutation(0.0);
        dt->set_k(4);
        for (double &mda : dt->age_specific_drug_absorption()) {
          mda = input.mean_drug_absorption[i];
        }
        //    Model::CONFIG->EC50_power_n_table()[0][i + start_drug_id] = pow(input.EC50[i],
        //    dt->n());
      }

      // ======= override therapy 0 ==========
      auto* sc_therapy = dynamic_cast<SCTherapy*>(Model::get_therapy_db()[0].get());
      sc_therapy->drug_ids.clear();
      sc_therapy->dosing_day.clear();

      for (int i = 0; i < input.number_of_drugs_in_combination; i++) {
        sc_therapy->drug_ids.push_back(i + start_drug_id);
        sc_therapy->dosing_day.push_back(input.dosing_days[i]);
      }

      // ==========reset and override reporters ==================
      p_model->get_reporters().clear();
      p_model->add_reporter(std::make_unique<PkPdReporter>(&input));

      Model::get_genotype_db()->clear();
      std::vector<Genotype*> genotype_inputs;
      for (const auto &genotype_str : input.genotypes) {
        genotype_inputs.push_back(Model::get_genotype_db()->get_genotype(genotype_str));
      }

      // =========infect population with genotype 0================
      auto* genotype =
          Model::get_genotype_db()->get_genotype(genotype_inputs.front()->get_aa_sequence());

      for (auto &person : Model::get_population()->all_persons()->v_person()) {
        auto density = Model::get_config()
                           ->get_parasite_parameters()
                           .get_parasite_density_levels()
                           .get_log_parasite_density_from_liver();
        auto* blood_parasite = person->add_new_parasite_to_blood(genotype);

        person->get_immune_system()->set_increase(true);
        person->set_host_state(Person::EXPOSED);

        blood_parasite->set_gametocyte_level(
            Model::get_config()->get_epidemiological_parameters().get_gametocyte_level_full());
        blood_parasite->set_last_update_log10_parasite_density(density);

        const int days_to_clinical = (person->get_age() <= 5)
                                         ? Model::get_config()
                                               ->get_epidemiological_parameters()
                                               .get_days_to_clinical_under_five()
                                         : Model::get_config()
                                               ->get_epidemiological_parameters()
                                               .get_days_to_clinical_over_five();
        auto event = std::make_unique<ProgressToClinicalEvent>(person.get());
        event->set_time(Person::calculate_future_time(days_to_clinical));
        event->set_clinical_caused_parasite(blood_parasite);
        person->schedule_basic_event(std::move(event));
      }

      // run model
      p_model->run();

      const auto result = 1 - Model::get_mdc()->blood_slide_prevalence_by_location()[0];
      //            fmt::print(
      //                "pop\tdose\thalflife\tkmax\tec50\tslope\tis_art\tefficacy\n"
      //                "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{:f}\n",
      //                p_model->get_population()->all_persons()->size(),
      //                fmt::join(input.dosing_days, "\t"), fmt::join(input.half_life, "\t"),
      //                fmt::join(input.k_max, "\t"), fmt::join(input.EC50, "\t"),
      //                fmt::join(input.slope, "\t"), input.is_art ? 1 : 0, result
      //            );
    }
  } else {
    std::cout << "ID\tGenotype\t";
    if (input.therapy_list.empty()) {
      for (auto therapy_id = min_therapy_id; therapy_id <= max_therapy_id; therapy_id++) {
        std::cout << *Model::get_therapy_db()[therapy_id] << "\t";
      }
    } else {
      for (auto therapy_id : input.therapy_list) {
        std::cout << *Model::get_therapy_db()[therapy_id] << "\t";
      }
    }
    std::cout << '\n';
    for (int g_index = 0; g_index < input.genotypes.size(); g_index++) {
      std::stringstream ss;
      if (input.is_old_format) {
        ss << g_index << "\t"
           << Model::get_mosquito()->get_old_genotype_string2(input.genotypes[g_index]) << "\t";
      } else {
        ss << g_index << "\t"
           << Model::get_mosquito()->get_old_genotype_string(input.genotypes[g_index]) << "\t";
      }
      if (input.therapy_list.empty()) {
        for (auto therapy_id = min_therapy_id; therapy_id <= max_therapy_id; therapy_id++) {
          double efficacy =
              get_efficacy_for_therapy(input.genotypes[g_index], p_model, input, therapy_id);
          //                ss << efficacy << (therapy_id == max_therapy_id ? "" : "\t");
          ss << efficacy << "\t";
        }
      } else {
        for (int t_index = 0; t_index < input.therapy_list.size(); t_index++) {
          double efficacy = get_efficacy_for_therapy(input.genotypes[g_index], p_model, input,
                                                     input.therapy_list[t_index]);
          //                ss << efficacy << (input.therapy_list[t_index] ==
          //                input.therapy_list.size() - 1 ? "" : "\t");
          ss << efficacy << "\t";
        }
      }
      std::cout << ss.str() << '\n';
    }
  }
  return 0;
}

double get_efficacy_for_therapy(std::string g_str, Model* p_model, utils::DxGAppInput &input,
                                int therapy_id) {
  Therapy* main_therapy = Model::get_therapy_db()[therapy_id].get();
  dynamic_cast<SFTStrategy*>(Model::get_treatment_strategy())->get_therapy_list().clear();
  dynamic_cast<SFTStrategy*>(Model::get_treatment_strategy())->add_therapy(main_therapy);

  // reset reporter
  p_model->get_reporters().clear();

  if (!input.output_file.empty()) {
    p_model->add_reporter(std::make_unique<PkPdReporter>(&input));
  } else {
    p_model->add_reporter(std::make_unique<PkPdReporter>());
  }

  for (auto &person : Model::get_population()->all_persons()->v_person()) {
    auto density = Model::get_config()
                       ->get_parasite_parameters()
                       .get_parasite_density_levels()
                       .get_log_parasite_density_from_liver();
    auto* genotype = Model::get_genotype_db()->get_genotype(g_str);
    auto* blood_parasite = person->add_new_parasite_to_blood(genotype);

    person->get_immune_system()->set_increase(true);
    person->set_host_state(Person::EXPOSED);

    blood_parasite->set_gametocyte_level(
        Model::get_config()->get_epidemiological_parameters().get_gametocyte_level_full());
    blood_parasite->set_last_update_log10_parasite_density(density);

    const int days_to_clinical = (person->get_age() <= 5) ? Model::get_config()
                                                                ->get_epidemiological_parameters()
                                                                .get_days_to_clinical_under_five()
                                                          : Model::get_config()
                                                                ->get_epidemiological_parameters()
                                                                .get_days_to_clinical_over_five();
    auto event = std::make_unique<ProgressToClinicalEvent>(person.get());
    event->set_time(Person::calculate_future_time(days_to_clinical));
    event->set_clinical_caused_parasite(blood_parasite);
    person->schedule_basic_event(std::move(event));
  }

  p_model->run();
  const auto result = 1 - Model::get_mdc()->blood_slide_prevalence_by_location()[0];

  Model::set_population(std::make_unique<Population>());
  Model::set_scheduler(std::make_unique<Scheduler>());

  Model::get_scheduler()->initialize(
      Model::get_config()->get_simulation_timeframe().get_starting_date(),
      Model::get_config()->get_simulation_timeframe().get_ending_date());
  Model::get_population()->initialize();

  return result;
}

double get_efficacy_for_therapy_crt(Model* p_model, utils::DxGAppInput &input,
                                    int therapy_id) {
  Therapy* main_therapy = Model::get_therapy_db()[therapy_id].get();
  dynamic_cast<SFTStrategy*>(Model::get_treatment_coverage())->get_therapy_list().clear();
  dynamic_cast<SFTStrategy*>(Model::get_treatment_coverage())->add_therapy(main_therapy);

  // reset reporter
  p_model->get_reporters().clear();

  if (!input.output_file.empty()) {
    p_model->add_reporter(std::make_unique<PkPdReporter>(&input));
  } else {
    p_model->add_reporter(std::make_unique<PkPdReporter>());
  }

  for (auto &person : Model::get_population()->all_persons()->v_person()) {
    // The genotype distribution is from table 2 in http://dx.doi.org/10.1016/S1473-3099(19)30391-3
    // Run these 3 genotypes in population with and without F145I,T93S,H97Y and I218F
    // to get ec50 of WT and mutant genotypes
    std::string g_str;
    int infect_prob = Model::get_random()->random_uniform(1, 104);
    if (infect_prob < 74) {  // KEL1/PLA1
      g_str = input.genotypes[2];
    } else if (infect_prob < 91) {  // KEL1
      g_str = input.genotypes[1];
    } else {  // WT
      g_str = input.genotypes[0];
    }
    auto* genotype = Model::get_genotype_db()->get_genotype(g_str);
    auto* blood_parasite = person->add_new_parasite_to_blood(genotype);
    auto density = Model::get_config()
                       ->get_parasite_parameters()
                       .get_parasite_density_levels()
                       .get_log_parasite_density_from_liver();

    person->get_immune_system()->set_increase(true);
    person->set_host_state(Person::EXPOSED);

    blood_parasite->set_gametocyte_level(
        Model::get_config()->get_epidemiological_parameters().get_gametocyte_level_full());
    blood_parasite->set_last_update_log10_parasite_density(density);

    const int days_to_clinical = (person->get_age() <= 5) ? Model::get_config()
                                                                ->get_epidemiological_parameters()
                                                                .get_days_to_clinical_under_five()
                                                          : Model::get_config()
                                                                ->get_epidemiological_parameters()
                                                                .get_days_to_clinical_over_five();
    auto event = std::make_unique<ProgressToClinicalEvent>(person.get());
    event->set_time(Person::calculate_future_time(days_to_clinical));
    event->set_clinical_caused_parasite(blood_parasite);
    person->schedule_basic_event(std::move(event));
  }

  p_model->run();
  const auto result = 1 - Model::get_mdc()->blood_slide_prevalence_by_location()[0];

  Model::set_population(std::make_unique<Population>());
  Model::set_scheduler(std::make_unique<Scheduler>());

  Model::get_scheduler()->initialize(
      Model::get_config()->get_simulation_timeframe().get_starting_date(),
      Model::get_config()->get_simulation_timeframe().get_ending_date());
  Model::get_population()->initialize();

  return result;
}

bool validate_config_for_ee(utils::DxGAppInput &input) {
  input.number_of_drugs_in_combination = input.half_life.size();

  if (input.number_of_drugs_in_combination > 5) {
    std::cerr << "Error: Number of drugs in combination should not greater than 5" << '\n';
    return false;
  }

  if (input.k_max.size() != input.number_of_drugs_in_combination
      || input.ec50.size() != input.number_of_drugs_in_combination
      || input.slope.size() != input.number_of_drugs_in_combination
      || input.dosing_days.size() != input.number_of_drugs_in_combination) {
    std::cerr << "Error: Wrong number of drugs in combination" << '\n';
    return false;
  }

  for (auto k_m : input.k_max) {
    if (k_m >= 1 || k_m < 0) {
      std::cerr << "Error: k_max should be in range of (0,1]" << '\n';
      return false;
    }
  }

  for (auto ec50 : input.ec50) {
    if (ec50 < 0) {
      std::cerr << "Error: EC50 should be greater than 0." << '\n';
      return false;
    }
  }

  for (auto n : input.slope) {
    if (n < 0) {
      std::cerr << "Error: n should greater than 0." << '\n';
      return false;
    }
  }

  for (auto dosing : input.dosing_days) {
    if (dosing < 0) {
      std::cerr << "Error: dosing should greater than 0." << '\n';
      return false;
    }
  }

  if (input.mean_drug_absorption.empty()) {
    for (int i = 0; i < input.number_of_drugs_in_combination; ++i) {
      input.mean_drug_absorption.push_back(1.0);
    }
  } else if (input.mean_drug_absorption.size() != input.number_of_drugs_in_combination) {
    std::cerr << "Error: Wrong number of drugs in combination" << '\n';
    return false;
  }

  return true;
}

