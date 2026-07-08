/**
 * PopulationEventBuilder.cpp
 *
 * Implement the population event builder factory, and most of the functions
 * specific to producing event objects. More complex functions may be found in
 * separate files under the EventBuilders directory.
 */
#include "PopulationEventBuilder.h"

#include <algorithm>
#include <cstdint>
#include <iostream>
#include <vector>

#include "AnnualBetaUpdateEvent.hxx"
#include "AnnualCoverageUpdateEvent.hxx"
#include "ChangeCirculationPercentEvent.hxx"
#include "ChangeInterruptedFeedingRateEvent.h"
#include "ChangeMutationMaskEvent.h"
#include "ChangeMutationProbabilityPerLocusEvent.h"
#include "ChangeTreatmentCoverageEvent.h"
#include "ChangeTreatmentStrategyEvent.h"
#include "ChangeWithinHostInducedFreeRecombinationEvent.h"
#include "Configuration/Config.h"
#include "DistrictImportationDailyEvent.h"
#include "ImportationEvent.h"
#include "ImportationPeriodicallyEvent.h"
#include "ImportationPeriodicallyRandomEvent.h"
#include "Introduce580YMutantEvent.h"
#include "IntroduceAmodiaquineMutantEvent.h"
#include "IntroduceLumefantrineMutantEvent.h"
#include "IntroduceMutantEvent.hxx"
#include "IntroduceMutantRasterEvent.hxx"
#include "IntroducePlas2CopyParasiteEvent.h"
#include "IntroduceTripleMutantToDPMEvent.h"
#include "ModifyNestedMFTEvent.h"
#include "Parasites/Genotype.h"
#include "RotateStrategyEvent.h"
#include "Simulation/Model.h"
#include "SingleRoundMDAEvent.h"
#include "TurnOffMutationEvent.h"
#include "TurnOnMutationEvent.h"
#include "UpdateBetaRasterEvent.hxx"

// Disable data flow analysis (DFA) in CLion
#pragma clang diagnostic push
#pragma ide diagnostic ignored "OCDFAInspection"

std::vector<std::unique_ptr<WorldEvent>> PopulationEventBuilder::build_introduce_parasite_events(
    const YAML::Node &node, Config* config) {
  std::vector<std::unique_ptr<WorldEvent>> events;
  for (const auto &entry : node) {
    auto location = entry["location"].as<int>();
    if (static_cast<std::size_t>(location) < config->number_of_locations()) {
      for (std::size_t j = 0; j < entry["parasite_info"].size(); j++) {
        auto genotype_aa_sequence =
            entry["parasite_info"][j]["genotype_aa_sequence"].as<std::string>();
        auto genotype_id =
            Model::get_genotype_db()->get_genotype(genotype_aa_sequence)->genotype_id();
        auto num = entry["parasite_info"][j]["number_of_cases"].as<int>();

        const auto starting_date = entry["parasite_info"][j]["date"].as<date::year_month_day>();
        auto time = (date::sys_days{starting_date}
                     - date::sys_days{config->get_simulation_timeframe().get_starting_date()})
                        .count();

        auto event = std::make_unique<ImportationEvent>(location, time, genotype_id, num);
        events.push_back(std::move(event));
      }
    }
  }
  return events;
}

std::vector<std::unique_ptr<WorldEvent>>
PopulationEventBuilder::build_introduce_parasites_periodically_events(const YAML::Node &node,
                                                                      Config* config) {
  std::vector<std::unique_ptr<WorldEvent>> events;

  for (const auto &entry : node) {
    const auto location = entry["location"].as<uint64_t>();
    const uint64_t location_from = location;
    const auto location_to =
        std::min(location + 1, static_cast<uint64_t>(config->number_of_locations()));

    for (auto loc = location_from; loc < location_to; ++loc) {
      for (std::size_t j = 0; j < entry["parasite_info"].size(); j++) {
        //            InitialParasiteInfo ipi;
        //            ipi.location = location;
        auto genotype_aa_sequence =
            entry["parasite_info"][j]["genotype_aa_sequence"].as<std::string>();
        auto genotype_id =
            Model::get_genotype_db()->get_genotype(genotype_aa_sequence)->genotype_id();
        // TODO: implement new importation parasite genotype based on allele
        // distribution

        const auto dur = entry["parasite_info"][j]["duration"].as<int>();
        const auto num = entry["parasite_info"][j]["number_of_cases"].as<int>();

        const auto starting_date =
            entry["parasite_info"][j]["start_date"].as<date::year_month_day>();
        auto time = (date::sys_days{starting_date}
                     - date::sys_days{config->get_simulation_timeframe().get_starting_date()})
                        .count();

        auto event = std::make_unique<ImportationPeriodicallyEvent>(static_cast<int>(loc), dur,
                                                                    genotype_id, num, time);
        events.push_back(std::move(event));
      }
    }
  }
  return events;
}

std::vector<std::unique_ptr<WorldEvent>>
PopulationEventBuilder::build_introduce_parasites_periodically_events_v2(const YAML::Node &node,
                                                                         Config* config) {
  std::vector<std::unique_ptr<WorldEvent>> events;
  for (const auto &event_node : node) {
    const auto location = event_node["location"].as<uint64_t>();
    const auto location_from = location == -1 ? 0 : location;
    const auto location_to =
        location == -1
            ? config->number_of_locations()
            : std::min(location + 1, static_cast<uint64_t>(config->number_of_locations()));

    for (auto loc = location_from; loc < location_to; ++loc) {
      for (auto j = 0; j < event_node["parasite_info"].size(); j++) {
        const auto dur = event_node["parasite_info"][j]["duration"].as<int>();
        const auto num = event_node["parasite_info"][j]["number_of_cases"].as<int>();

        const auto starting_date =
            event_node["parasite_info"][j]["start_date"].as<date::year_month_day>();

        date::year_month_day end_date =
            Model::get_config()->get_simulation_timeframe().get_ending_date();

        if (event_node["parasite_info"][j]["end_date"]) {
          end_date = event_node["parasite_info"][j]["end_date"].as<date::year_month_day>();
        }

        auto time = (date::sys_days{starting_date}
                     - date::sys_days{config->get_simulation_timeframe().get_starting_date()})
                        .count();
        auto end_time = (date::sys_days{end_date}
                         - date::sys_days{config->get_simulation_timeframe().get_starting_date()})
                            .count();
        throw std::runtime_error("Not implemented");
        // TODO: rework this with new genotype implementation
        // std::vector<std::vector<double>>
        // allele_distributions(Model::CONFIG->genotype_info().loci_vector.size());
        // // generate default distributions
        // for (int k = 0; k <
        // Model::CONFIG->genotype_info().loci_vector.size(); ++k) {
        //  auto number_of_alleles =
        //  Model::CONFIG->genotype_info().loci_vector[k].alleles.size(); for
        //  (int l = 0; l < number_of_alleles; ++l) {
        //    allele_distributions[k].push_back(1.0/number_of_alleles);
        //  }
        // }
        //
        // // read and replace
        // for (auto m = 0; m <
        // node[i]["parasite_info"][j]["allele_distributions"].size(); m++) {
        //  auto pos =
        //  node[i]["parasite_info"][j]["allele_distributions"][m]["position"].as<int>();
        //  for (int n = 0; n <
        //  node[i]["parasite_info"][j]["allele_distributions"][m]["distribution"].size();
        //  ++n) {
        //    allele_distributions[pos][n] =
        //    node[i]["parasite_info"][j]["allele_distributions"][m]["distribution"][n].as<double>();
        //  }
        // }
        //
        // auto* event = new
        // IntroduceParasitesPeriodicallyEventV2(allele_distributions, loc, dur,
        // num, time, end_time); events.push_back(event);
      }
    }
  }
  return events;
}

std::vector<std::unique_ptr<WorldEvent>>
PopulationEventBuilder::build_change_treatment_coverage_event(const YAML::Node &node,
                                                              Config* config) {
  std::vector<std::unique_ptr<WorldEvent>> events;
  for (const auto &entry : node) {
    auto tcm = ITreatmentCoverageModel::build(entry, config);
    auto event = std::make_unique<ChangeTreatmentCoverageEvent>(std::move(tcm));
    events.push_back(std::move(event));
  }
  return events;
}

std::vector<std::unique_ptr<WorldEvent>>
PopulationEventBuilder::build_change_treatment_strategy_event(const YAML::Node &node,
                                                              Config* config) {
  std::vector<std::unique_ptr<WorldEvent>> events;
  for (const auto &entry : node) {
    const auto starting_date = entry["date"].as<date::year_month_day>();
    auto time = (date::sys_days{starting_date}
                 - date::sys_days{config->get_simulation_timeframe().get_starting_date()})
                    .count();
    auto strategy_id = entry["strategy_id"].as<int>();

    // Verify that the strategy id is valid, if not fail
    if (strategy_id >= Model::get_strategy_db().size()) {
      spdlog::error("Invalid strategy_id! {} supplied, but strategy_db size is {}", strategy_id,
                    Model::get_strategy_db().size());
      exit(-1);
    }

    auto event = std::make_unique<ChangeTreatmentStrategyEvent>(strategy_id,time);
    events.push_back(std::move(event));
  }

  return events;
}

std::vector<std::unique_ptr<WorldEvent>> PopulationEventBuilder::build_single_round_mda_event(
    const YAML::Node &node, Config* config) {
  std::vector<std::unique_ptr<WorldEvent>> events;
  for (const auto &entry : node) {
    const auto starting_date = entry["date"].as<date::year_month_day>();
    auto time = (date::sys_days{starting_date}
                 - date::sys_days{config->get_simulation_timeframe().get_starting_date()})
                    .count();
    auto event = std::make_unique<SingleRoundMDAEvent>(time);
    for (std::size_t loc = 0; loc < config->number_of_locations(); loc++) {
      auto input_loc =
          entry["fraction_population_targeted"].size() < config->number_of_locations() ? 0 : loc;
      event->add_fraction_population_targeted(
          entry["fraction_population_targeted"][input_loc].as<double>());
    }

    event->set_days_to_complete(entry["days_to_complete_all_treatments"].as<int>());
    events.push_back(std::move(event));
  }

  return events;
}

std::vector<std::unique_ptr<WorldEvent>>
PopulationEventBuilder::build_modify_nested_mft_strategy_event(const YAML::Node &node,
                                                               Config* config) {
  std::vector<std::unique_ptr<WorldEvent>> events;
  for (const auto &entry : node) {
    const auto starting_date = entry["date"].as<date::year_month_day>();
    auto time = (date::sys_days{starting_date}
                 - date::sys_days{config->get_simulation_timeframe().get_starting_date()})
                    .count();
    auto strategy_id = entry["strategy_id"].as<int>();

    auto event = std::make_unique<ModifyNestedMFTEvent>(time, strategy_id);
    events.push_back(std::move(event));
  }

  return events;
}

std::vector<std::unique_ptr<WorldEvent>>
PopulationEventBuilder::build_introduce_plas2_parasite_events(const YAML::Node &node,
                                                              Config* config) {
  std::vector<std::unique_ptr<WorldEvent>> events;
  for (const auto &entry : node) {
    int location = entry["location"].as<int>();
    if (static_cast<std::size_t>(location) < config->number_of_locations()) {
      auto fraction = entry["fraction"].as<double>();

      const auto starting_date = entry["date"].as<date::year_month_day>();
      auto time = (date::sys_days{starting_date}
                   - date::sys_days{config->get_simulation_timeframe().get_starting_date()})
                      .count();

      auto event = std::make_unique<IntroducePlas2CopyParasiteEvent>(location, time, fraction);
      events.push_back(std::move(event));
    }
  }
  return events;
}

std::vector<std::unique_ptr<WorldEvent>> PopulationEventBuilder::build_turn_on_mutation_event(
    const YAML::Node &node, Config* config) {
  std::vector<std::unique_ptr<WorldEvent>> events;
  for (const auto &event_node : node) {
    const auto starting_date = event_node["date"].as<date::year_month_day>();
    auto time = (date::sys_days{starting_date}
                 - date::sys_days{config->get_simulation_timeframe().get_starting_date()})
                    .count();
    double mutation_probability =
        event_node["mutation_probability"]
            ? event_node["mutation_probability"].as<double>()
            : Model::get_config()->get_genotype_parameters().get_mutation_probability_per_locus();

    auto event = std::make_unique<TurnOnMutationEvent>(time, mutation_probability);
    events.push_back(std::move(event));
  }

  return events;
}

std::vector<std::unique_ptr<WorldEvent>> PopulationEventBuilder::build_turn_off_mutation_event(
    const YAML::Node &node, Config* config) {
  std::vector<std::unique_ptr<WorldEvent>> events;
  for (const auto &event_node : node) {
    const auto starting_date = event_node["date"].as<date::year_month_day>();
    auto time = (date::sys_days{starting_date}
                 - date::sys_days{config->get_simulation_timeframe().get_starting_date()})
                    .count();
    auto event = std::make_unique<TurnOffMutationEvent>(time);
    events.push_back(std::move(event));
  }

  return events;
}

std::vector<std::unique_ptr<WorldEvent>>
PopulationEventBuilder::build_change_interrupted_feeding_rate_event(const YAML::Node &node,
                                                                    Config* config) {
  std::vector<std::unique_ptr<WorldEvent>> events;
  for (const auto &event_node : node) {
    auto location = event_node["location"].as<int>();
    if (location < config->number_of_locations()) {
      const auto starting_date = event_node["date"].as<date::year_month_day>();
      auto time = (date::sys_days{starting_date}
                   - date::sys_days{config->get_simulation_timeframe().get_starting_date()})
                      .count();
      auto ifr = event_node["interrupted_feeding_rate"].as<double>();
      auto event = std::make_unique<ChangeInterruptedFeedingRateEvent>(location, ifr, time);
      events.push_back(std::move(event));
    }
  }
  return events;
}

std::vector<std::unique_ptr<WorldEvent>>
PopulationEventBuilder::build_change_within_host_induced_free_recombination_events(
    const YAML::Node node, Config* config) {
  std::vector<std::unique_ptr<WorldEvent>> events;
  for (const auto &event_node : node) {
    const auto starting_date = event_node["date"].as<date::year_month_day>();
    auto time = (date::sys_days{starting_date}
                 - date::sys_days{config->get_simulation_timeframe().get_starting_date()})
                    .count();
    auto value = event_node["value"].as<bool>();
    auto event = std::make_unique<ChangeWithinHostInducedFreeRecombinationEvent>(value, time);
    events.push_back(std::move(event));
  }

  return events;
}

std::vector<std::unique_ptr<WorldEvent>>
PopulationEventBuilder::build_change_mutation_probability_per_locus_events(const YAML::Node node,
                                                                           Config* config) {
  std::vector<std::unique_ptr<WorldEvent>> events;
  for (const auto &event_node : node) {
    const auto starting_date = event_node["date"].as<date::year_month_day>();
    auto time = (date::sys_days{starting_date}
                 - date::sys_days{config->get_simulation_timeframe().get_starting_date()})
                    .count();
    auto value = event_node["mutation_probability_per_locus"].as<double>();
    auto event = std::make_unique<ChangeMutationProbabilityPerLocusEvent>(value, time);
    events.push_back(std::move(event));
  }

  return events;
}

std::vector<std::unique_ptr<WorldEvent>>
PopulationEventBuilder::build_introduce_amodiaquine_mutant_parasite_events(const YAML::Node &node,
                                                                           Config* config) {
  std::vector<std::unique_ptr<WorldEvent>> events;
  for (const auto &entry : node) {
    auto location = entry["location"].as<int>();
    if (static_cast<std::size_t>(location) < config->number_of_locations()) {
      auto fraction = entry["fraction"].as<double>();

      const auto starting_date = entry["date"].as<date::year_month_day>();
      auto time = (date::sys_days{starting_date}
                   - date::sys_days{config->get_simulation_timeframe().get_starting_date()})
                      .count();
      std::vector<std::tuple<int, int, char>> alleles;
      for (const auto &allele_node : entry["alleles"]) {
        if (allele_node["allele"].as<std::string>().size() > 1) {
          spdlog::error("Allele {} should be 1 character", allele_node["allele"].as<std::string>());
        } else {
          alleles.emplace_back(allele_node["chromosome"].as<int>(), allele_node["locus"].as<int>(),
                               allele_node["allele"].as<std::string>().front());
        }
      }
      for (auto &allele : alleles) {
        spdlog::info("Mutation at {}:{} {}", std::get<0>(allele), std::get<1>(allele),
                     std::get<2>(allele));
      }
      auto event =
          std::make_unique<IntroduceAmodiaquineMutantEvent>(location, time, fraction, alleles);
      events.push_back(std::move(event));
    }
  }
  return events;
}

std::vector<std::unique_ptr<WorldEvent>>
PopulationEventBuilder::build_introduce_lumefantrine_mutant_parasite_events(const YAML::Node &node,
                                                                            Config* config) {
  std::vector<std::unique_ptr<WorldEvent>> events;
  for (const auto &entry : node) {
    int location = entry["location"].as<int>();
    if (static_cast<std::size_t>(location) < config->number_of_locations()) {
      auto fraction = entry["fraction"].as<double>();

      const auto starting_date = entry["date"].as<date::year_month_day>();
      auto time = (date::sys_days{starting_date}
                   - date::sys_days{config->get_simulation_timeframe().get_starting_date()})
                      .count();

      std::vector<std::tuple<int, int, char>> alleles;
      for (const auto &allele_node : entry["alleles"]) {
        if (allele_node["allele"].as<std::string>().size() > 1) {
          spdlog::error("Allele {} should be 1 character", allele_node["allele"].as<std::string>());
        } else {
          alleles.emplace_back(allele_node["chromosome"].as<int>(), allele_node["locus"].as<int>(),
                               allele_node["allele"].as<std::string>().front());
        }
      }
      for (auto &allele : alleles) {
        spdlog::info("Mutation at {}:{} {}", std::get<0>(allele), std::get<1>(allele),
                     std::get<2>(allele));
      }
      auto event =
          std::make_unique<IntroduceLumefantrineMutantEvent>(location, time, fraction, alleles);
      events.push_back(std::move(event));
    }
  }
  return events;
}

std::vector<std::unique_ptr<WorldEvent>> PopulationEventBuilder::build_introduce_580Y_mutant_events(
    const YAML::Node &node, Config* config) {
  std::vector<std::unique_ptr<WorldEvent>> events;
  for (const auto &event_node : node) {
    auto location = event_node["location"].as<int>();
    if (location < config->number_of_locations()) {
      auto fraction = event_node["fraction"].as<double>();

      const auto starting_date = event_node["date"].as<date::year_month_day>();
      auto time = (date::sys_days{starting_date}
                   - date::sys_days{config->get_simulation_timeframe().get_starting_date()})
                      .count();

      std::vector<std::tuple<int, int, char>> alleles;
      for (const auto &allele_node : event_node["alleles"]) {
        if (allele_node["allele"].as<std::string>().size() > 1) {
          spdlog::error("Allele {} should be 1 character", allele_node["allele"].as<std::string>());
        } else {
          alleles.emplace_back(allele_node["chromosome"].as<int>(), allele_node["locus"].as<int>(),
                               allele_node["allele"].as<std::string>().front());
        }
      }
      for (auto &allele : alleles) {
        spdlog::info("Mutation at {}:{} {}", std::get<0>(allele), std::get<1>(allele),
                     std::get<2>(allele));
      }
      auto event = std::make_unique<Introduce580YMutantEvent>(location, time, fraction, alleles);
      events.push_back(std::move(event));
    }
  }
  return events;
}

std::vector<std::unique_ptr<WorldEvent>>
PopulationEventBuilder::build_introduce_triple_mutant_to_dpm_parasite_events(const YAML::Node &node,
                                                                             Config* config) {
  std::vector<std::unique_ptr<WorldEvent>> events;
  for (const auto &event_node : node) {
    auto location = event_node["location"].as<int>();
    if (location < config->number_of_locations()) {
      auto fraction = event_node["fraction"].as<double>();

      const auto starting_date = event_node["date"].as<date::year_month_day>();
      auto time = (date::sys_days{starting_date}
                   - date::sys_days{config->get_simulation_timeframe().get_starting_date()})
                      .count();

      std::vector<std::tuple<int, int, char>> alleles;
      for (const auto &allele_node : event_node["alleles"]) {
        if (allele_node["allele"].as<std::string>().size() > 1) {
          spdlog::error("Allele {} should be 1 character", allele_node["allele"].as<std::string>());
        } else {
          alleles.emplace_back(allele_node["chromosome"].as<int>(), allele_node["locus"].as<int>(),
                               allele_node["allele"].as<std::string>().front());
        }
      }
      for (auto &allele : alleles) {
        spdlog::info("Mutation at {}:{} {}", std::get<0>(allele), std::get<1>(allele),
                     std::get<2>(allele));
      }
      auto event =
          std::make_unique<IntroduceTrippleMutantToDPMEvent>(location, time, fraction, alleles);
      events.push_back(std::move(event));
    }
  }
  return events;
}

std::vector<std::unique_ptr<WorldEvent>>
PopulationEventBuilder::build_rotate_treatment_strategy_event(const YAML::Node &node,
                                                              Config* config) {
  try {
    std::vector<std::unique_ptr<WorldEvent>> events;
    for (const auto &entry : node) {
      // Load the values
      auto start_date = entry["date"].as<date::year_month_day>();
      auto time = (date::sys_days{start_date}
                   - date::sys_days{config->get_simulation_timeframe().get_starting_date()})
                      .count();
      auto years = entry["years"].as<int>();
      auto first_strategy_id = entry["first_strategy_id"].as<int>();
      auto second_strategy_id = entry["second_strategy_id"].as<int>();

      // Make sure the years are reasonable
      if (years < 1) {
        spdlog::error(
            "Strategy rotation must be at least one year (whole numbers) in "
            "RotateStrategyEvent");
        throw std::invalid_argument("Strategy rotation years less than one");
      }

      // Check to make sure the strategy ids are valid
      if (first_strategy_id < 0 || second_strategy_id < 0) {
        spdlog::error(
            "Strategy id cannot be less than zero for "
            "RotateStrategyEvent");
        throw std::invalid_argument("Strategy id cannot be less than zero");
      }
      if (first_strategy_id >= Model::get_genotype_db()->size()
          || second_strategy_id >= Model::get_genotype_db()->size()) {
        spdlog::error(
            "Strategy id should be less than the total number of "
            "strategies (zero-indexing) for RotateStrategyEvent");
        throw std::invalid_argument("Strategy id greater than strategy_db size");
      }

      // Log and add the event to the queue
      auto event =
          std::make_unique<RotateStrategyEvent>(time, years, first_strategy_id, second_strategy_id);
      spdlog::debug(
          "Adding {} start: {}, rotation schedule: {}, initial strategy: {}, "
          "next strategy: {}",
          event->name(), StringHelpers::date_as_string(start_date), years, first_strategy_id,
          second_strategy_id);
      events.push_back(std::move(event));
    }
    return events;
  } catch (YAML::BadConversion &error) {
    spdlog::error(
        "Unrecoverable error parsing YAML value in "
        "{} node: {}",
        RotateStrategyEvent::EventName, error.msg);
    exit(EXIT_FAILURE);
  }
}

// Generate a new annual event that adjusts the beta at each location within the
// model, assumes that the YAML node contains a rate of change and start date.
std::vector<std::unique_ptr<WorldEvent>> PopulationEventBuilder::build_annual_beta_update_event(
    const YAML::Node &node, Config* config) {
  try {
    // Check the node size
    verify_single_node(node, AnnualBetaUpdateEvent::EventName);

    // Build the event
    auto start_date = node[0]["date"].as<date::year_month_day>();
    auto rate = node[0]["rate"].as<float>();
    auto time = (date::sys_days{start_date}
                 - date::sys_days{config->get_simulation_timeframe().get_starting_date()})
                    .count();
    auto event = std::make_unique<AnnualBetaUpdateEvent>(rate, time);

    // Log and add the event to the queue, only one for the country
    spdlog::debug("Adding {} start: {}, rate: {}", event->name(),
                  StringHelpers::date_as_string(start_date), rate);
    std::vector<std::unique_ptr<WorldEvent>> events;
    events.push_back(std::move(event));
    return events;
  } catch (YAML::BadConversion &error) {
    spdlog::error(
        "Unrecoverable error parsing YAML value in "
        "{} node: {}",
        AnnualBetaUpdateEvent::EventName, error.msg);
    exit(EXIT_FAILURE);
  }
}

// Generate a new annual event that adjusts the coverage at each location within
// the model, assumes that the YAML node contains a rate of change and start
// date.
std::vector<std::unique_ptr<WorldEvent>> PopulationEventBuilder::build_annual_coverage_update_event(
    const YAML::Node &node, Config* config) {
  try {
    // Check the node size
    verify_single_node(node, AnnualCoverageUpdateEvent::EventName);

    // Build the event
    auto start_date = node[0]["date"].as<date::year_month_day>();
    auto rate = node[0]["rate"].as<float>();
    auto time = (date::sys_days{start_date}
                 - date::sys_days{config->get_simulation_timeframe().get_starting_date()})
                    .count();
    auto event = std::make_unique<AnnualCoverageUpdateEvent>(rate, time);

    // Log and add the event to the queue, only one for the country
    spdlog::debug("Adding {} start: {}, rate: {}", event->name(),
                  StringHelpers::date_as_string(start_date), rate);
    std::vector<std::unique_ptr<WorldEvent>> events;
    events.push_back(std::move(event));
    return events;
  } catch (YAML::BadConversion &error) {
    spdlog::error(
        "Unrecoverable error parsing YAML value in "
        "{} node: {}",
        AnnualCoverageUpdateEvent::EventName, error.msg);
    exit(EXIT_FAILURE);
  }
}

std::vector<std::unique_ptr<WorldEvent>>
PopulationEventBuilder::build_change_circulation_percent_event(const YAML::Node &node,
                                                               Config* config) {
  try {
    std::vector<std::unique_ptr<WorldEvent>> events;
    for (const auto &entry : node) {
      // Load the values
      auto start_date = entry["date"].as<date::year_month_day>();
      auto time = (date::sys_days{start_date}
                   - date::sys_days{config->get_simulation_timeframe().get_starting_date()})
                      .count();
      auto rate = entry["circulation_percent"].as<float>();

      // Make sure the rate makes sense
      if (rate < 0.0) {
        spdlog::error(
            "The daily population circulation percentage must be "
            "greater than zero for {}",
            ChangeCirculationPercentEvent::EventName);
        throw std::invalid_argument("Population circulation percentage must be greater than zero");
      }
      if (rate > 1.0) {
        spdlog::error(
            "The daily population circulation percentage must be "
            "less than one (i.e., 100%) for ",
            ChangeCirculationPercentEvent::EventName);
        throw std::invalid_argument("Population circulation percentage must be less than one");
      }

      // Log and add the event to the queue
      auto event = std::make_unique<ChangeCirculationPercentEvent>(rate, time);
      spdlog::debug("Adding {} start: {}, rate: {}", event->name(),
                    StringHelpers::date_as_string(start_date), rate);
      events.push_back(std::move(event));
    }
    return events;
  } catch (YAML::BadConversion &error) {
    spdlog::error(
        "Unrecoverable error parsing YAML value in "
        "{} node: {}",
        ChangeCirculationPercentEvent::EventName, error.msg);
    exit(EXIT_FAILURE);
  }
}

// Generate a new importation periodically random event that uses a weighted
// random selection to add a new malaria infection with a specific genotype
// to the model.
std::vector<std::unique_ptr<WorldEvent>>
PopulationEventBuilder::build_importation_periodically_random_event(const YAML::Node &node,
                                                                    Config* config) {
  try {
    std::vector<std::unique_ptr<WorldEvent>> events;
    for (const auto &entry : node) {
      // Load the values
      auto start_date = entry["date"].as<date::year_month_day>();
      auto time = (date::sys_days{start_date}
                   - date::sys_days{config->get_simulation_timeframe().get_starting_date()})
                      .count();
      auto genotype_id = entry["genotype_id"].as<int>();
      auto count = entry["count"].as<int>();
      auto log_parasite_density = entry["log_parasite_density"].as<double>();

      // Check to make sure the date is valid
      if (start_date.day() != date::day{1}) {
        spdlog::error("The event must start on the first of the month for {} ",
                      ImportationPeriodicallyRandomEvent::EventName);
        throw std::invalid_argument("Event must start on the first of the month");
      }

      // Double check that the genotype id is valid
      if (genotype_id < 0) {
        spdlog::error(
            "Invalid genotype id supplied for {} genotype id cannot be less "
            "than zero",
            ImportationPeriodicallyRandomEvent::EventName);
        throw std::invalid_argument("Genotype id cannot be less than zero");
      }
      if (genotype_id >= Model::get_genotype_db()->size()) {
        spdlog::error(
            "Invalid genotype id supplied for {} genotype id cannot be greater "
            "than genotype_db size",
            ImportationPeriodicallyRandomEvent::EventName);
        throw std::invalid_argument("Genotype id cannot be greater than genotype_db size");
      }

      // Make sure the count makes sense
      if (count < 1) {
        spdlog::error(
            "The count of importations must be greater than zero for "
            "ImportationPeriodicallyRandomEvent");
        throw std::invalid_argument("Count must be greater than zero");
      }

      // Make sure the log parasite density makes sense
      if (log_parasite_density == 0) {
        spdlog::error(
            "Log parasite density cannot be zero for "
            "ImportationPeriodicallyRandomEvent Log10 of zero is undefined.");
        throw std::invalid_argument("Log10 of zero is undefined.");
      }

      // Log and add the event to the queue
      auto event = std::make_unique<ImportationPeriodicallyRandomEvent>(genotype_id, time, count,
                                                                        log_parasite_density);
      spdlog::debug(
          "Adding {} start: {}, genotype_id: {}, count: {}, "
          "log_parasite_density: {}",
          event->name(), StringHelpers::date_as_string(start_date), genotype_id, count,
          log_parasite_density);
      events.push_back(std::move(event));
    }
    return events;
  } catch (YAML::BadConversion &error) {
    spdlog::error(
        "Unrecoverable error parsing YAML value in "
        "{} node: {}",
        ImportationPeriodicallyRandomEvent::EventName, error.msg);
    exit(EXIT_FAILURE);
  }
}

std::vector<std::unique_ptr<WorldEvent>> PopulationEventBuilder::build_update_beta_raster_event(
    const YAML::Node &node, Config* config) {
  try {
    std::vector<std::unique_ptr<WorldEvent>> events;
    for (const auto &entry : node) {
      // Load the values
      auto start_date = entry["date"].as<date::year_month_day>();
      auto time = (date::sys_days{start_date}
                   - date::sys_days{config->get_simulation_timeframe().get_starting_date()})
                      .count();
      auto filename = entry["beta_raster"].as<std::string>();

      // Make sure the file actually exists
      std::ifstream file;
      file.open(filename);
      if (!file) {
        spdlog::error(
            "The file indicated, {}, cannot be found for "
            "{} event. Please check the file path.",
            filename, UpdateBetaRasterEvent::EVENT_NAME);
        throw std::invalid_argument("File for " + UpdateBetaRasterEvent::EVENT_NAME
                                    + " does not appear to exist.");
      }
      file.close();

      // Log and add the event to the queue
      auto event = std::make_unique<UpdateBetaRasterEvent>(filename, time);
      spdlog::debug("Adding {} start: {}, filename: {}", event->name(),
                    StringHelpers::date_as_string(start_date), filename);
      events.push_back(std::move(event));
    }
    return events;
  } catch (YAML::BadConversion &error) {
    spdlog::error(
        "Unrecoverable error parsing YAML value in "
        "{} node: {}",
        UpdateBetaRasterEvent::EVENT_NAME, error.msg);
    exit(EXIT_FAILURE);
  }
}

std::vector<std::unique_ptr<WorldEvent>>
PopulationEventBuilder::build_import_district_mutant_daily_events(const YAML::Node &node,
                                                                  Config* config) {
  std::vector<std::unique_ptr<WorldEvent>> events;
  for (const auto &entry : node) {
    auto district = entry["district"].as<int>();
    auto daily_rate = entry["daily_rate"].as<double>();
    std::vector<std::tuple<int, int, char>> alleles;
    for (const auto &allele_node : entry["alleles"]) {
      if (allele_node["allele"].as<std::string>().size() > 1) {
        spdlog::error("Allele {} should be 1 character", allele_node["allele"].as<std::string>());
      } else {
        alleles.emplace_back(allele_node["chromosome"].as<int>(), allele_node["locus"].as<int>(),
                             allele_node["allele"].as<std::string>().front());
      }
    }
    for (auto &allele : alleles) {
      spdlog::info("Mutation at {}:{} {}", std::get<0>(allele), std::get<1>(allele),
                   std::get<2>(allele));
    }
    auto start_date = entry["start_date"].as<date::year_month_day>();
    auto start_day = (date::sys_days{start_date}
                      - date::sys_days{config->get_simulation_timeframe().get_starting_date()})
                         .count();
    auto event =
        std::make_unique<DistrictImportationDailyEvent>(district, daily_rate, start_day, alleles);
    events.push_back(std::move(event));
  }
  return events;
}

std::vector<std::unique_ptr<WorldEvent>> PopulationEventBuilder::build_change_mutation_mask_events(
    const YAML::Node &node, Config* config) {
  std::vector<std::unique_ptr<WorldEvent>> events;
  for (const auto &event_node : node) {
    const auto starting_date = event_node["date"].as<date::year_month_day>();
    auto time = (date::sys_days{starting_date}
                 - date::sys_days{config->get_simulation_timeframe().get_starting_date()})
                    .count();
    auto mutation_mask = event_node["mutation_mask"].as<std::vector<bool>>();

    auto event = std::make_unique<ChangeMutationMaskEvent>(mutation_mask, time);
    events.push_back(std::move(event));
  }

  return events;
}

std::vector<std::unique_ptr<WorldEvent>> PopulationEventBuilder::build(const YAML::Node &node) {
  Config* config = Model::get_config();
  std::vector<std::unique_ptr<WorldEvent>> events;
  const auto name = node["name"].as<std::string>();
  spdlog::info("Building events of type: {}", name);

  if (name == "introduce_parasites") {
    events = build_introduce_parasite_events(node["info"], config);
  }
  if (name == "introduce_parasites_periodically") {
    events = build_introduce_parasites_periodically_events(node["info"], config);
  }
  if (name == "introduce_parasites_periodically_v2") {
    events = build_introduce_parasites_periodically_events_v2(node["info"], config);
  }
  if (name == "change_treatment_coverage") {
    events = build_change_treatment_coverage_event(node["info"], config);
  }
  if (name == "change_treatment_strategy") {
    events = build_change_treatment_strategy_event(node["info"], config);
  }
  if (name == "single_round_MDA") { events = build_single_round_mda_event(node["info"], config); }
  if (name == "modify_nested_mft_strategy") {
    events = build_modify_nested_mft_strategy_event(node["info"], config);
  }
  if (name == "introduce_plas2_parasites") {
    events = build_introduce_plas2_parasite_events(node["info"], config);
  }
  if (name == "introduce_amodiaquine_mutant_parasites") {
    events = build_introduce_amodiaquine_mutant_parasite_events(node["info"], config);
  }
  if (name == "introduce_lumefantrine_mutant_parasites") {
    events = build_introduce_lumefantrine_mutant_parasite_events(node["info"], config);
  }
  if (name == "introduce_580Y_parasites") {
    events = build_introduce_580Y_mutant_events(node["info"], config);
  }
  if (name == "turn_on_mutation") { events = build_turn_on_mutation_event(node["info"], config); }
  if (name == "turn_off_mutation") { events = build_turn_off_mutation_event(node["info"], config); }
  if (name == "change_within_host_induced_free_recombination") {
    events = build_change_within_host_induced_free_recombination_events(node["info"], config);
  }
  if (name == "change_mutation_probability_per_locus") {
    events = build_change_mutation_probability_per_locus_events(node["info"], config);
  }
  if (name == "change_interrupted_feeding_rate") {
    events = build_change_interrupted_feeding_rate_event(node["info"], config);
  }
  if (name == "introduce_triple_mutant_to_dpm_parasites") {
    events = build_introduce_triple_mutant_to_dpm_parasite_events(node["info"], config);
  }

  if (name == AnnualBetaUpdateEvent::EventName) {
    events = build_annual_beta_update_event(node["info"], config);
  }
  if (name == AnnualCoverageUpdateEvent::EventName) {
    events = build_annual_coverage_update_event(node["info"], config);
  }
  if (name == ChangeCirculationPercentEvent::EventName) {
    events = build_change_circulation_percent_event(node["info"], config);
  }
  if (name == ImportationPeriodicallyRandomEvent::EventName) {
    events = build_importation_periodically_random_event(node["info"], config);
  }
  if (name == IntroduceMutantEvent::EVENT_NAME) {
    auto admin_level_name = node["admin_level"].as<std::string>();
    events = build_introduce_mutant_event(node["info"], config, admin_level_name);
  }
  if (name == IntroduceMutantRasterEvent::EventName) {
    events = build_introduce_mutant_raster_event(node["info"], config);
  }
  if (name == UpdateBetaRasterEvent::EVENT_NAME) {
    events = build_update_beta_raster_event(node["info"], config);
  }
  if (name == RotateStrategyEvent::EventName) {
    events = build_rotate_treatment_strategy_event(node["info"], config);
  }
  if (name == DistrictImportationDailyEvent::EventName) {
    events = build_import_district_mutant_daily_events(node["info"], config);
  }
  if (name == "change_mutation_mask") {
    events = build_change_mutation_mask_events(node["info"], config);
  }
  return events;
}

void PopulationEventBuilder::verify_single_node(const YAML::Node &node, const std::string &name) {
  if (node.size() > 1) {
    spdlog::error("More than one sub node found for " + name + " in the configuration file");
  }
}

