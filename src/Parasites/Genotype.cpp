#include "Genotype.h"

#include <algorithm>

#include "Configuration/Config.h"
#include "Core/Scheduler/Scheduler.h"
#include "Population/DrugsInBlood.h"
#include "Simulation/Model.h"
#include "Treatment/Therapies/DrugDatabase.h"
#include "Utils/Helpers/NumberHelpers.h"

namespace {

bool drug_selects_for_double_copy(const GenotypeParameters::GeneInfo &gene_info,
                                  DrugsInBlood* drugs_in_blood) {
  if (drugs_in_blood == nullptr || drugs_in_blood->size() == 0) { return false; }

  for (const auto &[drug_id, drug] : *drugs_in_blood) {
    if (drug == nullptr || drug->last_update_value() <= 0) { continue; }

    for (const auto &cnv_effect : gene_info.get_cnv_multiplicative_effect_on_EC50()) {
      if (cnv_effect.get_drug_id() != drug_id) { continue; }
      const auto &factors = cnv_effect.get_factors();
      if (factors.size() > 1 && factors[1] > factors[0] && factors[1] > 1.0) { return true; }
    }
  }

  return false;
}
}  // namespace

Genotype::Genotype(const std::string &in_aa_sequence) : aa_sequence{in_aa_sequence} {
  // create aa structure
  std::string chromosome_str;
  std::istringstream token_stream(in_aa_sequence);
  auto ii = 0;
  while (std::getline(token_stream, chromosome_str, '|')) {
    std::string gene_str;
    std::istringstream chromosome_token_stream(chromosome_str);
    auto jj = 0;
    while (std::getline(chromosome_token_stream, gene_str, ',')) {
      pf_genotype_str[ii].push_back(gene_str);
      jj++;
    }
    ii++;
  }
  resistant_recombinations_in_mosquito = std::vector<MosquitoRecombinedGenotypeInfo>();
}

Genotype::~Genotype() = default;

bool Genotype::resist_to(DrugType* dt) {
  return EC50_power_n[dt->id()] > pow(dt->base_EC50, dt->n());
}

Genotype* Genotype::combine_mutation_to(const int &locus, const int &value) {
  // TODO: remove this
  return this;
}

Genotype* Genotype::modify_genotype_allele(const std::vector<std::tuple<int, int, char>> &alleles,
                                           Config* p_config) const {
  std::string new_aa_sequence{aa_sequence};
  /*
   * allele_map_info is (chromosome_id,locus_pos,allele)
   */
  // for (auto & allele : alleles) {
  //   spdlog::info("Mutation at {}:{}
  //   {}",std::get<0>(allele),std::get<1>(allele),std::get<2>(allele));
  // }
  for (auto const allele_info : alleles) {
    int chromosome_counter = 0;
    int locus_counter = 0;
    for (int char_counter = 0; char_counter < new_aa_sequence.size(); char_counter++) {
      if (chromosome_counter == (std::get<0>(allele_info) - 1)) {
        if (locus_counter == (std::get<1>(allele_info) - 1)) {
          if (p_config->get_genotype_parameters().get_mutation_mask()[char_counter] == '1') {
            if (new_aa_sequence[char_counter] != std::get<2>(allele_info)) {
              // spdlog::trace("{}:{} Changing allele at {} ({} -> {}) new aa_sequence {} mask: {}",
              //   std::get<0>(allele_info),std::get<1>(allele_info),
              //   char_counter,
              //   new_aa_sequence[char_counter],std::get<2>(allele_info),new_aa_sequence,
              //   pConfig->get_genotype_parameters().get_mutation_mask()[char_counter]);
              new_aa_sequence[char_counter] = std::get<2>(allele_info);
            }
          } else {
            spdlog::error(
                "{}:{} Changing allele at {} ({} -> {}) of genotype aa_sequence {} but the mask is "
                "0",
                std::get<0>(allele_info), std::get<1>(allele_info), char_counter,
                new_aa_sequence[char_counter], std::get<2>(allele_info), new_aa_sequence);
          }
        }
        locus_counter++;
      }
      if (new_aa_sequence[char_counter] == '|') chromosome_counter++;
    }
  }
  return Model::get_genotype_db()->get_genotype(new_aa_sequence);
}

double Genotype::get_EC50_power_n(DrugType* dt) { return EC50_power_n[dt->id()]; }

std::ostream &operator<<(std::ostream &os, Genotype &genotype) {
  os << genotype.genotype_id_ << "\t";
  os << genotype.get_aa_sequence();
  return os;
}

std::string Genotype::get_aa_sequence() const { return aa_sequence; }

bool Genotype::is_valid(const GenotypeParameters::PfGenotypeInfo &gene_info) {
  for (int chromosome_i = 0; chromosome_i < 14; ++chromosome_i) {
    auto chromosome_info = gene_info.chromosome_infos[chromosome_i];
    if (chromosome_info.get_genes().size() != pf_genotype_str[chromosome_i].size()) {
      spdlog::error("Error {} {}", pf_genotype_str[chromosome_i].size(),
                    chromosome_info.get_genes().size());
      return false;
    }

    for (int gene_i = 0; gene_i < pf_genotype_str[chromosome_i].size(); ++gene_i) {
      auto gene_info = chromosome_info.get_genes()[gene_i];
      auto max_aa_pos = gene_info.get_max_copies() > 1
                            ? pf_genotype_str[chromosome_i][gene_i].size() - 1
                            : pf_genotype_str[chromosome_i][gene_i].size();

      // check same size with aa postions info
      if (gene_info.get_aa_positions().size() != max_aa_pos) {
        spdlog::error("Error {} {}", pf_genotype_str[chromosome_i][gene_i],
                      gene_info.get_aa_positions().size());
        return false;
      }

      for (int aa_i = 0; aa_i < max_aa_pos; ++aa_i) {
        auto aa_pos_info = gene_info.get_aa_positions()[aa_i];
        auto element = pf_genotype_str[chromosome_i][gene_i][aa_i];
        if (std::find(aa_pos_info.get_amino_acids().begin(), aa_pos_info.get_amino_acids().end(),
                      std::string(1, element))
            == aa_pos_info.get_amino_acids().end()) {
          spdlog::error("Incorrect amino acid in aa sequence");
          return false;
        }
      }

      // check number copy valid or not
      if (gene_info.get_max_copies() > 1) {
        auto copy_number = NumberHelpers::char_to_single_digit_number(
            pf_genotype_str[chromosome_i][gene_i].back());
        if (copy_number > gene_info.get_max_copies()) {
          spdlog::error("Incorrect copy number");
          return false;
        }
      }
    }
  }
  // spdlog::info("Genotype {} is valid", aa_sequence);
  return true;
}

void Genotype::calculate_daily_fitness(const GenotypeParameters::PfGenotypeInfo &gene_info) {
  daily_fitness_multiple_infection = 1.0;

  spdlog::trace("Genotype: {}", aa_sequence);
  for (int chromosome_i = 0; chromosome_i < pf_genotype_str.size(); ++chromosome_i) {
    auto chromosome_info = gene_info.chromosome_infos[chromosome_i];

    for (int gene_i = 0; gene_i < pf_genotype_str[chromosome_i].size(); ++gene_i) {
      auto res_gene_info = chromosome_info.get_genes()[gene_i];
      auto max_aa_pos = res_gene_info.get_max_copies() > 1
                            ? pf_genotype_str[chromosome_i][gene_i].size() - 1
                            : pf_genotype_str[chromosome_i][gene_i].size();

      for (int aa_i = 0; aa_i < max_aa_pos; ++aa_i) {
        // calculate cost of resistance
        auto aa_pos_info = res_gene_info.get_aa_positions()[aa_i];
        auto element = pf_genotype_str[chromosome_i][gene_i][aa_i];

        auto it = std::find(aa_pos_info.get_amino_acids().begin(),
                            aa_pos_info.get_amino_acids().end(), std::string(1, element));
        auto element_id = it - aa_pos_info.get_amino_acids().begin();

        auto cr = aa_pos_info.get_daily_crs()[element_id];

        if (res_gene_info.get_average_daily_crs() > 0) {
          daily_fitness_multiple_infection *= (1 - res_gene_info.get_average_daily_crs() * cr);
          spdlog::trace(
              "\tUsing average CRS chromosome_i: {} gene_i: {} aa_i: {} cr: {} average_daily_crs: "
              "{} cr: {} (1 - res_gene_info.average_daily_crs*cr): {}",
              chromosome_i + 1, gene_i, aa_i, cr, res_gene_info.get_average_daily_crs(), cr,
              (1 - (res_gene_info.get_average_daily_crs() * cr)));
        } else {
          daily_fitness_multiple_infection *= (1 - cr);
        }
        spdlog::trace(
            "Genotype: {} chromosome_i: {} gene_i: {} aa_i: {} cr: {} "
            "daily_fitness_multiple_infection: {}",
            aa_sequence, chromosome_i + 1, gene_i, aa_i, cr, daily_fitness_multiple_infection);
      }

      // calculate for number copy variation
      if (res_gene_info.get_max_copies() > 1) {
        auto copy_number = static_cast<int>(pf_genotype_str[chromosome_i][gene_i].back()) - 48;
        if (copy_number > 1) {
          daily_fitness_multiple_infection *=
              1 - res_gene_info.get_cnv_daily_crs()[copy_number - 1];
          spdlog::trace(
              "Genotype: {} chromosome_i: {} gene_i: {} copy_number: {} "
              "daily_fitness_multiple_infection: {}",
              aa_sequence, chromosome_i + 1, gene_i, copy_number, daily_fitness_multiple_infection);
        }
      }
    }
  }
  //  std::cout << "\n";
}

void Genotype::calculate_EC50_power_n(const GenotypeParameters::PfGenotypeInfo &gene_info,
                                      DrugDatabase* drug_db) {
  EC50_power_n.resize(drug_db->size());
  for (const auto &dt : *drug_db) { EC50_power_n[dt->id()] = dt->base_EC50; }

  for (int chromosome_i = 0; chromosome_i < pf_genotype_str.size(); ++chromosome_i) {
    auto chromosome_info = gene_info.chromosome_infos[chromosome_i];

    for (int gene_i = 0; gene_i < pf_genotype_str[chromosome_i].size(); ++gene_i) {
      auto res_gene_info = chromosome_info.get_genes()[gene_i];
      auto max_aa_pos = res_gene_info.get_max_copies() > 1
                            ? pf_genotype_str[chromosome_i][gene_i].size() - 1
                            : pf_genotype_str[chromosome_i][gene_i].size();
      std::vector<int> number_of_effective_mutations_in_same_genes(drug_db->size(), 0);

      for (int aa_i = 0; aa_i < max_aa_pos; ++aa_i) {
        // calculate cost of resistance
        auto aa_pos_info = res_gene_info.get_aa_positions()[aa_i];
        auto element = pf_genotype_str[chromosome_i][gene_i][aa_i];
        auto it = std::find(aa_pos_info.get_amino_acids().begin(),
                            aa_pos_info.get_amino_acids().end(), std::string(1, element));
        if (it == aa_pos_info.get_amino_acids().end()) {
          spdlog::error("Incorrect AA in aa sequence");
        }
        auto element_id = it - aa_pos_info.get_amino_acids().begin();

        for (const auto &dt : *drug_db) {
          for (auto const &ec50s : aa_pos_info.get_multiplicative_effect_on_EC50()) {
            if (ec50s.get_drug_id() == dt->id()) {
              auto multiplicative_effect_factor = ec50s.get_factors()[element_id];
              if (multiplicative_effect_factor > 1) {
                // encounter resistant aa
                number_of_effective_mutations_in_same_genes[dt->id()] += 1;
                if (number_of_effective_mutations_in_same_genes[dt->id()] > 1) {
                  // if multiplicative effect can apply to this drug
                  for (auto const &ec50s_2_or_more :
                       res_gene_info.get_multiplicative_effect_on_ec50_for_2_or_more_mutations()) {
                    if (ec50s_2_or_more.get_drug_id() == dt->id()) {
                      multiplicative_effect_factor = ec50s_2_or_more.get_factor();
                      spdlog::trace(
                          "aa_sequence: {} DOUBLE MUT drug_id: {} chr: {} gene: {} aa: {} "
                          "EC50_power_n: {} * multiplicative_effect_factor: {}  = {}",
                          aa_sequence, dt->id(), chromosome_i + 1, gene_i, aa_i,
                          EC50_power_n[dt->id()], multiplicative_effect_factor,
                          EC50_power_n[dt->id()] * multiplicative_effect_factor);
                    }
                  }
                }
                spdlog::trace(
                    "aa_sequence: {} SINGLE MUT drug_id: {} chr: {} gene: {} aa: {} "
                    "EC50_power_n: {} * multiplicative_effect_factor: {}  = {}",
                    aa_sequence, dt->id(), chromosome_i + 1, gene_i, aa_i, EC50_power_n[dt->id()],
                    multiplicative_effect_factor,
                    EC50_power_n[dt->id()] * multiplicative_effect_factor);
              }
              EC50_power_n[dt->id()] *= multiplicative_effect_factor;
            }
          }
        }
      }

      // calculate for number copy variation
      if (res_gene_info.get_max_copies() > 1) {
        auto copy_number = static_cast<int>(pf_genotype_str[chromosome_i][gene_i].back()) - 48;
        if (copy_number > 1) {
          for (const auto &dt : *drug_db) {
            for (auto const &ec50s : res_gene_info.get_cnv_multiplicative_effect_on_EC50()) {
              if (ec50s.get_drug_id() == dt->id()) {
                spdlog::trace(
                    "aa_sequence: {} CNV drug_id: {} chr: {} gene: {} EC50_power_n: {} * "
                    "multiplicative_effect_factor: {}  = {}",
                    aa_sequence, dt->id(), chromosome_i + 1, gene_i, EC50_power_n[dt->id()],
                    ec50s.get_factors()[copy_number - 1],
                    EC50_power_n[dt->id()] * ec50s.get_factors()[copy_number - 1]);
                EC50_power_n[dt->id()] *= ec50s.get_factors()[copy_number - 1];
              }
            }
          }
        }
      }
    }
  }

  // power n
  for (const auto &dt : *drug_db) {
    EC50_power_n[dt->id()] = pow(EC50_power_n[dt->id()], dt->n());
  }
}

Genotype* Genotype::perform_mutation_by_drug(Config* p_config, utils::Random* p_random,
                                             DrugType* p_drug_type,
                                             double mutation_probability_by_locus) const {
  std::string new_aa_sequence{aa_sequence};
  for (auto aa_pos : p_drug_type->resistant_aa_locations) {
    // get aa position info (aa index in aa string, is copy number)
    if (p_config->get_genotype_parameters().get_mutation_mask()[aa_pos.aa_index_in_aa_string]
        == '1') {
      const auto p_mutation = p_random->random_flat(0.0, 1.0);
      // std::cout << "p: " << p << " Mutation probability: " << mutation_probability_by_locus <<
      // std::endl; std::cout << "aa_pos_id: " << aa_pos_id << " aa_pos.aa_index_in_aa_string: "
      // << aa_pos.aa_index_in_aa_string << " aa_pos.is_copy_number: " << aa_pos.is_copy_number <<
      // std::endl;
      if (p_mutation < mutation_probability_by_locus) {
        if (aa_pos.is_copy_number) {
          // Drug-driven CNV mutation can step one copy up or down, but must remain within the
          // configured [1, max_copies] range.
          auto old_copy_number =
              NumberHelpers::char_to_single_digit_number(aa_sequence[aa_pos.aa_index_in_aa_string]);
          const auto max_copy_number = p_config->get_genotype_parameters()
                                           .get_pf_genotype_info()
                                           .chromosome_infos[aa_pos.chromosome_id]
                                           .get_genes()[aa_pos.gene_id]
                                           .get_max_copies();
          const auto proposed_copy_number =
              p_random->random_uniform() < 0.5 ? old_copy_number - 1 : old_copy_number + 1;
          const auto new_copy_number = std::max(1, std::min(max_copy_number, proposed_copy_number));
          new_aa_sequence[aa_pos.aa_index_in_aa_string] =
              NumberHelpers::single_digit_number_to_char(new_copy_number);
        } else {
          const std::vector<std::string> aa_list = p_config->get_genotype_parameters()
                                                       .get_pf_genotype_info()
                                                       .chromosome_infos[aa_pos.chromosome_id]
                                                       .get_genes()[aa_pos.gene_id]
                                                       .get_aa_positions()[aa_pos.aa_id]
                                                       .get_amino_acids();
          // std::cout << "aa_list: [" << aa_list[0] << "," << aa_list[1] << "]" << std::endl;
          // draw random aa id
          auto new_aa_id = p_random->random_uniform(aa_list.size() - 1);
          // std::cout << "pRandom->random_uniform(aa_list.size() - 1) " <<
          // pRandom->random_uniform(aa_list.size() - 1) << std::endl;
          auto old_aa = aa_sequence[aa_pos.aa_index_in_aa_string];
          auto new_aa = aa_list[new_aa_id][0];
          // std::cout << "Mutation old_aa: " << old_aa << " -> new_aa: " << new_aa << std::endl;
          //                if (new_aa == old_aa) {
          //                    std::cout << "old_aa: " << old_aa << " == new_aa: " << new_aa;
          //                    new_aa = aa_list[new_aa_id + 1];
          //                }
          if (new_aa == old_aa) {
            if (new_aa_id + 1 < aa_list.size()) {
              new_aa = aa_list[new_aa_id + 1][0];
            } else {
              new_aa = aa_list[0][0];
            }
          }
          new_aa_sequence[aa_pos.aa_index_in_aa_string] = new_aa;
          // if (old_aa == 'C' && new_aa == 'Y'){
          //   spdlog::info("{} p: {} < {} select new_aa_id: {} from [0,{}] aa_list[new_aa_id] =
          //   aa_list[{}] = {}",
          //             Model::get_scheduler()->current_time(), p, mutation_probability_by_locus,
          //             new_aa_id, aa_list.size() - 1, new_aa_id, aa_list[new_aa_id]);
          //   spdlog::info("{} Mutation {} -> {} old: {} new: {} aa_pos_id: {} aa_pos: {}",
          //             Model::get_scheduler()->current_time(), p, mutation_probability_by_locus,
          //             old_aa, new_aa, aa_pos_id, aa_pos.aa_index_in_aa_string);
          // }
        }
      }
    }
  }
  // get genotype pointer from gene database based on aa sequence
  return Model::get_genotype_db()->get_genotype(new_aa_sequence);
}

Genotype* Genotype::perform_cnv_reversion(Config* p_config, utils::Random* p_random,
                                          DrugsInBlood* drugs_in_blood) const {
  auto new_pf_genotype_str = pf_genotype_str;
  auto mutated = false;

  const auto &pf_genotype_info = p_config->get_genotype_parameters().get_pf_genotype_info();
  const auto mutation_probability_by_locus =
      p_config->get_genotype_parameters().get_mutation_probability_per_locus();
  const auto default_cnv_reversion_multiplier =
      p_config->get_genotype_parameters().get_default_cnv_reversion_multiplier();

  for (const auto &cnv_gene_index : pf_genotype_info.get_cnv_gene_indices()) {
    const auto chromosome_i = cnv_gene_index.chromosome_index;
    const auto gene_i = cnv_gene_index.gene_index;
    const auto &gene_info = pf_genotype_info.chromosome_infos[chromosome_i].get_genes()[gene_i];
    // A gene-level multiplier overrides the global default; negative means "not configured".
    const auto reversion_multiplier = gene_info.get_cnv_reversion_multiplier() >= 0
                                          ? gene_info.get_cnv_reversion_multiplier()
                                          : default_cnv_reversion_multiplier;
    if (reversion_multiplier < 0) { continue; }
    if (new_pf_genotype_str[chromosome_i].size() <= gene_i
        || new_pf_genotype_str[chromosome_i][gene_i].empty()) {
      continue;
    }

    const auto old_copy_number =
        NumberHelpers::char_to_single_digit_number(new_pf_genotype_str[chromosome_i][gene_i].back());
    if (old_copy_number <= 1) { continue; }
    // Keep the extra copy while any active drug exposure still makes copy 2 more favorable.
    if (drug_selects_for_double_copy(gene_info, drugs_in_blood)) { continue; }

    // Iterate only CNV-capable genes here; the genotype structure is static across the run, so
    // the configuration precomputes this target list once instead of rescanning all genes.
    const auto reversion_probability =
        mutation_probability_by_locus * reversion_multiplier;
    if (reversion_probability <= 0) { continue; }

    if (p_random->random_flat(0.0, 1.0) < reversion_probability) {
      // The trailing character in the gene token encodes CNV copy number, and reversion removes
      // one copy at a time without dropping below the single-copy baseline.
      new_pf_genotype_str[chromosome_i][gene_i].back() =
          NumberHelpers::single_digit_number_to_char(std::max(1, old_copy_number - 1));
      mutated = true;
    }
  }

  if (!mutated) { return const_cast<Genotype*>(this); }
  // Reuse the genotype database so reverted sequences continue to share canonical instances.
  return Model::get_genotype_db()->get_genotype(
      Genotype::convert_pf_genotype_str_to_string(new_pf_genotype_str));
}

void Genotype::override_EC50_power_n(
    const std::vector<GenotypeParameters::OverrideEC50Pattern> &override_patterns,
    DrugDatabase* drug_db) {
  if (EC50_power_n.size() != drug_db->size()) { EC50_power_n.resize(drug_db->size()); }

  for (const auto &pattern : override_patterns) {
    if (match_pattern(pattern.get_pattern())) {
      // override ec50 power n
      EC50_power_n[pattern.get_drug_id()] =
          pow(pattern.get_ec50(), drug_db->at(pattern.get_drug_id())->n());
      spdlog::trace(
          "aa_sequence: {} OVERRIDE drug_id: {} genotype: {} EC50:{} n: {} EC50_power_n: {}",
          aa_sequence, pattern.get_drug_id(), aa_sequence, pattern.get_ec50(),
          drug_db->at(pattern.get_drug_id())->n(), EC50_power_n[pattern.get_drug_id()]);
    }
  }
}

bool Genotype::match_pattern(const std::string &pattern) {
  int id = 0;
  while (id < aa_sequence.length() && (aa_sequence[id] == pattern[id] || pattern[id] == '.')) {
    id++;
  }
  return id >= aa_sequence.length();
}

Genotype* Genotype::free_recombine_with(Config* p_config, utils::Random* p_random,
                                        Genotype* other) {
  // TODO: this function is not optimized 100%, use with care
  PfGenotypeStr new_pf_genotype_str = std::vector<ChromosomalGenotypeStr>(14);
  // for each chromosome
  for (int chromosome_id = 0; chromosome_id < pf_genotype_str.size(); ++chromosome_id) {
    if (pf_genotype_str[chromosome_id].empty()) continue;
    if (pf_genotype_str[chromosome_id].size() == 1) {
      // if single gene
      // draw random
      auto top_or_bottom = p_random->random_uniform();
      // if < 0.5 take from current, otherwise take from other
      auto gene_str = top_or_bottom < 0.5 ? pf_genotype_str[chromosome_id][0]
                                          : other->pf_genotype_str[chromosome_id][0];
      new_pf_genotype_str[chromosome_id].push_back(gene_str);
    } else {
      // if multiple genes
      // draw random to determine whether
      // within chromosome recombination happens
      auto with_chromosome_recombination = p_random->random_uniform();
      if (with_chromosome_recombination < p_config->get_parasite_parameters()
                                              .get_recombination_parameters()
                                              .get_within_chromosome_recombination_rate()) {
        // if happen draw a random crossover point based on ','
        auto cutting_gene_id =
            p_random->random_uniform(pf_genotype_str[chromosome_id].size() - 1) + 1;
        // draw another random to do top-bottom or bottom-top cross over
        auto top_or_bottom = p_random->random_uniform();
        for (auto gene_id = 0; gene_id < cutting_gene_id; ++gene_id) {
          auto gene_str = top_or_bottom < 0.5 ? pf_genotype_str[chromosome_id][gene_id]
                                              : other->pf_genotype_str[chromosome_id][gene_id];
          new_pf_genotype_str[chromosome_id].push_back(gene_str);
        }
        for (auto gene_id = cutting_gene_id; gene_id < pf_genotype_str[chromosome_id].size();
             ++gene_id) {
          auto gene_str = top_or_bottom < 0.5 ? other->pf_genotype_str[chromosome_id][gene_id]
                                              : pf_genotype_str[chromosome_id][gene_id];
          new_pf_genotype_str[chromosome_id].push_back(gene_str);
        }
      } else {
        // if not do the same with single gene
        auto top_or_bottom = p_random->random_uniform();
        for (int gene_id = 0; gene_id < pf_genotype_str[chromosome_id].size(); ++gene_id) {
          auto gene_str = top_or_bottom < 0.5 ? pf_genotype_str[chromosome_id][gene_id]
                                              : other->pf_genotype_str[chromosome_id][gene_id];

          new_pf_genotype_str[chromosome_id].push_back(gene_str);
        }
      }
    }
  }

  auto new_aa_sequence = convert_pf_genotype_str_to_string(new_pf_genotype_str);
  return Model::get_genotype_db()->get_genotype(new_aa_sequence);
}

std::string Genotype::convert_pf_genotype_str_to_string(const PfGenotypeStr &pf_genotype_str) {
  std::stringstream ss;

  for (const auto &chromosome : pf_genotype_str) {
    for (const auto &gene : chromosome) {
      ss << gene;
      if (gene != chromosome.back()) { ss << ","; }
    }

    if (&chromosome != &pf_genotype_str.back()) { ss << "|"; }
  }

  return ss.str();
}
Genotype* Genotype::free_recombine(Config* config, utils::Random* p_random, Genotype* female,
                                   Genotype* male) {
  PfGenotypeStr new_pf_genotype_str = std::vector<ChromosomalGenotypeStr>(14);
  // for each chromosome
  for (int chromosome_id = 0; chromosome_id < female->pf_genotype_str.size(); ++chromosome_id) {
    if (female->pf_genotype_str[chromosome_id].empty()) continue;
    if (female->pf_genotype_str[chromosome_id].size() == 1) {
      // if single gene
      // draw random
      auto top_or_bottom = p_random->random_uniform();
      // if < 0.5 take from current, otherwise take from other
      auto gene_str = top_or_bottom < 0.5 ? female->pf_genotype_str[chromosome_id][0]
                                          : male->pf_genotype_str[chromosome_id][0];
      new_pf_genotype_str[chromosome_id].push_back(gene_str);
    } else {
      // if multiple genes
      // draw random to determine whether
      // within chromosome recombination happens
      auto with_chromosome_recombination = p_random->random_uniform();
      if (with_chromosome_recombination < config->get_parasite_parameters()
                                              .get_recombination_parameters()
                                              .get_within_chromosome_recombination_rate()) {
        // if happen draw a random crossover point based on ','
        auto cutting_gene_id =
            p_random->random_uniform(female->pf_genotype_str[chromosome_id].size() - 1) + 1;
        // draw another random to do top-bottom or bottom-top cross over
        auto top_or_bottom = p_random->random_uniform();
        for (auto gene_id = 0; gene_id < cutting_gene_id; ++gene_id) {
          auto gene_str = top_or_bottom < 0.5 ? female->pf_genotype_str[chromosome_id][gene_id]
                                              : male->pf_genotype_str[chromosome_id][gene_id];
          new_pf_genotype_str[chromosome_id].push_back(gene_str);
        }
        for (auto gene_id = cutting_gene_id;
             gene_id < female->pf_genotype_str[chromosome_id].size(); ++gene_id) {
          auto gene_str = top_or_bottom < 0.5 ? male->pf_genotype_str[chromosome_id][gene_id]
                                              : female->pf_genotype_str[chromosome_id][gene_id];
          new_pf_genotype_str[chromosome_id].push_back(gene_str);
        }
      } else {
        // if there is no within chromosome recombination
        // do the same with single gene
        auto top_or_bottom = p_random->random_uniform();
        for (int gene_id = 0; gene_id < female->pf_genotype_str[chromosome_id].size(); ++gene_id) {
          auto gene_str = top_or_bottom < 0.5 ? female->pf_genotype_str[chromosome_id][gene_id]
                                              : male->pf_genotype_str[chromosome_id][gene_id];

          new_pf_genotype_str[chromosome_id].push_back(gene_str);
        }
      }
    }
  }

  auto new_aa_sequence = convert_pf_genotype_str_to_string(new_pf_genotype_str);
  return Model::get_genotype_db()->get_genotype(new_aa_sequence);
}
