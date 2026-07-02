#ifndef Genotype_H
#define Genotype_H

#include "Configuration/GenotypeParameters.h"
#include "Utils/Random.h"

class GenotypeParameters;

class PfGenotypeInfo;

class OverrideEC50Pattern;

class DrugDatabase;

class DrugType;

class DrugsInBlood;

class Therapy;

class Config;

class Random;

using GeneStr = std::string;
using ChromosomalGenotypeStr = std::vector<GeneStr>;
using PfGenotypeStr = std::vector<ChromosomalGenotypeStr>;
using MosquitoRecombinedGenotypeInfo =
    std::pair<std::vector<std::pair<int, std::string>>, std::pair<int, int>>;

class Genotype {
public:
  Genotype(const Genotype &) = delete;
  Genotype(Genotype &&) = delete;
  Genotype &operator=(const Genotype &) = delete;
  Genotype &operator=(Genotype &&) = delete;

  explicit Genotype(const std::string &aa_sequence);
  virtual ~Genotype();

  int genotype_id_{-1};
  PfGenotypeStr pf_genotype_str = std::vector<ChromosomalGenotypeStr>(14);
  std::string aa_sequence;
  double daily_fitness_multiple_infection{1};
  std::vector<double> EC50_power_n;
  std::vector<MosquitoRecombinedGenotypeInfo> resistant_recombinations_in_mosquito;

  [[nodiscard]] int genotype_id() const { return genotype_id_; }
  void set_genotype_id(int genotype_id) { genotype_id_ = genotype_id; }

  double get_EC50_power_n(DrugType* dt);

  bool resist_to(DrugType* dt);

  Genotype* combine_mutation_to(const int &locus, const int &value);

  Genotype* modify_genotype_allele(const std::vector<std::tuple<int, int, char>> &alleles,
                                   Config* p_config) const;

  [[nodiscard]] std::string get_aa_sequence() const;

  bool is_valid(const GenotypeParameters::PfGenotypeInfo &gene_info);

  void calculate_daily_fitness(const GenotypeParameters::PfGenotypeInfo &gene_info);

  void calculate_EC50_power_n(const GenotypeParameters::PfGenotypeInfo &info,
                              DrugDatabase* p_database);

  Genotype* perform_mutation_by_drug(Config* p_config, utils::Random* p_random,
                                     DrugType* p_drug_type,
                                     double mutation_probability_by_locus) const;

  Genotype* perform_cnv_reversion(Config* p_config, utils::Random* p_random,
                                  DrugsInBlood* drugs_in_blood) const;

  friend std::ostream &operator<<(std::ostream &os, const Genotype &genotype);

  void override_EC50_power_n(
      const std::vector<GenotypeParameters::OverrideEC50Pattern> &override_patterns,
      DrugDatabase* drug_db);

  bool match_pattern(const std::string &pattern);

  Genotype* free_recombine_with(Config* config, utils::Random* p_random, Genotype* other);

  static Genotype* free_recombine(Config* config, utils::Random* p_random, Genotype* female,
                                  Genotype* mmale);

  static std::string convert_pf_genotype_str_to_string(const PfGenotypeStr &pf_genotype_str);
};

#endif /* Genotype_H */
