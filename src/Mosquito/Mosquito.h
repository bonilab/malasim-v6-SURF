#ifndef POMS_SRC_MOSQUITO_MOSQUITO_H
#define POMS_SRC_MOSQUITO_MOSQUITO_H
#include <vector>

#include "Configuration/Config.h"

class Genotype;
class Model;
class Config;
class Population;

using MosquitoRecombinedGenotypeInfo =
    std::pair<std::vector<std::pair<int, std::string>>, std::pair<int, int>>;
class Mosquito {
public:
  Mosquito(const Mosquito &) = delete;
  Mosquito &operator=(const Mosquito &) = delete;
  Mosquito(Mosquito &&) = delete;
  Mosquito &operator=(Mosquito &&) = delete;

  explicit Mosquito() = default;
  Mosquito(std::map<int, double> drug_id_min_ec50,
           std::vector<std::vector<std::vector<Genotype*>>> genotypes_table)
      : drug_id_min_ec50(std::move(drug_id_min_ec50)),
        genotypes_table(std::move(genotypes_table)) {}
  virtual ~Mosquito() = default;

  void initialize(Config* config);

  void infect_new_cohort_in_prmc(Config* config,
                                 utils::Random* random,
                                 Population* population,
                                 int tracking_index);

  void infect_new_cohort_at_location(Config* config,
                                     utils::Random* random,
                                     Population* population,
                                     int location,
                                     int tracking_index);

  std::vector<std::vector<std::vector<Genotype*>>> genotypes_table; /* Mosquito table */

  [[nodiscard]] static std::vector<unsigned int> build_interrupted_feeding_indices(
      utils::Random* random, const double &interrupted_feeding_rate, const int &prmc_size);

  int random_genotype(int location, int tracking_index);

  // this function will populate values for both parasite densities and genotypes that carried by a
  // person
  static void get_genotypes_profile_from_person(Person* person,
                                                std::vector<Genotype*> &sampling_genotypes,
                                                std::vector<double> &relative_infectivity_each_pp);

  static std::string get_old_genotype_string(std::string new_genotype);
  static std::string get_old_genotype_string2(std::string new_genotype);

  Model* model{nullptr};
  std::map<int, double> drug_id_min_ec50;
  using ResistantDrugInfo = std::pair<std::vector<std::string>, std::vector<int>>;
  std::vector<ResistantDrugInfo> resistant_drug_list = {
      ResistantDrugInfo{{"DHA-PPQ:2-2"}, {0, 3}},
      ResistantDrugInfo{{"ASAQ:2-2", "ASAQ:2-3", "ASAQ:2-4", "ASAQ:2"}, {0, 2}},
      ResistantDrugInfo{{"AL:2-2", "AL:2-3", "AL:2-4", "AL:2"}, {0, 1}},
      ResistantDrugInfo{{"DHA-PPQ-AQ:3-3", "DHA-PPQ-AQ:3-4", "DHA-PPQ-AQ:3-5", "DHA-PPQ-AQ:3"},
                        {0, 3, 2}},
      ResistantDrugInfo{{"DHA-PPQ-LUM:3-3", "DHA-PPQ-LUM:3-4", "DHA-PPQ-LUM:3-5", "DHA-PPQ-LUM:3"},
                        {0, 3, 1}},
  };
};

#endif  // POMS_SRC_MOSQUITO_MOSQUITO_H
