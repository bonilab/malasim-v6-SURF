#include "DrugType.h"

#include <cmath>

#include "Parasites/Genotype.h"
#include "Simulation/Model.h"

#ifndef LOG2_10
#define LOG2_10 3.32192809489
#endif

double DrugType::get_parasite_killing_rate_by_concentration(const double &concentration,
                                                            const double &ec50_power_n) {
  const auto con_power_n = pow(concentration, n_);
  // std::cout << "c: " << concentration << " n: " << n_ << " p_max: " <<
  // maximum_parasite_killing_rate_ << " EC50_power_n: " << EC50_power_n << "
  // con_power_n: " << con_power_n << std::endl;
  const auto killing_perday =
      maximum_parasite_killing_rate_ * (con_power_n / (con_power_n + ec50_power_n));
  // std::cout<< "c: " << concentration << " n: " << n_ << " ppr: "<<
  // killing_perday << std::endl;
  return killing_perday;
}

double DrugType::n() { return n_; }

void DrugType::set_n(const double &n) {
  n_ = n;
  //    set_EC50_power_n(pow(EC50_, n_));
}

int DrugType::get_total_duration_of_drug_activity(const int &dosing_days) const {
  // CutOffPercent is 10
  // log2(100.0 / 10.0) = 3.32192809489
  return static_cast<int>(dosing_days + ceil(drug_half_life_ * LOG2_10));
}

void DrugType::populate_resistant_aa_locations() {
  resistant_aa_locations_.clear();
  for (const auto &chromosome_info :
       Model::get_config()->get_genotype_parameters().get_pf_genotype_info().chromosome_infos) {
    for (int gene_id = 0; gene_id < chromosome_info.get_genes().size(); ++gene_id) {
      const auto &gene_info = chromosome_info.get_genes()[gene_id];

      auto aa_pos_in_sequence =
          Model::get_config()->get_genotype_parameters().get_pf_genotype_info().calculate_aa_pos(
              chromosome_info.get_chromosome_id() - 1, gene_id, 0);

      for (int aa_id = 0; aa_id < gene_info.get_aa_positions().size(); ++aa_id) {
        for (auto const &multiple_ec50 :
             gene_info.get_aa_positions()[aa_id].get_multiplicative_effect_on_EC50()) {
          if (multiple_ec50.get_drug_id() == id_) {
            resistant_aa_locations_.push_back(
                {.chromosome_id = chromosome_info.get_chromosome_id() - 1,
                 .gene_id = gene_id,
                 .aa_id = aa_id,
                 .aa_index_in_aa_string = aa_pos_in_sequence,
                 .is_copy_number = false});
          }
        }
        aa_pos_in_sequence++;
      }
      if (gene_info.get_max_copies() > 1) {
        for (auto const &cnv_multiple_ec50 : gene_info.get_cnv_multiplicative_effect_on_EC50()) {
          if (cnv_multiple_ec50.get_drug_id() == id_) {
            resistant_aa_locations_.push_back(
                {.chromosome_id = chromosome_info.get_chromosome_id() - 1,
                 .gene_id = gene_id,
                 .aa_id = static_cast<int>(gene_info.get_aa_positions().size()),
                 .aa_index_in_aa_string = aa_pos_in_sequence,
                 .is_copy_number = true});
          }
        }
      }
    }
  }
}

