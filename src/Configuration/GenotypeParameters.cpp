#include "GenotypeParameters.h"

#include <stdexcept>

#include "Simulation/Model.h"

void GenotypeParameters::process_config_with_number_of_locations(size_t number_of_locations) {
  spdlog::info("Processing GenotypeParameters");
  for (const auto &initial_genotype_info_raw : get_initial_parasite_info_raw()) {
    const auto location = initial_genotype_info_raw.get_location_id();
    const auto location_from = location == -1 ? 0 : location;
    const auto location_to =
        location == -1 ? number_of_locations : std::min<size_t>(location + 1, number_of_locations);

    // apply for all location
    for (auto loc = location_from; loc < location_to; ++loc) {
      for (const auto &parasite_node : initial_genotype_info_raw.get_parasite_info()) {
        const auto &aa_sequence = parasite_node.get_aa_sequence();
        auto parasite_type_id = Model::get_genotype_db()->get_id(aa_sequence);
        auto prevalence = parasite_node.get_prevalence();
        initial_parasite_info_.emplace_back(loc, parasite_type_id, prevalence);
      }
    }
    // for(auto &initial_genotype_info : get_initial_parasite_info()) {
    // spdlog::debug("Location: {} parasite_type_id: {} prevalence: {}",
    //              initial_genotype_info.location, initial_genotype_info.parasite_type_id,
    // initial_genotype_info.prevalence);
    // }
  }
}


void GenotypeParameters::validate_cnv_reversion_multipliers(const GenotypeParameters &params) {
  constexpr double kMaxCnvReversionMultiplier = 10.0;
  if (params.get_default_cnv_reversion_multiplier() > kMaxCnvReversionMultiplier) {
    throw std::invalid_argument("Default CNV reversion multiplier should be in range [0,10]");
  }
  for (const auto &chromosome_info : params.get_pf_genotype_info().chromosome_infos) {
    for (const auto &gene_info : chromosome_info.get_genes()) {
      if (gene_info.get_cnv_reversion_multiplier() < 0) { continue; }
      if (gene_info.get_max_copies() <= 1) {
        throw std::invalid_argument(
            "CNV reversion multiplier can only be set for genes with max_copies > 1");
      }
      if (gene_info.get_cnv_reversion_multiplier() > kMaxCnvReversionMultiplier) {
        throw std::invalid_argument("CNV reversion multiplier should be in range [0,10]");
      }
    }
  }
}
