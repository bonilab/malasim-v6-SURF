#include "DrugParameters.h"

#include "Simulation/Model.h"

void DrugParameters::process_config() {
  spdlog::info("Processing DrugParameters");
  for (auto drug_id = 0; drug_id < drug_infos_.size(); drug_id++) {
    auto dt = std::make_unique<DrugType>();
    dt->set_id(drug_id);

    const auto i_s = NumberHelpers::number_to_string<int>(drug_id);
    const auto &dt_node = drug_infos_.at(drug_id);

    dt->set_name(dt_node.get_name());
    dt->set_drug_half_life(dt_node.get_half_life());
    dt->set_maximum_parasite_killing_rate(dt_node.get_maximum_parasite_killing_rate());
    dt->set_n(dt_node.get_n());
    //    dt->set_EC50(node["EC50"].as<double>());

    //    std::cout <<dt->drug_half_life() << "-" << dt->maximum_parasite_killing_rate() << "-" <<
    //    dt->n() << "-" << dt->EC50() << std::endl;
    for (double value : dt_node.get_age_specific_drug_concentration_sd()) {
      dt->age_group_specific_drug_concentration_sd().push_back(value);
      dt->age_specific_drug_absorption().push_back(1.0);
    }
    //    assert(dt->age_group_specific_drug_concentration_sd().size() == 15);

    if (!dt_node.get_age_specific_drug_absorption().empty()) {
      dt->set_age_specific_drug_absorption(dt_node.get_age_specific_drug_absorption());
    }

    dt->set_k(dt_node.get_k());

    dt->set_base_ec50(dt_node.get_base_ec50());

    dt->populate_resistant_aa_locations();

    Model::get_drug_db()->add(std::move(dt));
  }
}
