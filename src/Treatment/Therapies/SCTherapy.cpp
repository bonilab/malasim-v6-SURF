#include "SCTherapy.h"

#include "Configuration/Config.h"
#include "Simulation/Model.h"

void SCTherapy::add_drug(int drug_id) {
  Therapy::add_drug(drug_id);
  if (drug_id == 0) {
    artemisinin_id = drug_id;
  }
}

void SCTherapy::calculate_max_dosing_day() {
  max_dosing_day = dosing_day[0];
  for (auto ndx = 1; ndx < dosing_day.size(); ndx++) {
    max_dosing_day =
        (dosing_day[ndx] > max_dosing_day) ? dosing_day[ndx] : max_dosing_day;
  }
}

int SCTherapy::get_max_dosing_day() const {
  // The calculate_max_dosing_day must be called before this function
  if(max_dosing_day == -1) {
    throw std::runtime_error("Max dosing day has not been calculated. Call calculate_max_dosing_day first.");
  }

  return max_dosing_day;
}

void SCTherapy::print(std::ostream &os) const {
  os << Model::get_drug_db()->at(drug_ids[0])->name();
  for (std::size_t i = 1; i < drug_ids.size(); ++i) {
    os << "+" + Model::get_drug_db()->at(drug_ids[i])->name();
  }
  os << "(" << dosing_day[0];
  for (auto ndx = 1; ndx < dosing_day.size(); ndx++) {
    os << "," << dosing_day[ndx];
  }
  os << ")";
}
