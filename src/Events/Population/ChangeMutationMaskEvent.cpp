#include "ChangeMutationMaskEvent.h"

#include <algorithm>

#include "Configuration/Config.h"

ChangeMutationMaskEvent::ChangeMutationMaskEvent(const std::vector<bool> &mask, const int &at_time)
    : mask_{mask} {
  set_time(at_time);
}

void ChangeMutationMaskEvent::do_execute() {
  Model::get_config()->get_genotype_parameters().set_mutation_mask(mask_);
  spdlog::info("{}: change mutation mask to {} enabled positions across {} entries",
               Model::get_scheduler()->get_current_date_string(),
               std::count(mask_.begin(), mask_.end(), true), mask_.size());
}
