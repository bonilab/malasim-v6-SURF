#include "NonInfantImmuneComponent.h"
#include "Simulation/Model.h"
#include "Configuration/Config.h"


//OBJECTPOOL_IMPL(NonInfantImmuneComponent)

NonInfantImmuneComponent::NonInfantImmuneComponent(ImmuneSystem *immune_system) : ImmuneComponent(immune_system) {
}

NonInfantImmuneComponent::~NonInfantImmuneComponent() {
}

double NonInfantImmuneComponent::get_acquire_rate(const int &age) const {
  // Config is singleton - already cached at Model level
  return (age > 80) ? Model::get_config()->get_immune_system_parameters().acquire_rate_by_age[80]
                    : Model::get_config()->get_immune_system_parameters().acquire_rate_by_age[age];
}

double NonInfantImmuneComponent::get_decay_rate(const int &age) const {
  // Config is singleton - already cached at Model level
  return Model::get_config()->get_immune_system_parameters().decay_rate;
}