#include "ImmuneSystem.h"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstddef>

#include "Configuration/Config.h"
#include "Core/Scheduler/Scheduler.h"
#include "Population/ImmuneSystem/ImmuneSystemConstants.h"
#include "Population/Person/Person.h"
#include "Simulation/Model.h"
#include "Utils/Random.h"

namespace {

std::size_t age_index(core::Age age) {
  return static_cast<std::size_t>(std::min(age, immune::K_MAX_IMMUNE_AGE_INDEX));
}

}  // namespace

ImmuneSystem::ImmuneSystem(Person* person) : person_(person) {}

ImmuneSystem::~ImmuneSystem() = default;

void ImmuneSystem::draw_random_immune() {
  const auto &parameters = Model::get_config()->get_immune_system_parameters();
  latest_value_ = Model::get_random()->random_beta(parameters.alpha_immune, parameters.beta_immune);
}

double ImmuneSystem::get_latest_immune_value() const { return latest_value_; }

void ImmuneSystem::set_latest_immune_value(double value) { latest_value_ = value; }

double ImmuneSystem::get_current_value() const {
  if (person_ == nullptr) { return 0.0; }

  const auto duration = Model::get_scheduler()->current_time() - person_->get_latest_update_time();
  if (duration == 0) { return latest_value_; }

  if (mode_ == ImmuneSystemMode::INFANT) {
    const auto factor = duration == 1 ? immune::K_ONE_DAY_INFANT_DECAY_FACTOR
                                      : std::exp(-immune::K_INFANT_IMMUNE_DECAY_RATE * duration);
    return latest_value_ * factor;
  }

  const auto age = person_->get_age();
  if (increase_) {
    const auto factor = duration == 1 ? get_one_day_acquire_factor(age)
                                      : std::exp(-get_acquire_rate(age) * duration);
    return 1 - ((1 - latest_value_) * factor);
  }

  const auto factor =
      duration == 1 ? get_one_day_decay_factor() : std::exp(-get_decay_rate(age) * duration);
  const auto value = latest_value_ * factor;
  return value < immune::K_IMMUNE_VALUE_CUTOFF ? 0.0 : value;
}

void ImmuneSystem::switch_to_non_infant() {
  if (mode_ == ImmuneSystemMode::NON_INFANT) { return; }
  assert(person_ == nullptr
         || person_->get_latest_update_time() == Model::get_scheduler()->current_time());
  mode_ = ImmuneSystemMode::NON_INFANT;
}

double ImmuneSystem::get_parasite_size_after_t_days(const int &duration,
                                                    const double &original_size,
                                                    const double &fitness) const {
  const auto last_immune_level = get_latest_immune_value();
  const auto temp =
      (Model::get_config()->get_immune_system_parameters().c_max * (1 - last_immune_level))
      + (Model::get_config()->get_immune_system_parameters().c_min * last_immune_level);
  // std::cout << "day: " << Model::get_scheduler()->current_time() << "\tc_max: " <<
  // Model::CONFIG->immune_system_information().c_max << "\tc_min: " <<
  // Model::CONFIG->immune_system_information().c_min << "\tlast_immune_level: " <<
  // last_immune_level << "\ttemp: " << temp << std::endl;
  //  std::cout << "Day: " << Model::get_scheduler()->current_time() << "\tImmune: old density: " <<
  //  originalSize << "\t duration: " << duration << "\tfitness: "
  //  << fitness << "\tlast immune level: " << last_immune_level << "\ttemp: " << temp;
  const auto value = original_size + (duration * (log10(temp) + log10(fitness)));
  //  std::cout << "\tnew density: " << value << std::endl;
  return value;
}

double ImmuneSystem::get_clinical_progression_probability() const {
  const auto immune = get_current_value();

  const auto isf = Model::get_config()->get_immune_system_parameters();

  //    double PClinical = (isf.min_clinical_probability - isf.max_clinical_probability) *
  //    pow(immune, isf.immune_effect_on_progression_to_clinical) + isf.max_clinical_probability;

  //    const double p_m = 0.99;

  const auto p_clinical =
      isf.max_clinical_probability
      / (1 + pow((immune / isf.midpoint), isf.immune_effect_on_progression_to_clinical));

  // spdlog::info(
  //     "ImmuneSystem::get_clinical_progression_probability: immune: {}, PClinical: {}, max
  //     clinical " "probability: {}, immune effect on progression to clinical: {}", immune,
  //     p_clinical, isf.max_clinical_probability, isf.immune_effect_on_progression_to_clinical);
  //    std::cout << immune << "\t" << PClinical<< std::endl;
  return p_clinical;
}

void ImmuneSystem::update() { latest_value_ = get_current_value(); }

double ImmuneSystem::get_one_day_decay_factor() const {
  if (mode_ == ImmuneSystemMode::INFANT) { return immune::K_ONE_DAY_INFANT_DECAY_FACTOR; }
  return Model::get_config()->get_immune_system_parameters().decay_rate_one_day_factor;
}

double ImmuneSystem::get_one_day_acquire_factor(core::Age age) const {
  if (mode_ == ImmuneSystemMode::INFANT) { return 1.0; }
  const auto &parameters = Model::get_config()->get_immune_system_parameters();
  return parameters.acquire_rate_by_age_one_day_factor[age_index(age)];
}

double ImmuneSystem::get_decay_rate(core::Age) const {
  if (mode_ == ImmuneSystemMode::INFANT) { return immune::K_INFANT_IMMUNE_DECAY_RATE; }
  return Model::get_config()->get_immune_system_parameters().decay_rate;
}

double ImmuneSystem::get_acquire_rate(core::Age age) const {
  if (mode_ == ImmuneSystemMode::INFANT) { return 0.0; }
  const auto &parameters = Model::get_config()->get_immune_system_parameters();
  return parameters.acquire_rate_by_age[age_index(age)];
}
