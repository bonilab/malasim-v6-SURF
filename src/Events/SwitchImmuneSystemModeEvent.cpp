#include "SwitchImmuneSystemModeEvent.h"

#include <stdexcept>

#include "Population/ImmuneSystem/ImmuneSystem.h"
#include "Population/Person/Person.h"
#include "spdlog/spdlog.h"

SwitchImmuneSystemModeEvent::SwitchImmuneSystemModeEvent(Person* person) : PersonEvent(person) {
  if (person == nullptr) {
    throw std::invalid_argument("SwitchImmuneSystemModeEvent: person is nullptr");
  }
}

SwitchImmuneSystemModeEvent::~SwitchImmuneSystemModeEvent() = default;

void SwitchImmuneSystemModeEvent::do_execute() {
  auto* person = get_person();
  if (person == nullptr) {
    spdlog::error("SwitchImmuneSystemModeEvent::do_execute: person is nullptr");
    throw std::invalid_argument("SwitchImmuneSystemModeEvent::do_execute: person is nullptr");
  }
  person->get_immune_system()->switch_to_non_infant();
}
