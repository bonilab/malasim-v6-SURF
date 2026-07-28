#include "BirthdayEvent.h"

#include <cassert>

#include "Core/Scheduler/Scheduler.h"
#include "Population/Person/Person.h"
#include "Simulation/Model.h"

// OBJECTPOOL_IMPL(BirthdayEvent)

void BirthdayEvent::do_execute() {
  // spdlog::info("Time: {}, BirthdayEvent::do_execute, person age: {}",
  //              Model::get_scheduler()->current_time(), get_person()->get_age());
  auto* person = get_person();
  if (person == nullptr) {
    spdlog::error("BirthdayEvent::do_execute, person is nullptr");
    throw std::invalid_argument("BirthdayEvent::do_execute, person is nullptr");
  }
  person->increase_age_by_1_year();

  const auto days_to_next_birthday = Model::get_scheduler()->get_days_until_next_year_anniversary();
  // spdlog::info("Time: {}, Schedule BirthdayEvent::do_execute, person age: {},
  // days_to_next_birthday: {}",
  //              Model::get_scheduler()->current_time(), person->get_age(), days_to_next_birthday);

  person->schedule_birthday_event(days_to_next_birthday);
}
