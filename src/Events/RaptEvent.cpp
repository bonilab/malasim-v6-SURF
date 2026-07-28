/*
 * RaptEvent.cpp
 *
 * Implement the RAPT event class. When this event is active we assume that it
 * will run at least once per year, similar to individual birthday events. It
 * may run more often depending upon the period (in months) that the event has
 * been configured for.
 */
#include "RaptEvent.h"

#include <date/date.h>

#include "Configuration/Config.h"
#include "Core/Scheduler/Scheduler.h"
#include "Events/TestTreatmentFailureEvent.h"
#include "Population/Person/Person.h"
#include "Simulation/Model.h"

void RaptEvent::do_execute() {
  auto* person = get_person();
  if (person == nullptr) { throw std::runtime_error("Person is nullptr"); }

  const auto &rapt_config = Model::get_config()->get_rapt_settings();

  // Check to see if we should receive a therapy: RAPT is currently active, the
  // person is the correct age, and they have not recently taken a treatment in
  // the past 28 days (based on testing for treatment failure).
  if (Model::get_scheduler()->current_time() >= rapt_config.get_start_date_as_day()
      && person->get_age() >= rapt_config.get_age_start()
      && !person->has_event<TestTreatmentFailureEvent>()) {
    // Is their base compliance over-5 or under-5 treatment rate, presume that
    // we are following the Malaria Indicator Survey convention of under-5 being
    // 0 - 59 months.
    auto pr_treatment =
        person->get_age() < 5
            ? Model::get_config()->location_db()[person->get_location()].p_treatment_under_5
            : Model::get_config()->location_db()[person->get_location()].p_treatment_over_5;

    // Adjust the probability based upon the configured compliance rate with
    // RAPT
    double pr_rapt = pr_treatment * rapt_config.get_compliance();
    const auto pr = Model::get_random()->random_flat(0.0, 1.0);

    if (pr <= pr_rapt) {
      person->receive_therapy(Model::get_therapy_db()[rapt_config.get_therapy_id()].get(), nullptr);
    }
  }

  // Determine when the next RAPT dose should take place based upon scheduling
  // period
  const auto ymd = Model::get_scheduler()->get_ymd_after_months(
      Model::get_config()->get_rapt_settings().get_period());

  // Find the first and last day of the month of the next dose
  const auto first_day = date::year_month_day(ymd.year(), ymd.month(), date::day(1));
  const auto last_day = date::year_month_day_last{ymd.year(), ymd.month() / date::last};

  // Following adjustment scheduler days, conduct a uniform draw across the next
  // dosing month to determine the actual next date when the RAPT dose may be
  // taken
  const auto from = Model::get_scheduler()->get_days_to_ymd(first_day);
  const auto to = Model::get_scheduler()->get_days_to_ymd(last_day);
  const auto days_to_next_event = Model::get_random()->random_uniform<int>(from, to + 1);

  // Schedule the event
  person->schedule_rapt_event(days_to_next_event);
}
