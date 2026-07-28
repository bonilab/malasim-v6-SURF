/*
 * AnnualBetaUpdateEvent.hxx
 *
 * This event class provides the ability to alter the beta parameter on an
 * annual basis.
 */
#ifndef ANNUALBETAUPDATEEVENT_HXX
#define ANNUALBETAUPDATEEVENT_HXX

#include <cmath>

#include "Configuration/Config.h"
#include "Core/Scheduler/Scheduler.h"
#include "Events/Event.h"
#include "Simulation/Model.h"

class AnnualBetaUpdateEvent : public WorldEvent {
public:
  AnnualBetaUpdateEvent(const AnnualBetaUpdateEvent &) = delete;
  AnnualBetaUpdateEvent(AnnualBetaUpdateEvent &&) = delete;
  AnnualBetaUpdateEvent &operator=(const AnnualBetaUpdateEvent &) = delete;
  AnnualBetaUpdateEvent &operator=(AnnualBetaUpdateEvent &&) = delete;
  AnnualBetaUpdateEvent(float rate, int start) : rate_(rate) { set_time(start); }
  ~AnnualBetaUpdateEvent() override = default;

  static constexpr std::string_view EVENT_NAME{"annual_beta_update"};
  [[nodiscard]] std::string_view name() const noexcept override { return EVENT_NAME; }

private:
  float rate_ = 0.0;

  // Execute the annual beta update event
  void do_execute() override {
    // Grab a reference to the location_db to work with
    auto &location_db = Model::get_config()->location_db();

    // Iterate through and adjust the betas
    auto count = Model::get_config()->number_of_locations();
    for (auto ndx = 0; ndx < count; ndx++) {
      location_db[ndx].beta = adjust(location_db[ndx].beta, rate_);
    }

    // Schedule for one year from now
    auto time = Model::get_scheduler()->get_days_until_next_year_anniversary();
    auto event = std::make_unique<AnnualBetaUpdateEvent>(rate_, time);
    Model::get_scheduler()->schedule_population_event(std::move(event));

    // Log on demand
    spdlog::debug("Annual beta update event: {} - {} {}",
                  Model::get_scheduler()->get_current_date_string(), rate_, location_db[0].beta);
  }

  // Update the bete by the given rate, round up ate the fifth decimal place but
  // clamp at 0.0
  static double adjust(double beta, double rate) {
    beta += (beta * rate);
    beta = static_cast<int>(beta * pow(10, 5)) / pow(10, 5);
    return (beta < 0.0) ? 0.0 : beta;
  }
};

#endif
