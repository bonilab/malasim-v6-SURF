#ifndef SCHEDULER_H
#define SCHEDULER_H

#include <memory>

#include "Core/Scheduler/EventManager.h"
#include "Core/types.h"
#include "Events/Event.h"
#include "Utils/Helpers/StringHelpers.h"
#include "date/date.h"

class Model;

class Scheduler {
public:
  // Disable copy and assignment
  Scheduler(const Scheduler &) = delete;
  Scheduler &operator=(const Scheduler &) = delete;
  Scheduler(Scheduler &&) = delete;
  Scheduler &operator=(Scheduler &&) = delete;

  explicit Scheduler();
  virtual ~Scheduler();

  // Getter and Setter for current_time
  [[nodiscard]] core::SimDay current_time() const { return current_time_; }
  void set_current_time(core::SimDay time) { current_time_ = time; }

  // Getter and Setter for is_force_stop
  [[nodiscard]] bool is_force_stop() const { return is_force_stop_; }
  void set_is_force_stop(bool force_stop) { is_force_stop_ = force_stop; }

  void extend_total_time(int new_total_time);

  // Event management methods
  void clear_all_events() { world_events_.get_events().clear(); }

  virtual void schedule_population_event(std::unique_ptr<WorldEvent> event) {
    if (event != nullptr) { world_events_.schedule_event(std::move(event)); }
  }

  virtual void cancel(WorldEvent* event) {
    if (event != nullptr) { world_events_.cancel_all_events_except(event); }
  }

  // Core simulation methods
  void initialize(const date::year_month_day &starting_date,
                  const date::year_month_day &ending_date);

  void run();
  void begin_time_step();
  void end_time_step();
  void daily_update();

  // Time-related query methods
  [[nodiscard]] bool can_stop();
  [[nodiscard]] int get_current_day_in_year();
  [[nodiscard]] unsigned int get_current_month_in_year();
  [[nodiscard]] bool is_today_last_day_of_month();
  [[nodiscard]] bool is_today_first_day_of_month();
  [[nodiscard]] bool is_today_first_day_of_year();
  [[nodiscard]] bool is_today_last_day_of_year();
  [[nodiscard]] int get_days_to_next_year() const;
  [[nodiscard]] int get_days_to_next_n_year(int n) const;
  [[nodiscard]] int get_days_in_current_month() const;
  [[nodiscard]] int get_current_year() const;
  [[nodiscard]] unsigned int get_current_day_of_month() const;
  [[nodiscard]] date::year_month_day get_ymd_after_months(int months) const;
  [[nodiscard]] int get_days_to_ymd(const date::year_month_day &ymd) const;
  [[nodiscard]] date::year_month_day get_ymd_after_days(int days) const;
  [[nodiscard]] int get_unix_time() const;
  [[nodiscard]] date::year_month_day get_calendar_date() const;
  [[nodiscard]] std::string get_current_date_string() const {
    return StringHelpers::date_as_string(date::year_month_day{calendar_date_});
  }
  [[nodiscard]] std::string get_current_date_seperated_string() const {
    return StringHelpers::date_as_seperated_string(date::year_month_day{calendar_date_});
  }
  // Access to event manager
  EventManager<WorldEvent> &get_world_events() { return world_events_; }
  [[nodiscard]] const EventManager<WorldEvent> &get_world_events() const { return world_events_; }

private:
  core::SimDay current_time_{core::K_INVALID_SIM_DAY};
  bool is_force_stop_{false};
  date::sys_days calendar_date_;
  EventManager<WorldEvent> world_events_;  // Use EventManager for world/population events

public:
};

#endif /* SCHEDULER_H */
