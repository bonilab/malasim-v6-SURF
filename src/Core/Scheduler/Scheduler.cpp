#include "Scheduler.h"

#include <Configuration/Config.h>

#include "Core/types.h"
#include "Population/Population.h"
#include "Simulation/Model.h"
#include "Utils/Helpers/TimeHelpers.h"
#include "spdlog/spdlog.h"

Scheduler::Scheduler() = default;

Scheduler::~Scheduler() = default;

void Scheduler::initialize(const date::year_month_day &starting_date,
                           const date::year_month_day &ending_date) {
  set_current_time(0);
  calendar_date_ = date::sys_days(starting_date);
  world_events_.initialize();
}

void Scheduler::run() {
  current_time_ = 0;
  for (current_time_ = 0; !can_stop(); current_time_++) {
    if (current_time_ % Model::get_config()->get_model_settings().get_days_between_stdout_output()
        == 0) {
      spdlog::info("Day: {}", current_time_);
    }
    begin_time_step();
    daily_update();
    end_time_step();
    calendar_date_ += date::days{1};
  }
}

void Scheduler::begin_time_step() {
  if (Model::get_instance() != nullptr) { Model::get_instance()->begin_time_step(); }
}

void Scheduler::daily_update() {
  if (Model::get_instance() != nullptr) {
    Model::get_instance()->daily_update();

    if (is_today_first_day_of_month()) { Model::get_instance()->monthly_update(); }

    if (is_today_first_day_of_year()) { Model::get_instance()->yearly_update(); }

    // Execute world/population events
    world_events_.execute_events(current_time_);

    // Update individual events through the population
    Model::get_population()->execute_all_individual_events(current_time_);
  }
}

void Scheduler::end_time_step() {
  if (Model::get_instance() != nullptr) { Model::get_instance()->end_time_step(); }
}

bool Scheduler::can_stop() const {
  return current_time_ > Model::get_config()->get_simulation_timeframe().get_total_time();
}

int Scheduler::get_current_day_in_year() { return TimeHelpers::day_of_year(calendar_date_); }

unsigned int Scheduler::get_current_month_in_year() {
  return TimeHelpers::month_of_year(calendar_date_);
}

bool Scheduler::is_today_last_day_of_year() {
  date::year_month_day ymd{calendar_date_};
  return ymd.month() == date::month{12} && ymd.day() == date::day{31};
}

bool Scheduler::is_today_first_day_of_month() {
  date::year_month_day ymd{calendar_date_};
  return ymd.day() == date::day{1};
}

bool Scheduler::is_today_first_day_of_year() {
  date::year_month_day ymd{calendar_date_};
  return ymd.month() == date::month{1} && ymd.day() == date::day{1};
}

bool Scheduler::is_today_last_day_of_month() {
  const auto next_date = calendar_date_ + date::days{1};
  date::year_month_day ymd{next_date};
  return ymd.day() == date::day{1};
}

int Scheduler::get_days_to_next_year() const {
  return TimeHelpers::number_of_days_to_next_year(calendar_date_);
}

int Scheduler::get_days_to_next_n_year(int n) const {
  auto from_ymd = date::year_month_day{calendar_date_};
  auto to_ymd = from_ymd + date::years{n};
  return TimeHelpers::days_between(from_ymd, to_ymd);
}

unsigned Scheduler::get_days_in_current_month() const {
  auto date = static_cast<date::year_month_day>(calendar_date_);
  return TimeHelpers::days_in_month(static_cast<int>(date.year()),
                                    static_cast<unsigned int>(date.month()));
}

int Scheduler::get_current_year() const {
  auto date = static_cast<date::year_month_day>(calendar_date_);
  return static_cast<int>(date.year());
}

unsigned int Scheduler::get_current_day_of_month() const {
  auto date = static_cast<date::year_month_day>(calendar_date_);
  return static_cast<unsigned int>(date.day());
}

date::year_month_day Scheduler::get_ymd_after_months(int months) const {
  return date::year_month_day(calendar_date_) + date::months(months);
}

int Scheduler::get_days_to_ymd(const date::year_month_day &ymd) const {
  return TimeHelpers::days_between(calendar_date_, ymd);
}

date::year_month_day Scheduler::get_ymd_after_days(int days) const {
  return calendar_date_ + date::days(days);
}

std::time_t Scheduler::get_unix_time() const {
  return std::chrono::system_clock::to_time_t(calendar_date_);
}

date::year_month_day Scheduler::get_calendar_date() const { return calendar_date_; }

