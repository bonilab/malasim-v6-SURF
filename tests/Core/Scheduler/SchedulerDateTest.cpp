#include <gtest/gtest.h>

#include "Core/Scheduler/Scheduler.h"
#include "date/date.h"

// Covers the pure calendar/date helpers of Scheduler. These do not touch the
// Model singleton, so they can run without any configuration or fixtures.
class SchedulerDateTest : public ::testing::Test {
protected:
  void init(const date::year_month_day &start) { scheduler_.initialize(start, start); }

  Scheduler scheduler_;
};

TEST_F(SchedulerDateTest, InitializeSetsCurrentTimeAndCalendarDate) {
  const date::year_month_day input_start{date::year{2020}, date::month{3}, date::day{15}};
  init(input_start);

  EXPECT_EQ(scheduler_.current_time(), 0);
  EXPECT_EQ(scheduler_.get_calendar_date(), input_start);
}

TEST_F(SchedulerDateTest, FirstDayOfMonthDetected) {
  init(date::year_month_day{date::year{2020}, date::month{6}, date::day{1}});
  EXPECT_TRUE(scheduler_.is_today_first_day_of_month());
  EXPECT_FALSE(scheduler_.is_today_last_day_of_month());
}

TEST_F(SchedulerDateTest, MidMonthIsNotFirstDay) {
  init(date::year_month_day{date::year{2020}, date::month{6}, date::day{15}});
  EXPECT_FALSE(scheduler_.is_today_first_day_of_month());
  EXPECT_FALSE(scheduler_.is_today_first_day_of_year());
}

TEST_F(SchedulerDateTest, LastDayOfMonthDetected) {
  init(date::year_month_day{date::year{2020}, date::month{6}, date::day{30}});
  EXPECT_TRUE(scheduler_.is_today_last_day_of_month());
  EXPECT_FALSE(scheduler_.is_today_first_day_of_month());
}

TEST_F(SchedulerDateTest, LastDayOfYearDetected) {
  init(date::year_month_day{date::year{2020}, date::month{12}, date::day{31}});
  EXPECT_TRUE(scheduler_.is_today_last_day_of_year());
  EXPECT_TRUE(scheduler_.is_today_last_day_of_month());
}

TEST_F(SchedulerDateTest, FirstDayOfYearDetected) {
  init(date::year_month_day{date::year{2021}, date::month{1}, date::day{1}});
  EXPECT_TRUE(scheduler_.is_today_first_day_of_year());
  EXPECT_TRUE(scheduler_.is_today_first_day_of_month());
  EXPECT_FALSE(scheduler_.is_today_last_day_of_year());
}

TEST_F(SchedulerDateTest, CurrentYearMonthDayAccessors) {
  init(date::year_month_day{date::year{2020}, date::month{2}, date::day{20}});
  EXPECT_EQ(scheduler_.get_current_year(), 2020);
  EXPECT_EQ(scheduler_.get_current_day_of_month(), 20u);
  EXPECT_EQ(scheduler_.get_current_month_in_year(), 2u);
}

TEST_F(SchedulerDateTest, DaysInCurrentMonthHandlesLeapFebruary) {
  init(date::year_month_day{date::year{2020}, date::month{2}, date::day{1}});
  EXPECT_EQ(scheduler_.get_days_in_current_month(), 29u);  // 2020 leap year

  init(date::year_month_day{date::year{2021}, date::month{2}, date::day{1}});
  EXPECT_EQ(scheduler_.get_days_in_current_month(), 28u);  // non-leap
}

TEST_F(SchedulerDateTest, DaysToNextYear) {
  init(date::year_month_day{date::year{2020}, date::month{12}, date::day{31}});
  EXPECT_EQ(scheduler_.get_days_until_next_year_anniversary(), 365);

  init(date::year_month_day{date::year{2021}, date::month{1}, date::day{1}});
  EXPECT_EQ(scheduler_.get_days_until_next_year_anniversary(), 365);
}

TEST_F(SchedulerDateTest, DaysToNextNYears) {
  init(date::year_month_day{date::year{2020}, date::month{1}, date::day{1}});
  EXPECT_EQ(scheduler_.get_days_until_n_years_from_now(1), 366);  // 2020 is leap
  EXPECT_EQ(scheduler_.get_days_until_n_years_from_now(0), 0);
}

TEST_F(SchedulerDateTest, DayInYearAccessor) {
  init(date::year_month_day{date::year{2021}, date::month{1}, date::day{1}});
  EXPECT_EQ(scheduler_.get_current_day_in_year(), 1);

  init(date::year_month_day{date::year{2021}, date::month{12}, date::day{31}});
  EXPECT_EQ(scheduler_.get_current_day_in_year(), 365);
}

TEST_F(SchedulerDateTest, YmdAfterMonthsAndDays) {
  init(date::year_month_day{date::year{2020}, date::month{1}, date::day{15}});
  const date::year_month_day expected_month{date::year{2020}, date::month{2}, date::day{15}};
  EXPECT_EQ(scheduler_.get_ymd_after_months(1), expected_month);

  const date::year_month_day expected_day{date::year{2020}, date::month{1}, date::day{16}};
  EXPECT_EQ(scheduler_.get_ymd_after_days(1), expected_day);
}

TEST_F(SchedulerDateTest, DaysToYmdIsSigned) {
  init(date::year_month_day{date::year{2020}, date::month{1}, date::day{10}});
  const date::year_month_day later{date::year{2020}, date::month{1}, date::day{20}};
  const date::year_month_day earlier{date::year{2020}, date::month{1}, date::day{5}};
  EXPECT_EQ(scheduler_.get_days_to_ymd(later), 10);
  EXPECT_EQ(scheduler_.get_days_to_ymd(earlier), -5);
}

TEST_F(SchedulerDateTest, UnixTimeMatchesCalendarDate) {
  const date::year_month_day input_start{date::year{2020}, date::month{1}, date::day{1}};
  init(input_start);
  const auto expected = std::chrono::system_clock::to_time_t(date::sys_days(input_start));
  EXPECT_EQ(scheduler_.get_unix_time(), expected);
}

TEST_F(SchedulerDateTest, ForceStopFlagRoundTrips) {
  EXPECT_FALSE(scheduler_.is_force_stop());
  scheduler_.set_is_force_stop(true);
  EXPECT_TRUE(scheduler_.is_force_stop());
}

TEST_F(SchedulerDateTest, CurrentTimeRoundTrips) {
  scheduler_.set_current_time(42);
  EXPECT_EQ(scheduler_.current_time(), 42);
}
