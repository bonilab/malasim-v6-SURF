/*
 * AgeBandReporter.h
 *
 * This reporter is intended to be used during model calibration and validation
 * and reports the age-banded PfPR during the last year of the given simulation
 * to an SQLite database.
 */
#ifndef AGEBANDREPORTER_H
#define AGEBANDREPORTER_H

#include <memory>
#include <sstream>
#include <vector>

#include "Reporters/Reporter.h"

namespace spdlog {
class logger;
}

class AgeBandReporter : public Reporter {
public:
  // Disallow copy
  AgeBandReporter(const AgeBandReporter &) = delete;
  AgeBandReporter &operator=(const AgeBandReporter &) = delete;

  // Disallow move
  AgeBandReporter(AgeBandReporter &&) = delete;
  AgeBandReporter &operator=(AgeBandReporter &&) = delete;

  AgeBandReporter() = default;

  ~AgeBandReporter() override = default;

  void before_run() override {}

  void begin_time_step() override {}

  void initialize(int job_number, const std::string &path) override;

  void after_run() override {}

  void monthly_report() override;

private:
  // When to start logging the data

  int start_recording_ = -1;

  // String streams to use when working with the loggers
  std::stringstream pfpr_;
  std::stringstream cases_;

  // Mapping of the locations to their districts
  std::vector<int> district_lookup_;

  std::shared_ptr<spdlog::logger> pfpr_logger_;
  std::shared_ptr<spdlog::logger> cases_logger_;
};

#endif
