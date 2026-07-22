#ifndef MMCREPORTER_H
#define MMCREPORTER_H

#include <sstream>

#include "Reporter.h"

class PersonIndexByLocationStateAgeClass;

class MMCReporter : public Reporter {
public:
  // Disallow copy
  MMCReporter(const MMCReporter &) = delete;
  MMCReporter &operator=(const MMCReporter &) = delete;

  // Disallow move
  MMCReporter(MMCReporter &&) = delete;
  MMCReporter &operator=(MMCReporter &&) = delete;

  MMCReporter();

  ~MMCReporter() override = default;

  void initialize(int job_number, const std::string &path) override;

  void before_run() override;

  void after_run() override;

  void begin_time_step() override;

  void print_treatment_failure_rate_by_therapy();

  void print_ntf_by_location();

  void print_genotype_frequency();

  void monthly_report() override;

  void print_eir_pfpr_by_location();

  std::stringstream ss;
  const std::string group_sep = "-1111\t";
  const std::string sep = "\t";
};

#endif  // MMCREPORTER_H
