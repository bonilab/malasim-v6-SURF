#ifndef POMS_BFREPORTER_H
#define POMS_BFREPORTER_H

#include <spdlog/sinks/basic_file_sink.h>
#include <spdlog/spdlog.h>

#include "Reporter.h"

class MonthlyReporter : public Reporter {
public:
  // Disable copy, assignment, and move
  MonthlyReporter(const MonthlyReporter &) = delete;
  MonthlyReporter &operator=(const MonthlyReporter &) = delete;
  MonthlyReporter(MonthlyReporter &&) = delete;
  MonthlyReporter &operator=(MonthlyReporter &&) = delete;

  std::shared_ptr<spdlog::logger> gene_freq_logger;
  std::shared_ptr<spdlog::logger> monthly_data_logger;
  std::shared_ptr<spdlog::logger> summary_data_logger;
  std::shared_ptr<spdlog::logger> gene_db_logger;

  MonthlyReporter();
  ~MonthlyReporter() override = default;

  void initialize(int job_number, const std::string &path) override;
  void before_run() override;
  void after_run() override;
  void begin_time_step() override;
  void monthly_report() override;
  void print_eir_pfpr_by_location(std::stringstream &ss);
};

#endif  // POMS_BFREPORTER_H
