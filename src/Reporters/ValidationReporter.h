#ifndef VALIDATIONREPORTER_H
#define VALIDATIONREPORTER_H

#include <spdlog/sinks/basic_file_sink.h>
#include <spdlog/spdlog.h>

#include <memory>

#include "Reporter.h"

class ValidationReporter : public Reporter {
public:
  // Disable copy, assignment, and move
  ValidationReporter(const ValidationReporter &) = delete;
  ValidationReporter &operator=(const ValidationReporter &) = delete;
  ValidationReporter(ValidationReporter &&) = delete;
  ValidationReporter &operator=(ValidationReporter &&) = delete;

  std::shared_ptr<spdlog::logger> monthly_data_logger;
  std::shared_ptr<spdlog::logger> summary_data_logger;
  std::shared_ptr<spdlog::logger> gene_db_logger;
  std::shared_ptr<spdlog::logger> gene_freq_logger;
  std::shared_ptr<spdlog::logger> monthly_mutation_logger;
  std::shared_ptr<spdlog::logger> mosquito_res_count_logger;

  ValidationReporter();
  ~ValidationReporter() override = default;

  void initialize(int job_number, const std::string &path) override;
  void before_run() override;
  void after_run() override;
  void begin_time_step() override;
  void monthly_report() override;
  void print_eir_pfpr_by_location(std::stringstream &ss);

  std::string monthly_mutation_path;
  std::string mosquito_res_count_path;
};

#endif  // VALIDATIONREPORTER_H
