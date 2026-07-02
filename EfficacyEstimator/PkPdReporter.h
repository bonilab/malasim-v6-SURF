/* 
 * File:   PkPdReporter.h
 * Author: Merlin
 *
 * Created on October 29, 2014, 12:56 PM
 */

#ifndef PKPDREPORTER_H
#define PKPDREPORTER_H

#include "Reporters/Reporter.h"
#include "Utils/Cli.h"
#include "Utils/TypeDef.h"

#include <fstream>
#include <sstream>

class PkPdReporter : public Reporter {
  DoubleVector yesterday_density;

public:
  std::stringstream ss;
  const std::string group_sep = "-1111\t";
  const std::string sep = ",";
  std::string prefix;

  explicit PkPdReporter(utils::DxGAppInput* appInput = nullptr);

  ~PkPdReporter() override;

  void initialize(int job_number, const std::string& path) override;

  void before_run() override;

  void after_run() override;

  void begin_time_step() override;

  void after_time_step();

  void monthly_report() override;

private:
  [[nodiscard]] bool is_recurrence_test() const;

  utils::DxGAppInput* appInput{nullptr};
  std::ofstream outputFStream;
};

#endif /* PKPDREPORTER_H */
