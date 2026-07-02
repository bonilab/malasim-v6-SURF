/*
* File:   PkPdReporter.h
 * Author: Merlin
 *
 * Created on October 29, 2014, 12:56 PM
 */

#ifndef PKPDREPORTER_H
#define    PKPDREPORTER_H

#include "Reporters/Reporter.h"
#include <fstream>
#include "Utils/Cli.h"

namespace utils {
  class Cli;
}

class PkPdReporter : public Reporter {

public:
  std::string prefix;
  std::ofstream parasitaemia_file;

  PkPdReporter(utils::Cli::DxGAppInput* appInput=nullptr);

  virtual ~PkPdReporter();

  void initialize(int job_number, const std::string& path) override;

  void before_run() override;

  void after_run() override;

  void begin_time_step() override;

  virtual void after_time_step();

  void monthly_report() override;

private:
  utils::Cli::DxGAppInput* appInput{nullptr};
};

#endif    /* PKPDREPORTER_H */
