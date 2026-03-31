/* 
 * File:   PkPdReporter.h
 * Author: Merlin
 *
 * Created on October 29, 2014, 12:56 PM
 */

#ifndef PKPDREPORTER_H
#define    PKPDREPORTER_H

#include "Reporters/Reporter.h"
#include <sstream>
#include <fstream>

#include "Utils/Cli.h"
#include "Utils/TypeDef.h"

namespace utils {
class Cli;
}

class PkPdReporter : public Reporter {

 DoubleVector yesterday_density;

 public:
    std::stringstream ss;
    const std::string group_sep = "-1111\t";
    const std::string sep = ",";

    PkPdReporter(utils::DxGAppInput* appInput=nullptr);

  virtual ~PkPdReporter();

  void initialize(int job_number, const std::string& path) override{}

  void before_run() override;

  void after_run() override;

  void begin_time_step() override;

  virtual void after_time_step();

  void monthly_report() override;

 private:
    utils::DxGAppInput* appInput{nullptr};
    std::ofstream outputFStream;

};

#endif    /* PKPDREPORTER_H */
