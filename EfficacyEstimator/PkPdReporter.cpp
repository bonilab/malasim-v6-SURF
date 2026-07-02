/* 
 * File:   PkPdReporter.cpp
 * Author: Merlin
 * 
 * Created on October 29, 2014, 12:56 PM
 */

#include "PkPdReporter.h"

#include "Population/ImmuneSystem/ImmuneSystem.h"
#include "Population/Person/Person.h"
#include "Simulation/Model.h"
#include "Utils/Index/PersonIndexAll.h"

PkPdReporter::PkPdReporter(utils::DxGAppInput* appInput) : appInput{appInput} {
  if (appInput && !appInput->output_file.empty() && !appInput->is_recurrence_test) {
    outputFStream.open(appInput->output_file);
  }
}

PkPdReporter::~PkPdReporter() {
  if (outputFStream.is_open()) {
    outputFStream.close();
  }
}

void PkPdReporter::initialize(int /*job_number*/, const std::string& path) {
  prefix = path;
}

void PkPdReporter::before_run() {
  if (!is_recurrence_test()) {
    return;
  }

  const std::string out_path =
      appInput && !appInput->output_file.empty() ? appInput->output_file
                                                 : fmt::format("{}_parasitaemia.csv", prefix);
  outputFStream.open(out_path, std::ios::out);
  outputFStream << "time,individual,recrudescence,parasitaemia" << '\n';
}

void PkPdReporter::begin_time_step() {
  if (is_recurrence_test()) {
    Model::get_mdc()->perform_population_statistic();
    Model::get_mdc()->blood_slide_prevalence_by_location()[0] = 0.1;

    if ((Model::get_scheduler()->current_time()
         % Model::get_config()->get_model_settings().get_days_between_stdout_output() == 0)
        && (Model::get_scheduler()->current_time()
            > Model::get_config()->get_simulation_timeframe().get_start_collect_data_day())) {
      auto current_time = Model::get_scheduler()->current_time();

      for (int i = 0; i < Model::get_population()->all_persons()->v_person().size(); i++) {
        auto* person = Model::get_population()->all_persons()->v_person()[i].get();
        const auto recrudescence_state = person->get_recurrence_status();

        double parasitaemia = 0.0;
        if (person->get_all_clonal_parasite_populations()->size() >= 1) {
          parasitaemia =
              person->get_all_clonal_parasite_populations()->at(0)->get_log10_infectious_density();
        } else {
          parasitaemia = Model::get_config()
                             ->get_parasite_parameters()
                             .get_parasite_density_levels()
                             .get_log_parasite_density_cured();
        }

        outputFStream << current_time << "," << i << ","
                      << static_cast<int>(recrudescence_state) << "," << parasitaemia << '\n';
      }
    }
    return;
  }

  ss << Model::get_scheduler()->current_time();

  for (int i = 0; i < Model::get_population()->all_persons()->v_person().size(); i++) {
    Person* p_person = Model::get_population()->all_persons()->v_person()[i].get();
    if (p_person->get_all_clonal_parasite_populations()->size() > 0) {
      ss << sep
         << p_person->get_all_clonal_parasite_populations()->at(0)->get_log10_infectious_density();
    } else {
      ss << sep << Model::get_config()
                       ->get_parasite_parameters()
                       .get_parasite_density_levels()
                       .get_log_parasite_density_cured();
    }
    if (appInput && appInput->is_print_immunity_level) {
      ss << sep << p_person->get_immune_system()->get_latest_immune_value();
    }
  }

  if (outputFStream.is_open()) {
    outputFStream << ss.str() << '\n';
  }
  ss.str("");
}

void PkPdReporter::after_time_step() {}

void PkPdReporter::monthly_report() {}

void PkPdReporter::after_run() {
  Model::get_mdc()->update_after_run();
  if (outputFStream.is_open()) {
    outputFStream.close();
  }
}

bool PkPdReporter::is_recurrence_test() const {
  return appInput && appInput->is_recurrence_test;
}
