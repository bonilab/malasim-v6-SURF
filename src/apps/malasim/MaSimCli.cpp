#include "apps/malasim/MaSimCli.h"

#include <spdlog/spdlog.h>

#include <fstream>

namespace utils {

void MaSimCli::add_options(CLI::App &app, Input &input) {
  app.add_option("-i,--input", input.input_path, "Input filename. Default: `input.yml`.");
  app.add_option("-o,--output", input.output_path, "Output path. Default: `./`.");
  app.add_option("-r,--reporter", input.reporter,
                 "Reporter type(s), comma separated, e.g. "
                 "`-r SQLiteMonthlyReporter,SQLiteValidationReporter`. "
                 "Default: `SQLiteMonthlyReporter`.");
  app.add_option("-v,--verbosity", input.verbosity,
                 "Sets the verbosity of the logging. Default: 0");
  app.add_option("-j,--job", input.job_number,
                 "Sets the study to associate with the configuration (or database id). Default: 0");
  app.add_flag("-d,--dump", input.dump_movement_matrix, "Dump the movement matrix as calculated.");
  app.add_flag("-l,--list", input.list_reporters, "List the possible reporters and exit.");
  app.add_flag("--im", input.record_individual_movement, "Record individual movement data.");
  app.add_flag("--mc", input.record_cell_movement, "Record the movement between cells.");
  app.add_flag("--md", input.record_district_movement, "Record the movement between districts.");
  app.add_option("--replicate", input.replicate, "Replicate number. Default: 1");
  app.add_flag("--memory-stats", input.print_memory_stats,
               "Print memory statistics for key classes and exit.");
}

bool MaSimCli::validate(Input &input) {
  if (input.list_reporters || input.print_memory_stats) { return true; }

  std::ifstream input_file(input.input_path);
  if (!input_file.good()) {
    spdlog::error("Err: Input file error or not found!");
    return false;
  }

  if (input.record_cell_movement && input.record_district_movement) {
    spdlog::error("--mc and --md are mutual exclusive and may not be run together.");
    return false;
  }

  input.record_movement = input.record_individual_movement || input.record_cell_movement
                          || input.record_district_movement;
  if (input.record_movement) { spdlog::info("Movement data will be recorded."); }

  switch (input.verbosity) {
    case 1:
      spdlog::set_level(spdlog::level::debug);
      spdlog::info("Verbosity level set to 1. Info, debug will be logged.");
      break;
    case 2:
      spdlog::set_level(spdlog::level::trace);
      spdlog::info("Verbosity level set to 2. Info, debug and trace will be logged.");
      break;
    case 0:
    default:
      spdlog::set_level(spdlog::level::info);
      spdlog::info("Verbosity level set to 0. Only info will be logged.");
      break;
  }

  return true;
}

}  // namespace utils
