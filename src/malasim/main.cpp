#include "Simulation/Model.h"
#include "Utils/Cli.h"
#include "Utils/Logger.h"
#include "spdlog/spdlog.h"
#include "version.h"

double get_memory_kb() {
  struct rusage usage{};
  getrusage(RUSAGE_SELF, &usage);
#ifdef __APPLE__
  // macOS: ru_maxrss in bytes
  return usage.ru_maxrss / 1024.0;
#else
  // Linux: ru_maxrss in kilobytes
  return usage.ru_maxrss;
#endif
}

int main(int argc, char** argv) {
  Logger::initialize(spdlog::level::info);
  spdlog::info("Malaria Simulation v{}", VERSION);
  spdlog::info("Starting...");
  try {
    auto cli_input = utils::Cli::parse_args(argc, argv);
    Model::set_cli_input(std::move(cli_input));
  } catch (...) {
    spdlog::error("Argument parsing failed. Exiting.");
    return 1;
  }
  if (Model::get_instance()->initialize()) {
    Model::get_instance()->run();

    double memory_kb = get_memory_kb();
    std::cout << "Memory Usage: " << memory_kb << " KB" << '\n';
    Model::get_instance()->release();
  } else {
    spdlog::get("default_logger")->error("Model initialization failed.");
  }
  spdlog::drop_all();
  return 0;
}
