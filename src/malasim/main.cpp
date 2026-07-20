#include "Simulation/Model.h"
#include "Utils/Cli.h"
#include "Utils/Logger.h"
#include "spdlog/spdlog.h"
#include "version.h"

#include "Population/ClonalParasitePopulation.h"
#include "Population/Person/Person.h"
#include "Population/Population.h"

#include <sys/resource.h>

double get_memory_kb() {
  struct rusage usage{};
  getrusage(RUSAGE_SELF, &usage);
#ifdef __APPLE__
  return usage.ru_maxrss / 1024.0;
#else
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

  if (Model::get_cli_input().print_memory_stats) {
    std::cout << "sizeof(Person) = " << sizeof(Person) << '\n';
    std::cout << "sizeof(ClonalParasitePopulation) = " << sizeof(ClonalParasitePopulation)
              << '\n';
    std::cout << "sizeof(Population) = " << sizeof(Population) << '\n';
    spdlog::drop_all();
    return 0;
  }

  if (Model::get_instance()->initialize()) {
    Model::get_instance()->run();

    const double memory_kb = get_memory_kb();

    std::cout << std::fixed << std::setprecision(2);

    if (memory_kb >= 1024.0 * 1024.0) {
      const double memory_gb = memory_kb / (1024.0 * 1024.0);
      std::cout << "Memory Usage: " << memory_gb << " GB\n";
    } else if (memory_kb >= 1024.0) {
      const double memory_mb = memory_kb / 1024.0;
      std::cout << "Memory Usage: " << memory_mb << " MB\n";
    } else {
      std::cout << "Memory Usage: " << memory_kb << " KB\n";
    }

    Model::get_instance()->release();
  } else {
    spdlog::get("default_logger")->error("Model initialization failed.");
  }
  spdlog::drop_all();
  return 0;
}
