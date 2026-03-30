#include <spdlog/spdlog-inl.h>

#include "gtest/gtest.h"

class GlobalEnvironment : public ::testing::Environment {
public:
  void SetUp() override {
    // Run once before all tests
    std::cout << "Global SetUp\n";
    spdlog::set_level(spdlog::level::off);
  }

  void TearDown() override {
    // Run once after all tests
    std::cout << "Global TearDown\n";
    // e.g., shutdown services, flush logs
  }
};

int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  ::testing::AddGlobalTestEnvironment(new GlobalEnvironment);
  // Create a console logger used by the whole test suite
  // auto console = spdlog::stderr_color_mt("gtest");
  // console->set_pattern("[%Y-%m-%d %T.%e] [%n] [%^%l%$] %v");
  // console->set_level(log_level);
  //
  // // Optional: set as default logger so `SPDLOG_LOGGER_*` macros use it
  // spdlog::set_default_logger(console);

  return RUN_ALL_TESTS();
}
