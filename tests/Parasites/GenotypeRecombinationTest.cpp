#include <gtest/gtest.h>

#include "Parasites/Genotype.h"
#include "Simulation/Model.h"
#include "Utils/Cli.h"
#include "fixtures/TestFileGenerators.h"

namespace {
class FixedUniformRandom final : public utils::Random {
 public:
  explicit FixedUniformRandom(double value) : utils::Random(nullptr, 1234), value_(value) {}
  double random_uniform() override { return value_; }

 private:
  double value_;
};
}  // namespace

class GenotypeRecombinationTest : public ::testing::Test {
 protected:
  void SetUp() override {
    test_fixtures::setup_test_environment();
    utils::Cli::MaSimAppInput cli_input;
    cli_input.input_path = "test_input.yml";
    Model::set_cli_input(cli_input);
    ASSERT_TRUE(Model::get_instance()->initialize());
  }

  void TearDown() override {
    Model::get_instance()->release();
    test_fixtures::cleanup_test_files();
  }
};

TEST_F(GenotypeRecombinationTest, FreeRecombinationProducesDatabaseGenotype) {
  const auto aa_sequence = Model::get_genotype_db()->at(0)->get_aa_sequence();
  Genotype female(aa_sequence);
  Genotype male(aa_sequence);
  FixedUniformRandom random(0.25);

  auto* result = Genotype::free_recombine(Model::get_config(), &random, &female, &male);

  ASSERT_NE(result, nullptr);
  EXPECT_EQ(result->get_aa_sequence(), aa_sequence);
}
