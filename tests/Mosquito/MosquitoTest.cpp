#include <gtest/gtest.h>

#include <algorithm>

#include "Configuration/Config.h"
#include "Utils/Random.h"
#include "Mosquito/Mosquito.h"
#include "Population/Person/Person.h"

class MosquitoTest : public ::testing::Test {
protected:
  void SetUp() override {}
  // void TearDown() override {}
};


TEST_F(MosquitoTest, PrmcSample) {
  Mosquito m;
  utils::Random r;
  r.set_seed(1);

  int prmc_size { 10 };
  int pop_size { 100 };

  // repeat the process 10 times
  for (int n = 0; n < 10; ++n) {
    std::vector<double> foi_distribution(pop_size, 1.0);

    // even id person will have no selection
    for (int i = 0; i < pop_size; i += 2) {
      foi_distribution[i] = 0;
    }

    std::vector<std::unique_ptr<Person>> all_person;
    std::vector<Person*> all_person_ptr;
    for (int i = 0; i < pop_size; ++i) {
      auto p = std::make_unique<Person>();
      p->set_last_therapy_id(i);
      all_person.push_back(std::move(p));
      all_person_ptr.push_back(all_person.back().get());
    }

    auto samples = r.multinomial_sampling<Person>(prmc_size, foi_distribution, all_person_ptr, true);

    EXPECT_EQ(samples.size(), prmc_size);

    for (int i = 0; i < prmc_size; ++i) {
      EXPECT_EQ(samples[i]->get_last_therapy_id() % 2, 1)
          << fmt::format("failed with p_id: {}", samples[i]->get_last_therapy_id());
    }

    all_person.clear();
    all_person_ptr.clear();
  }
}

TEST_F(MosquitoTest, PrmcSampleWhenThereIsNoFOI) {
  Mosquito m;
  utils::Random r;
  r.set_seed(1);

  int prmc_size { 10 };
  int pop_size { 100 };

  // repeat the process 10 times
  for (int n = 0; n < 1; ++n) {
    std::vector<double> foi_distribution(pop_size, 0);

    std::vector<std::unique_ptr<Person>> all_person;
    std::vector<Person*> all_person_ptr;
    for (int i = 0; i < pop_size; ++i) {
      auto p = std::make_unique<Person>();
      p->set_last_therapy_id(i);
      all_person.push_back(std::move(p));
      all_person_ptr.push_back(all_person.back().get());
    }

    auto samples = r.multinomial_sampling<Person>(prmc_size, foi_distribution, all_person_ptr, true);
    EXPECT_EQ(samples.size(), prmc_size);

    for (int i = 0; i < prmc_size; ++i) {
      EXPECT_EQ(samples[i], nullptr);
    }

    all_person.clear();
    all_person_ptr.clear();
  }
}

TEST_F(MosquitoTest, BuildsInterruptedFeedingIndicesWithinRequestedSize) {
  Mosquito mosquito;
  utils::Random random;
  random.set_seed(7);
  const auto indices = Mosquito::build_interrupted_feeding_indices(&random, 0.25, 20);
  ASSERT_EQ(indices.size(), 20u);
  EXPECT_LE(std::count(indices.begin(), indices.end(), 1U), 20);
}

TEST_F(MosquitoTest, ConvertsNewGenotypeStringsToLegacyFormats) {
  const std::string genotype = "aaaa|bbbb|cccc|dddd|eeee|ffff|gggggggX|hhhh|iiii|jjjj|kkkk|llll|mmmmmmmmmmmmm|n|o";
  EXPECT_EQ(Mosquito::get_old_genotype_string(genotype), "eeee|ggggggg|mmmmmmmmmmmmm|n");
  EXPECT_EQ(Mosquito::get_old_genotype_string2(genotype), "geeeemn");
}
