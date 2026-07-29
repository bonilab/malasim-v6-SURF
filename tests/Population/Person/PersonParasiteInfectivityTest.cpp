#include "PersonParasiteTestBase.h"
#include "Population/SingleHostClonalParasitePopulations.h"

TEST_F(PersonParasiteTest, RelativeInfectivity) {
  EXPECT_DOUBLE_EQ(
      Person::relative_infectivity(ClonalParasitePopulation::LOG_ZERO_PARASITE_DENSITY), 0.0);
}

TEST_F(PersonParasiteTest, RelativeInfectivityWithLowDensity) {
  EXPECT_CALL(*mock_random_, cdf_standard_normal_distribution(::testing::_))
      .WillOnce(::testing::Return(0.4));
  EXPECT_DOUBLE_EQ(Person::relative_infectivity(1.0), 0.17);
}

TEST_F(PersonParasiteTest, RelativeInfectivityWithMediumDensity) {
  EXPECT_CALL(*mock_random_, cdf_standard_normal_distribution(::testing::_))
      .WillOnce(::testing::Return(0.7));
  EXPECT_DOUBLE_EQ(Person::relative_infectivity(3.0), 0.5);
}

TEST_F(PersonParasiteTest, RelativeInfectivityWithHighDensity) {
  EXPECT_CALL(*mock_random_, cdf_standard_normal_distribution(::testing::_))
      .WillOnce(::testing::Return(0.95));
  const double infectivity = Person::relative_infectivity(5.0);
  EXPECT_DOUBLE_EQ(infectivity, 0.9125);
  EXPECT_GT(infectivity, 0.9);
}

TEST_F(PersonParasiteTest, RelativeInfectivityMonotonicIncreasing) {
  EXPECT_CALL(*mock_random_, cdf_standard_normal_distribution(::testing::_))
      .WillOnce(::testing::Return(0.3))
      .WillOnce(::testing::Return(0.6))
      .WillOnce(::testing::Return(0.9));

  const double infectivity1 = Person::relative_infectivity(-2.0);
  const double infectivity2 = Person::relative_infectivity(0.0);
  const double infectivity3 = Person::relative_infectivity(2.0);

  EXPECT_LT(infectivity1, infectivity2);
  EXPECT_LT(infectivity2, infectivity3);
  EXPECT_DOUBLE_EQ(infectivity1, 0.1);
  EXPECT_DOUBLE_EQ(infectivity2, 0.37);
  EXPECT_DOUBLE_EQ(infectivity3, 0.82);
}

TEST_F(PersonParasiteTest, RelativeInfectivityFormulaVerification) {
  constexpr double cdf_value = 0.5;
  EXPECT_CALL(*mock_random_, cdf_standard_normal_distribution(::testing::_))
      .WillOnce(::testing::Return(cdf_value));
  EXPECT_DOUBLE_EQ(Person::relative_infectivity(1.0), (cdf_value * cdf_value) + 0.01);
}

TEST_F(PersonParasiteTest, RelativeInfectivityFormulaWithDifferentCDFValues) {
  for (double cdf : {0.1, 0.3, 0.5, 0.7, 0.9, 0.99}) {
    EXPECT_CALL(*mock_random_, cdf_standard_normal_distribution(::testing::_))
        .WillOnce(::testing::Return(cdf));
    EXPECT_DOUBLE_EQ(Person::relative_infectivity(2.0), (cdf * cdf) + 0.01);
  }
}

TEST_F(PersonParasiteTest, RelativeInfectivityDnCalculation) {
  constexpr double log10_density = 2.0;
  constexpr double expected_d_n = (log10_density * 3.91) + 0.00031;
  EXPECT_CALL(*mock_random_, cdf_standard_normal_distribution(::testing::DoubleEq(expected_d_n)))
      .WillOnce(::testing::Return(0.8));
  EXPECT_DOUBLE_EQ(Person::relative_infectivity(log10_density), (0.8 * 0.8) + 0.01);
}

TEST_F(PersonParasiteTest, RelativeInfectivityWithCustomRoStar) {
  auto &epi_params = mock_config_->get_epidemiological_parameters();
  auto rel_inf = epi_params.get_relative_infectivity();
  rel_inf.set_ro_star(1.0);
  epi_params.set_relative_infectivity(rel_inf);

  EXPECT_CALL(*mock_random_, cdf_standard_normal_distribution(::testing::DoubleEq(4.91)))
      .WillOnce(::testing::Return(0.6));
  EXPECT_DOUBLE_EQ(Person::relative_infectivity(1.0), 0.37);
}

TEST_F(PersonParasiteTest, RelativeInfectivityWithCustomSigma) {
  auto &epi_params = mock_config_->get_epidemiological_parameters();
  auto rel_inf = epi_params.get_relative_infectivity();
  rel_inf.set_sigma(2.0);
  epi_params.set_relative_infectivity(rel_inf);

  EXPECT_CALL(*mock_random_, cdf_standard_normal_distribution(::testing::DoubleEq(6.00031)))
      .WillOnce(::testing::Return(0.95));
  EXPECT_DOUBLE_EQ(Person::relative_infectivity(3.0), 0.9125);
}

TEST_F(PersonParasiteTest, RelativeInfectivityMinimumValue) {
  EXPECT_CALL(*mock_random_, cdf_standard_normal_distribution(::testing::_))
      .WillOnce(::testing::Return(0.0));
  EXPECT_DOUBLE_EQ(Person::relative_infectivity(1.0), 0.01);
}

TEST_F(PersonParasiteTest, RelativeInfectivityMaximumValue) {
  EXPECT_CALL(*mock_random_, cdf_standard_normal_distribution(::testing::_))
      .WillOnce(::testing::Return(1.0));
  EXPECT_DOUBLE_EQ(Person::relative_infectivity(1.0), 1.0);
}
