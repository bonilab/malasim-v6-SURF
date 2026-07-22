#include <gtest/gtest.h>
#include <yaml-cpp/yaml.h>
#include "Configuration/ParasiteParameters.h"
#include "Utils/Constants.h"

TEST(ParasiteParametersDefaultsTest, UsesNamedCuredDensityConstant) {
    const ParasiteParameters::ParasiteDensityLevels density_levels;
    EXPECT_DOUBLE_EQ(density_levels.get_log_parasite_density_cured(),
                     Constants::DEFAULT_LOG10_PARASITE_DENSITY_CURED);
}

class ParasiteParametersYAMLTest : public ::testing::Test {
protected:
    ParasiteParameters parasite_parameters;

    void SetUp() override {
        // Initialize ParasiteDensityLevels
        ParasiteParameters::ParasiteDensityLevels density_levels;
        density_levels.set_log_parasite_density_cured(1.1);
        density_levels.set_log_parasite_density_from_liver(2.2);
        density_levels.set_log_parasite_density_asymptomatic(3.3);
        density_levels.set_log_parasite_density_clinical(4.4);
        density_levels.set_log_parasite_density_clinical_from(5.5);
        density_levels.set_log_parasite_density_clinical_to(6.6);
        density_levels.set_log_parasite_density_detectable(7.7);
        density_levels.set_log_parasite_density_detectable_pfpr(8.8);
        density_levels.set_log_parasite_density_pyrogenic(9.9);

        // Initialize RecombinationParameters
        ParasiteParameters::RecombinationParameters recombination;
        recombination.set_within_chromosome_recombination_rate(0.1);
        recombination.set_using_free_recombination(true);

        // Set values in ParasiteParameters
        parasite_parameters.set_parasite_density_levels(density_levels);
        parasite_parameters.set_recombination_parameters(recombination);
    }
};

// Test encoding functionality for ParasiteParameters
TEST_F(ParasiteParametersYAMLTest, EncodeParasiteParameters) {
    YAML::Node node = YAML::convert<ParasiteParameters>::encode(parasite_parameters);

    // Validate ParasiteDensityLevels encoding
    EXPECT_EQ(node["parasite_density_levels"]["log_parasite_density_cured"].as<double>(), 1.1);
    EXPECT_EQ(node["parasite_density_levels"]["log_parasite_density_from_liver"].as<double>(), 2.2);
    EXPECT_EQ(node["parasite_density_levels"]["log_parasite_density_asymptomatic"].as<double>(), 3.3);
    EXPECT_EQ(node["parasite_density_levels"]["log_parasite_density_clinical"].as<double>(), 4.4);

    // Validate RecombinationParameters encoding
    EXPECT_EQ(node["recombination_parameters"]["within_chromosome_recombination_rate"].as<double>(), 0.1);
    EXPECT_EQ(node["recombination_parameters"]["using_free_recombination"].as<bool>(), true);
}

// Test decoding functionality for ParasiteParameters
TEST_F(ParasiteParametersYAMLTest, DecodeParasiteParameters) {
    YAML::Node node;
    node["parasite_density_levels"]["log_parasite_density_cured"] = 1.1;
    node["parasite_density_levels"]["log_parasite_density_from_liver"] = 2.2;
    node["parasite_density_levels"]["log_parasite_density_asymptomatic"] = 3.3;
    node["parasite_density_levels"]["log_parasite_density_clinical"] = 4.4;
    node["parasite_density_levels"]["log_parasite_density_clinical_from"] = 5.5;
    node["parasite_density_levels"]["log_parasite_density_clinical_to"] = 6.6;
    node["parasite_density_levels"]["log_parasite_density_detectable"] = 7.7;
    node["parasite_density_levels"]["log_parasite_density_detectable_pfpr"] = 8.8;
    node["parasite_density_levels"]["log_parasite_density_pyrogenic"] = 9.9;

    node["recombination_parameters"]["within_chromosome_recombination_rate"] = 0.1;
    node["recombination_parameters"]["using_free_recombination"] = true;

    ParasiteParameters decoded_parameters;
    EXPECT_NO_THROW(YAML::convert<ParasiteParameters>::decode(node, decoded_parameters));

    // Validate ParasiteDensityLevels decoding
    EXPECT_EQ(decoded_parameters.get_parasite_density_levels().get_log_parasite_density_cured(), 1.1);
    EXPECT_EQ(decoded_parameters.get_parasite_density_levels().get_log_parasite_density_from_liver(), 2.2);
    EXPECT_EQ(decoded_parameters.get_parasite_density_levels().get_log_parasite_density_asymptomatic(), 3.3);
    EXPECT_EQ(decoded_parameters.get_parasite_density_levels().get_log_parasite_density_clinical(), 4.4);

    // Validate RecombinationParameters decoding
    EXPECT_EQ(decoded_parameters.get_recombination_parameters().get_within_chromosome_recombination_rate(), 0.1);
    EXPECT_EQ(decoded_parameters.get_recombination_parameters().get_using_free_recombination(), true);
}

// Test for decoding with missing fields
TEST_F(ParasiteParametersYAMLTest, DecodeParasiteParametersMissingField) {
    YAML::Node node;
    node["parasite_density_levels"]["log_parasite_density_cured"] = 1.1;  // Missing other fields

    ParasiteParameters decoded_parameters;
    EXPECT_THROW(YAML::convert<ParasiteParameters>::decode(node, decoded_parameters), std::runtime_error);
}
