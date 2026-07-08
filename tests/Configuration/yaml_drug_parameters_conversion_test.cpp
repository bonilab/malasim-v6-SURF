#include <gtest/gtest.h>
#include <yaml-cpp/yaml.h>

#include "Configuration/DrugParameters.h"

class DrugParametersTest : public ::testing::Test {
protected:
  DrugParameters drug_parameters;
  std::vector<double> age_specific_drug_absorption_vec;

  void SetUp() override {
    age_specific_drug_absorption_vec = std::vector<double>(15, 0.0);
    for (int i = 0; i < 15; i++) { age_specific_drug_absorption_vec[i] = i * 0.1; }
    DrugParameters::DrugInfo drug_info_art;
    drug_info_art.set_name("Artemisinin");
    drug_info_art.set_half_life(0.0);
    drug_info_art.set_maximum_parasite_killing_rate(0.999);
    drug_info_art.set_n(25);
    drug_info_art.set_age_specific_drug_concentration_sd(std::vector<double>(15, 0.4));
    drug_info_art.set_k(4);
    drug_info_art.set_base_EC50(0.75);

    // Insert the drug_info into the drug_parameters map
    std::map<int, DrugParameters::DrugInfo> drug_db;
    drug_db[0] = drug_info_art;
    drug_parameters.set_drug_db_raw(drug_db);
  }
};

TEST_F(DrugParametersTest, EncodeDrugParameters) {
  YAML::Node node = YAML::convert<DrugParameters>::encode(drug_parameters);

  // Manually add lumefantrine to the map in the node
  node["drug_db"]["1"]["name"] = "Lumefantrine";
  node["drug_db"]["1"]["half_life"] = 4.5;
  node["drug_db"]["1"]["maximum_parasite_killing_rate"] = 0.99;
  node["drug_db"]["1"]["n"] = 20;
  node["drug_db"]["1"]["age_specific_drug_concentration_sd"] = std::vector<double>(15, 0.4);
  node["drug_db"]["1"]["age_specific_drug_absorption"] = age_specific_drug_absorption_vec;
  node["drug_db"]["1"]["k"] = 4;
  node["drug_db"]["1"]["base_EC50"] = 0.6;

  EXPECT_EQ(node["drug_db"]["1"]["name"].as<std::string>(), "Lumefantrine");
  EXPECT_EQ(node["drug_db"]["1"]["half_life"].as<double>(), 4.5);
  EXPECT_EQ(node["drug_db"]["1"]["maximum_parasite_killing_rate"].as<double>(), 0.99);
  EXPECT_EQ(node["drug_db"]["1"]["n"].as<int>(), 20);
  EXPECT_EQ(node["drug_db"]["1"]["age_specific_drug_concentration_sd"].as<std::vector<double>>(),
            std::vector<double>(15, 0.4));
  EXPECT_EQ(node["drug_db"]["1"]["age_specific_drug_absorption"].as<std::vector<double>>(),
            age_specific_drug_absorption_vec);
  EXPECT_EQ(node["drug_db"]["1"]["k"].as<int>(), 4);
  EXPECT_EQ(node["drug_db"]["1"]["base_EC50"].as<double>(), 0.6);
}

TEST_F(DrugParametersTest, DecodeDrugParameters) {
  YAML::Node node;
  node["drug_db"]["0"]["name"] = "Artemisinin";
  node["drug_db"]["0"]["half_life"] = 0.0;
  node["drug_db"]["0"]["maximum_parasite_killing_rate"] = 0.999;
  node["drug_db"]["0"]["n"] = 25;
  node["drug_db"]["0"]["age_specific_drug_concentration_sd"] = std::vector<double>(15, 0.4);
  node["drug_db"]["0"]["k"] = 4;
  node["drug_db"]["0"]["base_EC50"] = 0.75;

  DrugParameters decoded_parameters;
  EXPECT_NO_THROW(YAML::convert<DrugParameters>::decode(node, decoded_parameters));

  EXPECT_EQ(decoded_parameters.get_drug_db_raw().at(0).get_name(), "Artemisinin");
  EXPECT_EQ(decoded_parameters.get_drug_db_raw().at(0).get_half_life(), 0.0);
  EXPECT_EQ(decoded_parameters.get_drug_db_raw().at(0).get_maximum_parasite_killing_rate(), 0.999);
  EXPECT_EQ(decoded_parameters.get_drug_db_raw().at(0).get_n(), 25);
  EXPECT_EQ(decoded_parameters.get_drug_db_raw().at(0).get_age_specific_drug_concentration_sd(),
            std::vector<double>(15, 0.4));
  EXPECT_EQ(decoded_parameters.get_drug_db_raw().at(0).get_k(), 4);
  EXPECT_EQ(decoded_parameters.get_drug_db_raw().at(0).get_base_ec50(), 0.75);
}

// Test for decoding with missing fields
TEST_F(DrugParametersTest, DecodeDrugParametersMissingField) {
  // This node has missing required fields such as `n` and `base_EC50`
  YAML::Node node;
  node["drug_db"]["0"]["name"] = "Artemisinin";
  node["drug_db"]["0"]["half_life"] = 0.0;
  node["drug_db"]["0"]["maximum_parasite_killing_rate"] = 0.999;

  DrugParameters decoded_parameters;
  // Expect an exception due to missing required fields
  EXPECT_THROW(YAML::convert<DrugParameters>::decode(node, decoded_parameters), std::runtime_error);
}
