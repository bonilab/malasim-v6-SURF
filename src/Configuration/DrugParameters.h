#ifndef DRUGPARAMETERS_H
#define DRUGPARAMETERS_H
#include <yaml-cpp/yaml.h>

#include <map>
#include <stdexcept>
#include <string>
#include <vector>

#include "IConfigData.h"

class DrugParameters : public IConfigData {
public:
  // Inner class: DrugInfo
  class DrugInfo {
  public:
    // Getters and Setters
    [[nodiscard]] const std::string &get_name() const { return name_; }
    void set_name(const std::string &value) { name_ = value; }

    [[nodiscard]] double get_half_life() const { return half_life_; }
    void set_half_life(const double value) { half_life_ = value; }

    [[nodiscard]] double get_maximum_parasite_killing_rate() const {
      return maximum_parasite_killing_rate_;
    }
    void set_maximum_parasite_killing_rate(const double value) {
      maximum_parasite_killing_rate_ = value;
    }

    [[nodiscard]] int get_n() const { return n_; }
    void set_n(const int value) { n_ = value; }

    [[nodiscard]] const std::vector<double> &get_age_specific_drug_concentration_sd() const {
      return age_specific_drug_concentration_sd_;
    }
    void set_age_specific_drug_concentration_sd(const std::vector<double> &value) {
      age_specific_drug_concentration_sd_ = value;
    }

    [[nodiscard]] const std::vector<double> &get_age_specific_drug_absorption() const {
      return age_specific_drug_absorption_;
    }
    void set_age_specific_drug_absorption(const std::vector<double> &value) {
      age_specific_drug_absorption_ = value;
    }

    [[nodiscard]] int get_k() const { return k_; }
    void set_k(const int value) { k_ = value; }

    [[nodiscard]] double get_base_ec50() const { return base_EC50_; }
    void set_base_EC50(const double value) { base_EC50_ = value; }

  private:
    std::string name_;
    double half_life_ = -1;
    double maximum_parasite_killing_rate_ = -1;
    int n_ = -1;
    std::vector<double> age_specific_drug_concentration_sd_;
    std::vector<double> age_specific_drug_absorption_;
    int k_ = -1;
    double base_EC50_ = -1;
  };

  // Getters and Setters for DrugParameters
  [[nodiscard]] const std::map<int, DrugInfo> &get_drug_db_raw() const { return drug_infos_; }
  void set_drug_db_raw(const std::map<int, DrugInfo> &value) { drug_infos_ = value; }

  // process config data
  void process_config() override;

private:
  std::map<int, DrugInfo> drug_infos_;
};

namespace YAML {

// DrugParameters::DrugInfo YAML conversion
template <>
struct convert<DrugParameters::DrugInfo> {
  static Node encode(const DrugParameters::DrugInfo &rhs) {
    Node node;
    node["name"] = rhs.get_name();
    node["half_life"] = rhs.get_half_life();
    node["maximum_parasite_killing_rate"] = rhs.get_maximum_parasite_killing_rate();
    node["n"] = rhs.get_n();
    node["age_specific_drug_concentration_sd"] = rhs.get_age_specific_drug_concentration_sd();
    node["age_specific_drug_absorption"] = rhs.get_age_specific_drug_absorption();
    node["k"] = rhs.get_k();
    node["base_EC50"] = rhs.get_base_ec50();
    return node;
  }

  static bool decode(const Node &node, DrugParameters::DrugInfo &rhs) {
    if (!node["name"] || !node["half_life"] || !node["maximum_parasite_killing_rate"] || !node["n"]
        || !node["age_specific_drug_concentration_sd"] || !node["k"] || !node["base_EC50"]) {
      throw std::runtime_error("Missing fields in DrugParameters::DrugInfo");
    }
    rhs.set_name(node["name"].as<std::string>());
    rhs.set_half_life(node["half_life"].as<double>());
    rhs.set_maximum_parasite_killing_rate(node["maximum_parasite_killing_rate"].as<double>());
    rhs.set_n(node["n"].as<int>());
    rhs.set_age_specific_drug_concentration_sd(
        node["age_specific_drug_concentration_sd"].as<std::vector<double>>());
    if (node["age_specific_drug_absorption"])
      rhs.set_age_specific_drug_absorption(
          node["age_specific_drug_absorption"].as<std::vector<double>>());
    rhs.set_k(node["k"].as<int>());
    rhs.set_base_EC50(node["base_EC50"].as<double>());
    return true;
  }
};

// DrugParameters YAML conversion
template <>
struct convert<DrugParameters> {
  static Node encode(const DrugParameters &rhs) {
    Node node;
    for (const auto &[key, value] : rhs.get_drug_db_raw()) { node["drug_db"][key] = value; }
    return node;
  }

  static bool decode(const Node &node, DrugParameters &rhs) {
    if (!node["drug_db"]) { throw std::runtime_error("Missing 'drug_db' field in DrugParameters"); }
    std::map<int, DrugParameters::DrugInfo> drug_db;
    for (const auto &element : node["drug_db"]) {
      int key = element.first.as<int>();
      drug_db[key] = element.second.as<DrugParameters::DrugInfo>();
    }
    rhs.set_drug_db_raw(drug_db);
    return true;
  }
};

}  // namespace YAML

#endif  // DRUGPARAMETERS_H
