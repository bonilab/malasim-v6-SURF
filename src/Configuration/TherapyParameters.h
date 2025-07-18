#ifndef THERAPYPARAMETERS_H
#define THERAPYPARAMETERS_H

#include <yaml-cpp/yaml.h>
#include <stdexcept>
#include <vector>
#include <map>
#include <spdlog/spdlog.h>
#include "IConfigData.h"
#include "Treatment/Therapies/SCTherapy.h"
#include "Treatment/Therapies/TherapyBuilder.h"
#include "Utils/Helpers/NumberHelpers.h"

class TherapyParameters: public IConfigData {
public:
  // Inner class: TherapyInfo
  class TherapyInfo {
  public:
    // Getters and Setters
    [[nodiscard]] const std::vector<int>& get_drug_ids() const { return drug_ids_; }
    void set_drug_ids(const std::vector<int>& value) { drug_ids_ = value; }

    [[nodiscard]] const std::vector<int>& get_dosing_days() const { return dosing_days_; }
    void set_dosing_days(const std::vector<int>& value) { dosing_days_ = value; }

    [[nodiscard]] const std::vector<int>& get_therapy_ids() const { return therapy_ids_; }
    void set_therapy_ids(const std::vector<int>& value) { therapy_ids_ = value; }

    [[nodiscard]] const std::vector<int>& get_regiment() const { return regiment_; }
    void set_regiment(const std::vector<int>& value) { regiment_ = value; }

    [[nodiscard]] const std::string& get_name() const { return name_; }
    void set_name(const std::string& value) { name_ = value; }

  private:
    std::string name_;
    std::vector<int> drug_ids_;
    std::vector<int> dosing_days_;
    std::vector<int> therapy_ids_;
    std::vector<int> regiment_;
  };

  // Getters and Setters for TherapyParameters
  [[nodiscard]] int get_tf_testing_day() const { return tf_testing_day_; }
  void set_tf_testing_day(const int value) { tf_testing_day_ = value; }

  [[nodiscard]] double get_tf_rate() const { return tf_rate_; }
  void set_tf_rate(const double value) { tf_rate_ = value; }

  [[nodiscard]] int get_recurrence_therapy_id() const { return recurrence_therapy_id_; }
  void set_recurrence_therapy_id(const int value) { recurrence_therapy_id_ = value; }

  [[nodiscard]] const std::map<int, TherapyInfo>& get_therapy_db_raw() const { return therapy_db_raw_; }
  void set_therapy_db_raw(const std::map<int, TherapyInfo>& value) { therapy_db_raw_ = value; }

  [[nodiscard]] YAML::Node get_node() const { return node_; }
  void set_node(const YAML::Node& value) { node_ = value; }

  std::unique_ptr<Therapy> read_therapy(const YAML::Node &n, const int &therapy_id) {
    const auto t_id = NumberHelpers::number_to_string<int>(therapy_id);
    auto therapy = TherapyBuilder::build(n[t_id], therapy_id);
    return therapy;
  }

  void process_config() override; 

private:
  int tf_testing_day_ = 28;
  double tf_rate_ = 0.1;
  int recurrence_therapy_id_ = -1;
  std::map<int, TherapyInfo> therapy_db_raw_;  // Changed from vector to map
  YAML::Node node_;
};

namespace YAML {

// TherapyParameters::TherapyInfo YAML conversion
template<>
struct convert<TherapyParameters::TherapyInfo> {
  static Node encode(const TherapyParameters::TherapyInfo& rhs) {
    Node node;
    if(!rhs.get_drug_ids().empty()) {
      node["drug_ids"] = rhs.get_drug_ids();
      node["dosing_days"] = rhs.get_dosing_days();
    }
    return node;
  }

  static bool decode(const Node& node, TherapyParameters::TherapyInfo& rhs) {
    if (!node["drug_ids"] && !node["dosing_days"] && !node["therapy_ids"] && !node["regiment"]) {
      throw std::runtime_error("Missing fields in TherapyParameters::TherapyInfo");
    }
    if(node["drug_ids"]) {
      rhs.set_drug_ids(node["drug_ids"].as<std::vector<int>>());
    }
    if(node["dosing_days"]) {
      rhs.set_dosing_days(node["dosing_days"].as<std::vector<int>>());
    }
    if(node["therapy_ids"]) {
      rhs.set_therapy_ids(node["therapy_ids"].as<std::vector<int>>());
    }
    if(node["regiment"]) {
      rhs.set_regiment(node["regiment"].as<std::vector<int>>());
    }
    if(node["name"]) {
      rhs.set_name(node["name"].as<std::string>());
    }
    return true;
  }
};

// TherapyParameters YAML conversion
template<>
struct convert<TherapyParameters> {
  static Node encode(const TherapyParameters& rhs) {
    Node node;
    node["tf_testing_day"] = rhs.get_tf_testing_day();
    node["tf_rate"] = rhs.get_tf_rate();

    node["recurrence_therapy_id"] = rhs.get_recurrence_therapy_id();
    // Encode therapy_db as a map
    Node therapy_db_node;
    for (const auto& [key, value] : rhs.get_therapy_db_raw()) {
      therapy_db_node[key] = value;
    }
    node["therapy_db"] = therapy_db_node;
    return node;
  }

  static bool decode(const Node& node, TherapyParameters& rhs) {
    if (!node["tf_testing_day"] || !node["tf_rate"] || !node["therapy_db"]) {
      throw std::runtime_error("Missing fields in TherapyParameters");
    }
    rhs.set_tf_testing_day(node["tf_testing_day"].as<int>());
    rhs.set_tf_rate(node["tf_rate"].as<double>());

    if (node["recurrence_therapy_id"]) {
      rhs.set_recurrence_therapy_id(node["recurrence_therapy_id"].as<int>());
    } else {
      rhs.set_recurrence_therapy_id(-1);
    }

    // Decode therapy_db as a map
    std::map<int, TherapyParameters::TherapyInfo> therapy_db;
    for (const auto& element : node["therapy_db"]) {
      int key = element.first.as<int>();
      therapy_db[key] = element.second.as<TherapyParameters::TherapyInfo>();
    }
    rhs.set_therapy_db_raw(therapy_db);
    rhs.set_node(node["therapy_db"]);
    return true;
  }
};

}  // namespace YAML

#endif //THERAPYPARAMETERS_H
