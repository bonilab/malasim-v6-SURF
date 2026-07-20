#ifndef IMMUNESYSTEMPARAMETEROVERRIDES_H
#define IMMUNESYSTEMPARAMETEROVERRIDES_H

#include <map>
#include <stdexcept>
#include <string>

#include <spdlog/spdlog.h>
#include <yaml-cpp/yaml.h>

// Full config-path keys used in override maps (colon-separated YAML path segments).
// These match the layout of the config file so users can identify what is being overridden.
namespace ImmuneSystemOverridePaths {
  constexpr const char* K_P_CI_SYMP =
      "epidemiological_parameters:allow_new_coinfection_to_cause_symptoms:probability";
  constexpr const char* K_Z =
      "immune_system_parameters:immune_effect_on_progression_to_clinical";
  constexpr const char* K_KAPPA =
      "immune_system_parameters:factor_effect_age_mature_immunity";
  constexpr const char* K_MIDPOINT =
      "immune_system_parameters:midpoint";
  constexpr const char* K_P_SEEK_BASE =
      "epidemiological_parameters:age_based_probability_of_seeking_treatment:power:base";
  constexpr const char* K_MUTATION_PROB =
      "genotype_parameters:mutation_probability_per_locus";
  constexpr const char* K_DEFAULT_CNV_REVERSION_MULTIPLIER =
      "genotype_parameters:default_cnv_reversion_multiplier";
}  // namespace ImmuneSystemOverridePaths

// A single override set: a map from full config-path key to the override value.
struct ImmuneSystemParametercalibration_id {
  std::map<std::string, double> overrides;

  [[nodiscard]] bool has(const std::string &path) const {
    return overrides.count(path) > 0;
  }

  [[nodiscard]] double get(const std::string &path, double default_val) const {
    const auto it = overrides.find(path);
    return (it != overrides.end()) ? it->second : default_val;
  }
};

// Container for all override calibration_ids and the index selecting which is used in simulation.
class ImmuneSystemParameterOverrides {
public:
  [[nodiscard]] int get_chosen_calibration_id() const { return chosen_calibration_id_; }
  void set_chosen_calibration_id(const int value) { chosen_calibration_id_ = value; }

  [[nodiscard]] bool get_random_selection() const { return random_selection_; }
  void set_random_selection(const bool value) { random_selection_ = value; }

  [[nodiscard]] const std::map<int, ImmuneSystemParametercalibration_id> &get_calibration_ids() const {
    return calibration_ids_;
  }
  void set_calibration_ids(const std::map<int, ImmuneSystemParametercalibration_id> &value) {
    calibration_ids_ = value;
  }

  [[nodiscard]] bool has_selected_calibration_id() const {
    return calibration_ids_.count(chosen_calibration_id_) > 0;
  }

  [[nodiscard]] const ImmuneSystemParametercalibration_id &get_selected_calibration_id() const {
    auto it = calibration_ids_.find(chosen_calibration_id_);
    if (it == calibration_ids_.end()) {
      throw std::runtime_error(
          fmt::format("version6_pfpr_incidence_calibrations: chosen_calibration_id={} not found in calibration_ids",
                      chosen_calibration_id_));
    }
    return it->second;
  }

  void log_all() const {
    spdlog::info("ImmuneSystemParameterOverrides: chosen_calibration_id={}, random_selection={}",
                 chosen_calibration_id_, random_selection_);
    for (const auto &[idx, c] : calibration_ids_) {
      for (const auto &[path, val] : c.overrides) {
        spdlog::info("  calibration_id[{}]: {}={}", idx, path, val);
      }
    }
  }

private:
  int chosen_calibration_id_ = 0;
  bool random_selection_ = false;
  std::map<int, ImmuneSystemParametercalibration_id> calibration_ids_;
};

// YAML conversion for ImmuneSystemParameterOverrides
template <>
struct YAML::convert<ImmuneSystemParameterOverrides> {
  static Node encode(const ImmuneSystemParameterOverrides &rhs) {
    Node node;
    node["chosen_calibration_id"] = rhs.get_chosen_calibration_id();
    node["random_selection"] = rhs.get_random_selection();
    Node calibration_ids_node;
    for (const auto &[idx, c] : rhs.get_calibration_ids()) {
      Node cnode;
      for (const auto &[path, val] : c.overrides) {
        cnode[path] = val;
      }
      calibration_ids_node[idx] = cnode;
    }
    node["calibration_ids"] = calibration_ids_node;
    return node;
  }

  static bool decode(const Node &node, ImmuneSystemParameterOverrides &rhs) {
    if (!node["chosen_calibration_id"]) {
      throw std::runtime_error(
          "version6_pfpr_incidence_calibrations: missing 'chosen_calibration_id'");
    }
    if (!node["calibration_ids"]) {
      throw std::runtime_error(
          "version6_pfpr_incidence_calibrations: missing 'calibration_ids'");
    }

    rhs.set_chosen_calibration_id(node["chosen_calibration_id"].as<int>());
    if (node["random_selection"]) {
      rhs.set_random_selection(node["random_selection"].as<bool>());
    }

    std::map<int, ImmuneSystemParametercalibration_id> calibration_ids;
    YAML::Node cmap = node["calibration_ids"];
    for (auto it = cmap.begin(); it != cmap.end(); ++it) {
      const int idx = it->first.as<int>();
      YAML::Node cnode = it->second;

      ImmuneSystemParametercalibration_id c;
      for (auto kv = cnode.begin(); kv != cnode.end(); ++kv) {
        c.overrides[kv->first.as<std::string>()] = kv->second.as<double>();
      }
      calibration_ids[idx] = c;
    }
    rhs.set_calibration_ids(calibration_ids);
    return true;
  }
};

#endif  // IMMUNESYSTEMPARAMETEROVERRIDES_H