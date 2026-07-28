#ifndef RAPTSETTINGS_H
#define RAPTSETTINGS_H
#include <spdlog/spdlog.h>
#include <yaml-cpp/yaml.h>

#include "IConfigData.h"
#include "date/date.h"

class RaptSettings : public IConfigData {
public:
  [[nodiscard]] bool get_is_enabled() const { return is_enabled_; }
  void set_is_enabled(const bool value) { is_enabled_ = value; }
  [[nodiscard]] int get_period() const { return period_; }
  void set_period(const int value) { period_ = value; }
  [[nodiscard]] int get_therapy_id() const { return therapy_id_; }
  void set_therapy_id(const int value) { therapy_id_ = value; }
  [[nodiscard]] double get_compliance() const { return compliance_; }
  void set_compliance(const double value) { compliance_ = value; }
  [[nodiscard]] int get_age_start() const { return age_start_; }
  void set_age_start(const int value) { age_start_ = value; }
  [[nodiscard]] date::year_month_day get_start_date() const { return start_date_; }
  [[nodiscard]] int get_start_date_as_day() const { return start_date_as_day_; }
  void set_start_date(const date::year_month_day value) { start_date_ = value; }

private:
  bool is_enabled_ = true;
  int period_ = 12;
  int therapy_id_ = 1;
  double compliance_ = 0.7;
  int age_start_ = 18;
  int start_date_as_day_ = 0;
  date::year_month_day start_date_ =
      date::year_month_day(date::year(2000), date::month(1), date::day(1));

public:
  void process_config() override {}

  void process_config_with_starting_date(const date::year_month_day starting_date) {
    spdlog::info("Processing RaptSettings");
    start_date_as_day_ = (date::sys_days{start_date_} - date::sys_days(starting_date)).count();
  }
};

// Specialization of convert for the SimulationTimeframe class
template <>
struct YAML::convert<RaptSettings> {
  static Node encode(const RaptSettings &rhs) {
    Node node;
    node["enabled"] = rhs.get_is_enabled();
    node["period"] = rhs.get_period();
    node["therapy_id"] = rhs.get_therapy_id();
    node["compliance"] = rhs.get_compliance();
    node["age_start"] = rhs.get_age_start();
    node["start_date"] = format("%Y/%m/%d", rhs.get_start_date());
    return node;
  }

  static bool decode(const Node &node, RaptSettings &rhs) {
    if (!node["enabled"]) { throw std::runtime_error("Missing 'enabled' field in RaptSettings."); }
    rhs.set_is_enabled(node["enabled"].as<bool>());
    if (rhs.get_is_enabled()) {
      if (!node["period"]) { throw std::runtime_error("Missing 'period' field in RaptSettings."); }
      if (!node["therapy_id"]) {
        throw std::runtime_error("Missing 'therapy_id' field in RaptSettings.");
      }
      if (!node["compliance"]) {
        throw std::runtime_error("Missing 'compliance' field in RaptSettings.");
      }
      if (!node["age_start"]) {
        throw std::runtime_error("Missing 'age_start' field in RaptSettings.");
      }
      if (!node["start_date"]) {
        throw std::runtime_error("Missing 'start_date' field in RaptSettings.");
      }
      rhs.set_period(node["period"].as<int>());
      rhs.set_therapy_id(node["therapy_id"].as<int>());
      rhs.set_compliance(node["compliance"].as<double>());
      rhs.set_age_start(node["age_start"].as<int>());
      rhs.set_start_date(node["start_date"].as<date::year_month_day>());
    }
    return true;
  }
};
#endif  // RAPTSETTINGS_H
