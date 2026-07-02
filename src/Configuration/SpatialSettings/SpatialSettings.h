#ifndef SPATIALSETTINGS_H
#define SPATIALSETTINGS_H
#include <date/date.h>
#include <spdlog/spdlog.h>
#include <yaml-cpp/yaml.h>

#include <cstddef>
#include <stdexcept>
#include <string>
#include <vector>

#include "Configuration/IConfigData.h"
#include "Spatial/GIS/SpatialData.h"
#include "Spatial/Location/Location.h"

// Class for SpatialSettings
class SpatialSettings : public IConfigData {
public:
  // Class for GridBased settings
  struct GridBased {
    struct AdministativeBoundaries {
      std::string name;
      std::string raster;
    };
    std::vector<Spatial::Location> locations;
    std::string population_raster;
    std::vector<AdministativeBoundaries> administrative_boundaries;
    std::string p_treatment_under_5_raster;
    std::string p_treatment_over_5_raster;
    std::string beta_raster;
    std::string ecoclimatic_raster;
    double cell_size{0.0};
    std::vector<std::vector<double>> age_distribution_by_location;
    int number_of_location{0};
  };

  // Class for LocationBased settings
  struct LocationBased {
    std::vector<Spatial::Location> locations;
    std::vector<std::vector<double>> age_distribution_by_location;
    std::vector<double> p_treatment_under_5_by_location;
    std::vector<double> p_treatment_over_5_by_location;
    std::vector<double> beta_by_location;
    std::vector<int> population_size_by_location;
    int number_of_locations{0};
  };

  // Getters and Setters for mode
  [[nodiscard]] const std::string &get_mode() const { return mode_; }
  void set_mode(const std::string &value) { mode_ = value; }

  [[nodiscard]] std::vector<std::vector<double>> &get_spatial_distance_matrix() {
    return spatial_distance_matrix_;
  }
  void set_spatial_distance_matrix(const std::vector<std::vector<double>> &value) {
    spatial_distance_matrix_ = value;
  }

  [[nodiscard]] size_t get_number_of_locations() const { return number_of_location_; }
  void set_number_of_locations(const size_t value) { number_of_location_ = value; }

  void set_node(const YAML::Node &value) { node_ = value; }
  [[nodiscard]] const YAML::Node &get_node() const { return node_; }

  [[nodiscard]] std::vector<Spatial::Location> &location_db() { return location_db_; }
  void set_location_db(const std::vector<Spatial::Location> &value) { location_db_ = value; }

  [[nodiscard]] SpatialData* spatial_data() { return spatial_data_.get(); }
  void set_spatial_data(std::unique_ptr<SpatialData> value) { spatial_data_ = std::move(value); }

  void process_config() override;

  void cross_validate();

  static constexpr std::string GRID_BASED_MODE = "grid_based";
  static constexpr std::string LOCATION_BASED_MODE = "location_based";

private:
  std::string mode_;  // "grid_based" or "location_based"
  // we have to store the yaml node as the process_config method will used it
  // and combine with other data in the model to populate the right data
  YAML::Node node_;

  std::vector<std::vector<double>> spatial_distance_matrix_;
  size_t number_of_location_{0};
  std::vector<Spatial::Location> location_db_;
  std::unique_ptr<SpatialData> spatial_data_{nullptr};
};

namespace YAML {
// Conversion for GridBased class
template <>
struct convert<SpatialSettings::GridBased> {
  static Node encode(const SpatialSettings::GridBased &rhs) {
    Node node;
    node["population_raster"] = rhs.population_raster;
    node["p_treatment_under_5_raster"] = rhs.p_treatment_under_5_raster;
    node["p_treatment_over_5_raster"] = rhs.p_treatment_over_5_raster;
    node["beta_raster"] = rhs.beta_raster;
    node["ecoclimatic_raster"] = rhs.ecoclimatic_raster;
    node["cell_size"] = rhs.cell_size;
    node["age_distribution_by_location"] = rhs.age_distribution_by_location;
    for (int i = 0; i < rhs.administrative_boundaries.size(); i++) {
      node["administrative_boundaries"][i]["name"] = rhs.administrative_boundaries[i].name;
      node["administrative_boundaries"][i]["raster"] = rhs.administrative_boundaries[i].raster;
    }
    return node;
  }

  static bool decode(const Node &node, SpatialSettings::GridBased &rhs) {
    // make error message more specific
    if (!node["population_raster"]) {
      throw std::runtime_error("Missing 'population_raster' field in grid-based settings.");
    }
    if (!node["administrative_boundaries"]) {
      throw std::runtime_error("Missing 'administrative_boundaries' field in grid-based settings.");
    }
    if (!node["p_treatment_under_5_raster"]) {
      throw std::runtime_error(
          "Missing 'p_treatment_under_5_raster' field in grid-based settings.");
    }
    if (!node["p_treatment_over_5_raster"]) {
      throw std::runtime_error("Missing 'p_treatment_over_5_raster' field in grid-based settings.");
    }
    if (!node["beta_raster"]) {
      throw std::runtime_error("Missing 'beta_raster' field in grid-based settings.");
    }
    if (!node["cell_size"]) {
      throw std::runtime_error("Missing 'cell_size' field in grid-based settings.");
    }
    if (!node["age_distribution_by_location"]) {
      throw std::runtime_error(
          "Missing 'age_distribution_by_location' field in grid-based settings.");
    }

    rhs.population_raster = node["population_raster"].as<std::string>();
    rhs.p_treatment_under_5_raster = node["p_treatment_under_5_raster"].as<std::string>();
    rhs.p_treatment_over_5_raster = node["p_treatment_over_5_raster"].as<std::string>();
    rhs.beta_raster = node["beta_raster"].as<std::string>();
    rhs.cell_size = node["cell_size"].as<double>();
    if (node["ecoclimatic_raster"]) {
      rhs.ecoclimatic_raster = node["ecoclimatic_raster"].as<std::string>();
    }
    std::vector<SpatialSettings::GridBased::AdministativeBoundaries> admin_boundaries;
    for (auto i = 0; i < node["administrative_boundaries"].size(); i++) {
      SpatialSettings::GridBased::AdministativeBoundaries admin_boundary;
      admin_boundary.name = node["administrative_boundaries"][i]["name"].as<std::string>();
      admin_boundary.raster = node["administrative_boundaries"][i]["raster"].as<std::string>();
      admin_boundaries.push_back(admin_boundary);
    }
    rhs.administrative_boundaries = admin_boundaries;
    /* use one age distribution for all locations */
    rhs.age_distribution_by_location =
        node["age_distribution_by_location"].as<std::vector<std::vector<double>>>();
    return true;
  }
};

// Conversion for Location class
template <>
struct convert<Spatial::Location> {
  static Node encode(const Spatial::Location &rhs) {
    Node node;
    node.push_back(rhs.id);
    node.push_back(rhs.coordinate.latitude);
    node.push_back(rhs.coordinate.longitude);
    return node;
  }

  static bool decode(const Node &node, Spatial::Location &rhs) {
    if (!node.IsSequence() || node.size() != 3) {
      throw std::runtime_error("Invalid location info.");
    }
    rhs.id = node[0].as<int>();
    rhs.coordinate = Spatial::Coordinate(node[1].as<float>(), node[2].as<float>());
    return true;
  }
};

// Conversion for LocationBased class
template <>
struct convert<SpatialSettings::LocationBased> {
  static Node encode(const SpatialSettings::LocationBased &rhs) {
    Node node;
    node["location_info"] = rhs.locations;
    node["age_distribution_by_location"] = rhs.age_distribution_by_location;
    node["p_treatment_under_5_by_location"] = rhs.p_treatment_under_5_by_location;
    node["p_treatment_over_5_by_location"] = rhs.p_treatment_over_5_by_location;
    node["beta_by_location"] = rhs.beta_by_location;
    node["population_size_by_location"] = rhs.population_size_by_location;
    return node;
  }

  static bool decode(const Node &node, SpatialSettings::LocationBased &rhs) {
    // make error message more specific
    if (!node["location_info"]) {
      throw std::runtime_error("Missing 'location_info' field in location-based settings.");
    }
    if (!node["age_distribution_by_location"]) {
      throw std::runtime_error(
          "Missing 'age_distribution_by_location' field in location-based settings.");
    }
    if (!node["p_treatment_under_5_by_location"]) {
      throw std::runtime_error(
          "Missing 'p_treatment_under_5_by_location' field in location-based settings.");
    }
    if (!node["p_treatment_over_5_by_location"]) {
      throw std::runtime_error(
          "Missing 'p_treatment_over_5_by_location' field in location-based settings.");
    }
    if (!node["beta_by_location"]) {
      throw std::runtime_error("Missing 'beta_by_location' field in location-based settings.");
    }
    if (!node["population_size_by_location"]) {
      throw std::runtime_error(
          "Missing 'population_size_by_location' field in location-based settings.");
    }

    rhs.locations = std::move(node["location_info"].as<std::vector<Spatial::Location>>());
    rhs.age_distribution_by_location =
        std::move(node["age_distribution_by_location"].as<std::vector<std::vector<double>>>());
    rhs.p_treatment_under_5_by_location =
        std::move(node["p_treatment_under_5_by_location"].as<std::vector<double>>());
    rhs.p_treatment_over_5_by_location =
        std::move(node["p_treatment_over_5_by_location"].as<std::vector<double>>());
    rhs.beta_by_location = std::move(node["beta_by_location"].as<std::vector<double>>());
    rhs.population_size_by_location =
        std::move(node["population_size_by_location"].as<std::vector<int>>());

    spdlog::info("Location based settings decoded successfully with {} locations", rhs.locations.size());
    spdlog::info("Total locations: {}", rhs.population_size_by_location.size());
    for (int loc = 0; loc < rhs.locations.size(); loc++) {
      spdlog::info("\tLocation {}: id={}, lat={}, lon={}, pop_size={}", loc, rhs.locations[loc].id,
                   rhs.locations[loc].coordinate.latitude, rhs.locations[loc].coordinate.longitude,
                   rhs.population_size_by_location[loc]);
    }
    return true;
  }
};

// Conversion for SpatialSettings class
template <>
struct convert<SpatialSettings> {
  static Node encode(const SpatialSettings &rhs) {
    Node node;
    node["mode"] = rhs.get_mode();
    node[rhs.get_mode()] = rhs.get_node();
    return node;
  }

  static bool decode(const Node &node, SpatialSettings &rhs) {
    if (!node["mode"]) { throw std::runtime_error("Missing 'mode' field."); }

    auto mode = node["mode"].as<std::string>();
    rhs.set_mode(mode);

    /* Here we only get the node for the mode specified in the configuration file
     * the specific processor will convert the node to the appropriate class
     * and process the configuration
     */
    if (mode != SpatialSettings::GRID_BASED_MODE && mode != SpatialSettings::LOCATION_BASED_MODE) {
      throw std::runtime_error("Unknown mode in 'spatial_settings'.");
    }

    auto node_name = mode == SpatialSettings::GRID_BASED_MODE
                         ? SpatialSettings::GRID_BASED_MODE
                         : SpatialSettings::LOCATION_BASED_MODE;
    if (!node[node_name]) { throw std::runtime_error("Missing " + node_name + " settings."); }
    rhs.set_node(node[node_name]);
    return true;
  }
};
}  // namespace YAML

#endif  // SPATIALSETTINGS_H
