#include "SpatialSettings.h"

#include <spdlog/spdlog.h>

#include "GridBasedProcessor.h"
#include "LocationBasedProcessor.h"
#include "Simulation/Model.h"

void SpatialSettings::process_config() {
  spdlog::info("Processing SpatialSettings");
  std::unique_ptr<ISpatialSettingsProcessor> processor = nullptr;
  if (mode_ == GRID_BASED_MODE) {
    processor = std::make_unique<GridBasedProcessor>(this);
  } else if (mode_ == LOCATION_BASED_MODE) {
    processor = std::make_unique<LocationBasedProcessor>(this);
  } else {
    throw std::runtime_error("Unknown mode in 'spatial_settings'.");
  }
  processor->process_config();
}

void SpatialSettings::cross_validate() {
  // Check if mode is either grid_based or location_based
  if (mode_ != "grid_based" && mode_ != "location_based") {
    throw std::invalid_argument("Spatial mode should be either grid_based or location_based");
  }
  // If mode is grid_based, check if all raster file paths are provided
  if (mode_ == "grid_based") {
    const auto grid_based = node_.as<SpatialSettings::GridBased>();
    if (grid_based.population_raster.empty() || grid_based.beta_raster.empty()
        || grid_based.p_treatment_over_5_raster.empty()
        || grid_based.p_treatment_under_5_raster.empty()) {
      throw std::invalid_argument(
          "All raster file paths should be provided for grid based spatial mode");
    }
    // Check if age_distribution_by_location size is different from 1
    if (grid_based.age_distribution_by_location.size() != 1) {
      throw std::invalid_argument(
          "Age distribution using raster must be 1 location (to distribute equally)");
    }
    // Check if age_distribution_by_location size matched initial_age_structure size
    if (grid_based.age_distribution_by_location[0].size()
        != Model::get_config()->get_population_demographic().get_initial_age_structure().size()) {
      throw std::invalid_argument(
          "Age distribution by raster must have size match initial age structure size");
    }
  }
  // Location based
  if (mode_ == "location_based") {
    const auto location_based = node_.as<SpatialSettings::LocationBased>();
    if (location_based.population_size_by_location.empty() || location_based.locations.empty()
        || location_based.age_distribution_by_location.empty()
        || location_based.p_treatment_under_5_by_location.empty()
        || location_based.p_treatment_over_5_by_location.empty()
        || location_based.beta_by_location.empty()) {
      throw std::invalid_argument(
          "All locations should be provided for location based spatial mode");
    }
    if (location_based.age_distribution_by_location.size() != location_based.locations.size()
        && location_based.age_distribution_by_location.size() != 1) {
      throw std::invalid_argument(
          "Age distribution by location size should be equal to number of locations or 1");
    }

    // Check if age_distribution_by_location size matched initial_age_structure size
    for (const auto &age_dist : location_based.age_distribution_by_location) {
      if (age_dist.size()
          != Model::get_config()->get_population_demographic().get_initial_age_structure().size()) {
        spdlog::info("Age distribution by location size: {}", age_dist.size());
        spdlog::info(
            "Initial age structure size: {}",
            Model::get_config()->get_population_demographic().get_initial_age_structure().size());
        throw std::invalid_argument(
            "Age distribution by location size should match initial age structure size");
      }
    }
  }
}

