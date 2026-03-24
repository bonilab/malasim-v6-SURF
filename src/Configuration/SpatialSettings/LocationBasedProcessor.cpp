#include "LocationBasedProcessor.h"

#include <spdlog/spdlog.h>

#include "Simulation/Model.h"

void LocationBasedProcessor::process_config() {
  spdlog::info("Processing using LocationBasedProcessor");
  auto location_based = get_spatial_settings()->get_node().as<SpatialSettings::LocationBased>();
  // process location_based
  std::vector<Spatial::Location> location_db;
  location_db.clear();
  // loop index based on location_based_.get_location_info()
  for (const auto &location : location_based.locations) {
    auto new_location = location;
    auto pop_index =
        (location.id >= location_based.population_size_by_location.size()) ? 0 : location.id;
    new_location.population_size = location_based.population_size_by_location[pop_index];

    auto age_index =
        (location.id >= location_based.age_distribution_by_location.size()) ? 0 : location.id;
    new_location.age_distribution = location_based.age_distribution_by_location[age_index];

    auto beta_index = (location.id >= location_based.beta_by_location.size()) ? 0 : location.id;
    new_location.beta = location_based.beta_by_location[beta_index];

    auto p_treatment_under_5_index =
        (location.id >= location_based.p_treatment_under_5_by_location.size()) ? 0 : location.id;
    new_location.p_treatment_under_5 =
        location_based.p_treatment_under_5_by_location[p_treatment_under_5_index];

    auto p_treatment_over_5_index =
        (location.id >= location_based.p_treatment_over_5_by_location.size()) ? 0 : location.id;
    new_location.p_treatment_over_5 =
        location_based.p_treatment_over_5_by_location[p_treatment_over_5_index];

    location_db.push_back(new_location);
  }

  // populate the distance matrix
  auto number_of_location = location_db.size();

  auto spatial_distance_matrix =
      std::vector<std::vector<double>>(static_cast<uint64_t>(number_of_location));
  for (auto from_location = 0; from_location < number_of_location; from_location++) {
    spatial_distance_matrix[from_location].resize(static_cast<uint64_t>(number_of_location));
    for (auto to_location = 0; to_location < number_of_location; to_location++) {
      spatial_distance_matrix[from_location][to_location] =
          Spatial::Coordinate::calculate_distance_in_km(location_db[from_location].coordinate,
                                                        location_db[to_location].coordinate);
    }
  }

  // show error if not size equal to number of locations or 1
  if (location_based.age_distribution_by_location.size() != number_of_location
      && location_based.age_distribution_by_location.size() != 1) {
    throw std::runtime_error(
        "Age distribution by location size should be equal to number of locations or 1");
  }
  // show error if not size equal to initial age structure size
  for (const auto &age_distribution : location_based.age_distribution_by_location) {
    if (age_distribution.size() <= 0) {
      throw std::runtime_error("Number of age distribution by location should be greater than 0");
    }
  }

  // assign back to spatial settings
  get_spatial_settings()->set_location_db(location_db);
  get_spatial_settings()->set_spatial_distance_matrix(spatial_distance_matrix);
  get_spatial_settings()->set_number_of_locations(number_of_location);

  // ensure spatial_data and its admin_level_manager are prepared
  get_spatial_settings()->set_spatial_data(std::make_unique<SpatialData>(get_spatial_settings()));
  get_spatial_settings()->spatial_data()->initialize_admin_boundaries();
}

