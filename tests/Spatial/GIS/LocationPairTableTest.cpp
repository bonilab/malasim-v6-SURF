#include <gtest/gtest.h>

#include <cmath>
#include <vector>

#include "Spatial/GIS/GridPairTable.h"
#include "Spatial/Location/Location.h"
#include "Spatial/Movement/MovementKernel.h"

TEST(GridPairTableTest, LookupMatchesDenseEuclideanDistances) {
  std::vector<Spatial::Location> locations(4);
  locations[0].coordinate = {.latitude = 0, .longitude = 0};
  locations[1].coordinate = {.latitude = 0, .longitude = 2};
  locations[2].coordinate = {.latitude = 3, .longitude = 0};
  locations[3].coordinate = {.latitude = 3, .longitude = 2};

  constexpr float cell_size = 5.0F;
  const auto table = GridPairTable::make_distances(locations, cell_size);

  ASSERT_EQ(table.size(), locations.size());
  for (size_t from = 0; from < locations.size(); ++from) {
    for (size_t to = 0; to < locations.size(); ++to) {
      const double expected = std::sqrt(
          std::pow(
              cell_size * (locations[from].coordinate.latitude - locations[to].coordinate.latitude),
              2)
          + std::pow(
              cell_size
                  * (locations[from].coordinate.longitude - locations[to].coordinate.longitude),
              2));
      EXPECT_DOUBLE_EQ(table.at(static_cast<int>(from), static_cast<int>(to)), expected);
    }
  }
}

TEST(GridPairTableTest, DerivedTableKeepsZeroDistanceSentinel) {
  std::vector<Spatial::Location> locations(2);
  locations[0].coordinate = {.latitude = 0, .longitude = 0};
  locations[1].coordinate = {.latitude = 0, .longitude = 2};
  const auto distances = GridPairTable::make_distances(locations, 1.0F);
  const auto squared = Spatial::make_grid_movement_kernel(
      distances, [](double distance) { return distance * distance; }, "test squared distances");

  EXPECT_DOUBLE_EQ(squared.at(0, 0), 0.0);
  EXPECT_DOUBLE_EQ(squared.at(0, 1), 4.0);
  EXPECT_DOUBLE_EQ(squared.at(1, 0), 4.0);
  EXPECT_DOUBLE_EQ(squared.at(1, 1), 0.0);
}

TEST(GridPairTableTest, ReadsCoordinatesFromLocationDatabase) {
  std::vector<Spatial::Location> locations(2);
  locations[0].coordinate = {.latitude = 0, .longitude = 0};
  locations[1].coordinate = {.latitude = 0, .longitude = 2};
  const auto table = GridPairTable::make_distances(locations, 1.0F);

  EXPECT_DOUBLE_EQ(table.at(0, 1), 2.0);
  locations[1].coordinate.longitude = 1;
  EXPECT_DOUBLE_EQ(table.at(0, 1), 1.0);
}

TEST(GridPairTableTest, StorageIsSmallerThanDenseStorage) {
  std::vector<Spatial::Location> locations;
  for (int row = 0; row < 10; ++row) {
    for (int col = 0; col < 12; ++col) {
      Spatial::Location location;
      location.coordinate = {.latitude = static_cast<float>(row),
                             .longitude = static_cast<float>(col)};
      locations.push_back(location);
    }
  }

  const auto table = GridPairTable::make_distances(locations, 5.0F);
  const size_t dense_bytes = locations.size() * locations.size() * sizeof(double);
  EXPECT_LT(table.memory_bytes(), dense_bytes);
}
