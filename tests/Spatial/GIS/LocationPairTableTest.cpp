#include <gtest/gtest.h>

#include <cmath>
#include <vector>

#include "Spatial/GIS/LocationPairTable.h"
#include "Spatial/Location/Location.h"

TEST(LocationPairTableTest, GridLookupMatchesDenseEuclideanDistances) {
  std::vector<Spatial::Location> locations(4);
  locations[0].coordinate = {.latitude = 0, .longitude = 0};
  locations[1].coordinate = {.latitude = 0, .longitude = 2};
  locations[2].coordinate = {.latitude = 3, .longitude = 0};
  locations[3].coordinate = {.latitude = 3, .longitude = 2};

  constexpr float cell_size = 5.0F;
  const auto table = LocationPairTable::make_grid_distances(locations, cell_size);

  ASSERT_TRUE(table.is_grid());
  ASSERT_EQ(table.size(), locations.size());
  for (size_t from = 0; from < locations.size(); ++from) {
    for (size_t to = 0; to < locations.size(); ++to) {
      const double expected = std::sqrt(
          std::pow(cell_size * (locations[from].coordinate.latitude
                                - locations[to].coordinate.latitude),
                   2)
          + std::pow(cell_size * (locations[from].coordinate.longitude
                                  - locations[to].coordinate.longitude),
                     2));
      EXPECT_DOUBLE_EQ(table.at(static_cast<int>(from), static_cast<int>(to)), expected);
    }
  }
}

TEST(LocationPairTableTest, DenseFallbackPreservesRows) {
  const std::vector<std::vector<double>> matrix = {
      {0.0, 1.5, 2.5}, {1.5, 0.0, 3.5}, {2.5, 3.5, 0.0}};
  const auto table = LocationPairTable::make_dense(matrix);

  ASSERT_FALSE(table.is_grid());
  ASSERT_EQ(table.size(), matrix.size());
  for (size_t row = 0; row < matrix.size(); ++row) {
    EXPECT_EQ(table.row(static_cast<int>(row)), matrix[row]);
  }
}

TEST(LocationPairTableTest, DerivedTableKeepsZeroDistanceSentinel) {
  const auto distances = LocationPairTable::make_dense({{0.0, 2.0}, {2.0, 0.0}});
  const auto squared = distances.map_with_zero_sentinel(
      [](double distance) { return distance * distance; }, "test squared distances");

  EXPECT_DOUBLE_EQ(squared.at(0, 0), 0.0);
  EXPECT_DOUBLE_EQ(squared.at(0, 1), 4.0);
  EXPECT_DOUBLE_EQ(squared.at(1, 0), 4.0);
  EXPECT_DOUBLE_EQ(squared.at(1, 1), 0.0);
}

TEST(LocationPairTableTest, GridStorageIsSmallerThanDenseStorage) {
  std::vector<Spatial::Location> locations;
  for (int row = 0; row < 10; ++row) {
    for (int col = 0; col < 12; ++col) {
      Spatial::Location location;
      location.coordinate = {.latitude = static_cast<float>(row),
                             .longitude = static_cast<float>(col)};
      locations.push_back(location);
    }
  }

  const auto table = LocationPairTable::make_grid_distances(locations, 5.0F);
  const size_t dense_bytes = locations.size() * locations.size() * sizeof(double);
  EXPECT_LT(table.memory_bytes(), dense_bytes);
}
