/*
 * GridPairTable.h
 *
 * Compact storage for pair values on a raster grid.
 *
 * A value such as distance depends only on the absolute coordinate difference
 * between two cells, not on their location IDs:
 *
 *   value(from, to) = lookup[abs(from.row - to.row)][abs(from.col - to.col)]
 *
 * The table therefore stores one value per possible grid offset instead of an
 * N x N matrix. The location database owns the coordinates; this class keeps a
 * non-owning reference to it and owns only the offset lookup values.
 */
#ifndef SPATIAL_GIS_GRIDPAIRTABLE_H
#define SPATIAL_GIS_GRIDPAIRTABLE_H

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <vector>

#include "Spatial/Location/Location.h"

class GridPairTable {
public:
  /*
   * Lazy access to all values from one source location.
   *
   * No N-element row is created. The source coordinate is cached once, and
   * operator[] calculates the offset to each requested destination before
   * reading the shared lookup table.
   */
  class ValuesFromLocation {
  public:
    ValuesFromLocation() = default;

    [[nodiscard]] double operator[](size_t to_location) const noexcept {
      const auto &to_coordinate = (*locations_)[to_location].coordinate;
      const auto to_row = static_cast<int32_t>(to_coordinate.latitude);
      const auto to_col = static_cast<int32_t>(to_coordinate.longitude);
      const int32_t d_row = std::abs(from_row_ - to_row);
      const int32_t d_col = std::abs(from_col_ - to_col);
      return values_[(static_cast<size_t>(d_row) * value_cols_) + d_col];
    }

    [[nodiscard]] size_t size() const noexcept { return size_; }

  private:
    friend class GridPairTable;

    ValuesFromLocation(const std::vector<Spatial::Location>* locations,
                       const double* values,
                       int32_t value_cols,
                       int32_t from_row,
                       int32_t from_col,
                       size_t size)
        : locations_(locations),
          values_(values),
          value_cols_(value_cols),
          from_row_(from_row),
          from_col_(from_col),
          size_(size) {}

    const std::vector<Spatial::Location>* locations_{nullptr};
    const double* values_{nullptr};
    int32_t value_cols_{0};
    int32_t from_row_{0};
    int32_t from_col_{0};
    size_t size_{0};
  };

  GridPairTable() = default;

  static GridPairTable make_distances(const std::vector<Spatial::Location> &locations,
                                      float cell_size) {
    GridPairTable table;
    table.locations_ = &locations;
    table.size_ = locations.size();
    if (table.empty()) { return table; }

    int32_t max_row = 0;
    int32_t max_col = 0;
    for (const auto &location : locations) {
      max_row = std::max(static_cast<int32_t>(location.coordinate.latitude), max_row);
      max_col = std::max(static_cast<int32_t>(location.coordinate.longitude), max_col);
    }

    // Raster coordinates are zero-based, so the largest possible difference
    // is also the largest row or column index.
    const int32_t value_rows = max_row + 1;
    table.value_cols_ = max_col + 1;
    table.values_.resize(static_cast<size_t>(value_rows) * table.value_cols_);
    for (int32_t d_row = 0; d_row < value_rows; ++d_row) {
      for (int32_t d_col = 0; d_col < table.value_cols_; ++d_col) {
        table.values_[(static_cast<size_t>(d_row) * table.value_cols_) + d_col] =
            std::sqrt(std::pow(cell_size * static_cast<float>(d_row), 2)
                      + std::pow(cell_size * static_cast<float>(d_col), 2));
      }
    }
    return table;
  }

  template <typename Func>
  [[nodiscard]] GridPairTable map(Func func) const {
    // Derived tables, such as movement kernels, reuse the same locations and
    // transform only the compact offset values.
    GridPairTable table;
    table.locations_ = locations_;
    table.size_ = size_;
    table.value_cols_ = value_cols_;
    table.values_.resize(values_.size());
    for (size_t i = 0; i < values_.size(); ++i) { table.values_[i] = func(values_[i]); }
    return table;
  }

  [[nodiscard]] bool empty() const noexcept { return size_ == 0; }
  [[nodiscard]] size_t size() const noexcept { return size_; }

  [[nodiscard]] ValuesFromLocation values_from(size_t from_location) const noexcept {
    const auto &coordinate = (*locations_)[from_location].coordinate;
    return {locations_,
            values_.data(),
            value_cols_,
            static_cast<int32_t>(coordinate.latitude),
            static_cast<int32_t>(coordinate.longitude),
            size_};
  }

  // Convenient single-pair lookup. Repeated lookups from the same source
  // should retain values_from(from_location) and index that object instead.
  [[nodiscard]] double at(size_t from_location, size_t to_location) const noexcept {
    return values_from(from_location)[to_location];
  }

  [[nodiscard]] const std::vector<double> &stored_values() const noexcept { return values_; }

  [[nodiscard]] size_t memory_bytes() const noexcept { return values_.size() * sizeof(double); }

private:
  // SpatialSettings finalizes its location database before constructing the
  // table and must not change its shape while this non-owning reference is used.
  const std::vector<Spatial::Location>* locations_{nullptr};
  size_t size_{0};
  int32_t value_cols_{0};
  std::vector<double> values_;
};

#endif  // SPATIAL_GIS_GRIDPAIRTABLE_H
