/*
 * LocationPairTable.h
 *
 * Compact storage for any quantity that is a pure function of a pair of
 * locations (the spatial distance matrix, and the movement-model kernels
 * derived from it).
 *
 * Rationale
 * ---------
 * A grid-based configuration creates exactly one location per non-NODATA raster
 * cell, and SpatialData::generate_locations() stores that cell's (row, col)
 * directly into Location::coordinate. The Euclidean distance between two
 * locations therefore depends only on (|dRow|, |dCol|), so the full n*n matrix
 * collapses to an (nRows x nCols) lookup table with *bit-identical* values.
 *
 * For a 277x308 raster with 50,745 valid cells that is 20.6 GB -> ~635 KB, and
 * the table then fits in cache instead of missing on every access.
 *
 * Location-based configurations use arbitrary lat/lon with a haversine distance
 * (Coordinate::calculate_distance_in_km), which does not quantize onto a grid.
 * Those fall back to the original dense n*n matrix; such configurations have few
 * locations, so the memory cost is irrelevant.
 */
#ifndef SPATIAL_LOCATIONPAIRTABLE_H
#define SPATIAL_LOCATIONPAIRTABLE_H

#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <stdexcept>
#include <vector>

#include "Spatial/Location/Location.h"

class LocationPairTable {
public:
  LocationPairTable() = default;

  /*
   * Build the grid distance table.
   *
   * The per-entry expression is character-for-character the one previously used
   * by SpatialData::generate_distances(), including the `float` type of
   * cell_size and of the coordinate delta, so results are bit-identical and RNG
   * streams are unaffected.
   */
  static LocationPairTable make_grid_distances(const std::vector<Spatial::Location> &location_db,
                                               float cell_size) {
    LocationPairTable table;
    table.grid_mode_ = true;
    table.size_ = location_db.size();
    if (table.size_ == 0) { return table; }

    table.grid_row_.resize(table.size_);
    table.grid_col_.resize(table.size_);

    int32_t max_row = 0;
    int32_t max_col = 0;
    for (size_t i = 0; i < table.size_; ++i) {
      // Coordinates are raster (row, col) indices held in a float; exact for any
      // realistic raster (< 2^24 rows/cols).
      table.grid_row_[i] = static_cast<int32_t>(location_db[i].coordinate.latitude);
      table.grid_col_[i] = static_cast<int32_t>(location_db[i].coordinate.longitude);
      if (table.grid_row_[i] > max_row) { max_row = table.grid_row_[i]; }
      if (table.grid_col_[i] > max_col) { max_col = table.grid_col_[i]; }
    }

    table.lut_rows_ = max_row + 1;
    table.lut_cols_ = max_col + 1;
    table.lut_.resize(static_cast<size_t>(table.lut_rows_) * table.lut_cols_);
    for (int32_t d_row = 0; d_row < table.lut_rows_; ++d_row) {
      for (int32_t d_col = 0; d_col < table.lut_cols_; ++d_col) {
        table.lut_[(static_cast<size_t>(d_row) * table.lut_cols_) + d_col] =
            std::sqrt(std::pow(cell_size * static_cast<float>(d_row), 2)
                      + std::pow(cell_size * static_cast<float>(d_col), 2));
      }
    }
    return table;
  }

  // Dense fallback, used by the location-based path.
  static LocationPairTable make_dense(std::vector<std::vector<double>> matrix) {
    LocationPairTable table;
    table.grid_mode_ = false;
    table.size_ = matrix.size();
    table.dense_ = std::move(matrix);
    return table;
  }

  /*
   * Derive a second table whose values are func(distance), preserving the
   * compact representation. Used for the Marshall / BurkinaFaso kernels, which
   * are pure functions of distance and so collapse identically.
   */
  template <typename Func>
  [[nodiscard]] LocationPairTable map(Func func) const {
    LocationPairTable table;
    table.grid_mode_ = grid_mode_;
    table.size_ = size_;
    if (grid_mode_) {
      table.grid_row_ = grid_row_;
      table.grid_col_ = grid_col_;
      table.lut_rows_ = lut_rows_;
      table.lut_cols_ = lut_cols_;
      table.lut_.resize(lut_.size());
      for (size_t i = 0; i < lut_.size(); ++i) { table.lut_[i] = func(lut_[i]); }
    } else {
      table.dense_.resize(dense_.size());
      for (size_t i = 0; i < dense_.size(); ++i) {
        table.dense_[i].resize(dense_[i].size());
        for (size_t j = 0; j < dense_[i].size(); ++j) { table.dense_[i][j] = func(dense_[i][j]); }
      }
    }
    return table;
  }

  [[nodiscard]] bool empty() const { return size_ == 0; }
  [[nodiscard]] size_t size() const { return size_; }
  [[nodiscard]] bool is_grid() const { return grid_mode_; }

  [[nodiscard]] double at(int from, int to) const {
    if (grid_mode_) {
      const int32_t d_row = std::abs(grid_row_[from] - grid_row_[to]);
      const int32_t d_col = std::abs(grid_col_[from] - grid_col_[to]);
      return lut_[(static_cast<size_t>(d_row) * lut_cols_) + d_col];
    }
    return dense_[from][to];
  }

  /*
   * Materialise row `from`.
   *
   * WARNING: in grid mode the returned reference points at an internal scratch
   * buffer and is invalidated by the next call to row() on this object. Callers
   * only ever need one row at a time (Model::daily_update drives circulation
   * location-by-location, serially), but do not retain the reference. Not
   * thread-safe; if circulation is ever parallelised, give each worker its own
   * buffer.
   */
  [[nodiscard]] const std::vector<double> &row(int from) const {
    if (!grid_mode_) { return dense_[from]; }
    row_buffer_.resize(size_);
    const int32_t from_row = grid_row_[from];
    const int32_t from_col = grid_col_[from];
    for (size_t to = 0; to < size_; ++to) {
      const int32_t d_row = std::abs(from_row - grid_row_[to]);
      const int32_t d_col = std::abs(from_col - grid_col_[to]);
      row_buffer_[to] = lut_[(static_cast<size_t>(d_row) * lut_cols_) + d_col];
    }
    return row_buffer_;
  }

  // Resident footprint, excluding the scratch row buffer.
  [[nodiscard]] size_t memory_bytes() const {
    if (grid_mode_) {
      return (lut_.size() * sizeof(double)) + (grid_row_.size() * sizeof(int32_t))
             + (grid_col_.size() * sizeof(int32_t));
    }
    size_t bytes = 0;
    for (const auto &row : dense_) { bytes += row.size() * sizeof(double); }
    return bytes;
  }

private:
  bool grid_mode_{false};
  size_t size_{0};

  // Grid mode: values indexed by (|dRow|, |dCol|), plus each location's cell.
  int32_t lut_rows_{0};
  int32_t lut_cols_{0};
  std::vector<double> lut_;
  std::vector<int32_t> grid_row_;
  std::vector<int32_t> grid_col_;

  // Location-based mode: the original dense matrix.
  std::vector<std::vector<double>> dense_;

  mutable std::vector<double> row_buffer_;
};

#endif  // SPATIAL_LOCATIONPAIRTABLE_H
