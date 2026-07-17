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
 * collapses to an (nRows x nCols) lookup table with bit-identical values.
 *
 * For a 277x308 raster with 50,745 valid cells, the distance values shrink from
 * about 20.6 GB to about 683 KB. The location-to-grid mappings add about 406 KB.
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
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "Spatial/Location/Location.h"

class LocationPairTable {
public:
  /*
   * Lightweight, non-owning access to one source row.
   *
   * In grid mode this stores the source cell once, then computes only the
   * destination offset for operator[]. This avoids materialising a temporary
   * N-element row and then reading it again in the movement model.
   *
   * A RowView remains valid while its LocationPairTable remains alive and is not
   * modified. It owns no memory and is safe for the current serial circulation
   * loop.
   */
  class RowView {
  public:
    RowView() = default;

    [[nodiscard]] double operator[](size_t to) const noexcept {
      // Discriminates on the mode, not on dense_row_ != nullptr: a dense row of
      // length zero yields data() == nullptr, which would silently divert to the
      // grid branch and dereference a null lut_.
      if (!grid_) { return dense_row_[to]; }

      const int32_t d_row = std::abs(from_row_ - grid_row_[to]);
      const int32_t d_col = std::abs(from_col_ - grid_col_[to]);
      return lut_[(static_cast<size_t>(d_row) * lut_cols_) + d_col];
    }

    [[nodiscard]] size_t size() const noexcept { return size_; }

  private:
    friend class LocationPairTable;

    RowView(const double* lut, const int32_t* grid_row, const int32_t* grid_col,
            int32_t lut_cols, int32_t from_row, int32_t from_col, size_t size)
        : lut_(lut),
          grid_row_(grid_row),
          grid_col_(grid_col),
          lut_cols_(lut_cols),
          from_row_(from_row),
          from_col_(from_col),
          size_(size),
          grid_(true) {}

    RowView(const double* dense_row, size_t size) : dense_row_(dense_row), size_(size) {}

    const double* lut_{nullptr};
    const int32_t* grid_row_{nullptr};
    const int32_t* grid_col_{nullptr};
    const double* dense_row_{nullptr};
    int32_t lut_cols_{0};
    int32_t from_row_{0};
    int32_t from_col_{0};
    size_t size_{0};
    bool grid_{false};
  };

  LocationPairTable() = default;

  /*
   * Build the grid distance table.
   *
   * The per-entry expression is character-for-character the one previously used
   * by SpatialData::generate_distances(), including the `float` type of
   * cell_size and of the coordinate delta, so results are bit-identical and RNG
   * streams are unaffected.
   */
  static LocationPairTable make_grid_distances(const std::vector<Spatial::Location>& location_db,
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
   * compact representation. Used to precompute movement-model distance terms.
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

  /*
   * Derive func(distance) with every zero-distance pair forced to exactly 0.0,
   * so a movement model can use `NumberHelpers::is_zero(value)` as its skip test
   * and never touch the distance row at all.
   *
   * The sentinel is only sound while no pair at a genuinely non-zero distance
   * maps to a value that is_zero() would also catch. If that ever happened the
   * model would silently drop a destination the unpatched code kept, feeding
   * random_multinomial a different weight vector and desynchronising the RNG
   * stream from that day on -- a divergence that would surface as unexplained
   * trajectory differences, not as a crash. So check it here, once, at startup,
   * and refuse to run rather than corrupt the stream.
   *
   * The check is conservative in grid mode: it covers every (dRow, dCol) the
   * table can express, including deltas no actual pair of locations realises.
   */
  template <typename Func>
  [[nodiscard]] LocationPairTable map_with_zero_sentinel(Func func,
                                                         const std::string &context) const {
    constexpr double kEpsilon = std::numeric_limits<double>::epsilon();
    const auto reads_as_zero = [](double value) { return std::fabs(value) < kEpsilon; };

    LocationPairTable table = map([&](double distance) -> double {
      if (reads_as_zero(distance)) { return 0.0; }
      return func(distance);
    });

    const auto verify = [&](double distance, double value) {
      if (reads_as_zero(distance) || !reads_as_zero(value)) { return; }
      std::ostringstream message;
      message << context << ": the zero sentinel is ambiguous. A pair at distance " << distance
              << " maps to " << value
              << ", which is_zero() cannot distinguish from the zero-distance sentinel (epsilon "
              << kEpsilon
              << "). The movement model would silently drop that destination and desynchronise "
                 "the RNG stream. Check the model parameters against the raster extent.";
      throw std::runtime_error(message.str());
    };

    if (grid_mode_) {
      for (size_t i = 0; i < lut_.size(); ++i) { verify(lut_[i], table.lut_[i]); }
    } else {
      for (size_t i = 0; i < dense_.size(); ++i) {
        for (size_t j = 0; j < dense_[i].size(); ++j) { verify(dense_[i][j], table.dense_[i][j]); }
      }
    }
    return table;
  }

  [[nodiscard]] bool empty() const { return size_ == 0; }
  [[nodiscard]] size_t size() const { return size_; }
  [[nodiscard]] bool is_grid() const { return grid_mode_; }

  [[nodiscard]] RowView row_view(int from) const noexcept {
    if (grid_mode_) {
      return RowView(lut_.data(), grid_row_.data(), grid_col_.data(), lut_cols_, grid_row_[from],
                     grid_col_[from], size_);
    }
    return RowView(dense_[from].data(), dense_[from].size());
  }

  [[nodiscard]] double at(int from, int to) const noexcept {
    return row_view(from)[static_cast<size_t>(to)];
  }

  /*
   * Materialise row `from` for compatibility and testing.
   *
   * Performance-sensitive movement code should use row_view(), which avoids the
   * extra full-row generation, write, and reread.
   */
  [[nodiscard]] const std::vector<double>& row(int from) const {
    if (!grid_mode_) { return dense_[from]; }
    row_buffer_.resize(size_);
    const auto view = row_view(from);
    for (size_t to = 0; to < size_; ++to) { row_buffer_[to] = view[to]; }
    return row_buffer_;
  }

  // Resident footprint, excluding the compatibility scratch row buffer.
  [[nodiscard]] size_t memory_bytes() const {
    if (grid_mode_) {
      return (lut_.size() * sizeof(double)) + (grid_row_.size() * sizeof(int32_t))
             + (grid_col_.size() * sizeof(int32_t));
    }
    size_t bytes = 0;
    for (const auto& row : dense_) { bytes += row.size() * sizeof(double); }
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
