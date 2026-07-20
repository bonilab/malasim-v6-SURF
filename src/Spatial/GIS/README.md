# Geographic Information System (GIS)

This module implements the geographic information system components for spatial data management and analysis in the malaria simulation.

## Overview

The GIS module provides:

- Administrative level management
- Spatial data handling
- Geographic data processing
- ASC file operations
- Location mapping
- Interchangeable distance storage

## Components

### Core Systems

- `AdminLevelManager`: Administrative hierarchy
  - Level definitions
  - Boundary management
  - Region relationships
  - Location mapping
  - Hierarchy traversal

### Spatial Data

- `SpatialData`: Geographic data management
  - Coordinate systems
  - Region definitions
  - Distance calculations
  - Area computations
  - Spatial relationships
- `DistanceProvider`: Common distance interface for both spatial modes
  - `DenseDistanceProvider` owns location-based haversine distances
  - `DenseGridDistanceProvider` stores the reference raster distance matrix
  - `GridLutDistanceProvider` stores distances by absolute row/column delta
  - All providers expose pair lookup and their number of locations
- `GridPairTable`: Compact raster-offset storage used by `GridLutDistanceProvider` and grid movement
  kernels

### File Operations

- `AscFile`: ASC file handling
  - Grid data reading
  - Raster processing
  - Data validation
  - Format conversion
  - File I/O

## Implementation

### Admin Level Management

```cpp
class AdminLevelManager {
    void add_level(int level, string name);
    void map_locations(int level);
    vector<int> get_locations(int level);
    bool is_ancestor_of(int location1, int location2);
};
```

### Distance Providers

```cpp
class DistanceProvider {
public:
  virtual ~DistanceProvider() = default;
  virtual double distance(std::size_t from_location,
                          std::size_t to_location) const noexcept = 0;
  virtual std::size_t size() const noexcept = 0;

  // Non-null only for a provider backed by a dense matrix.
  virtual const std::vector<double>* dense_row(
      std::size_t from_location) const noexcept;

  // Non-null only for the raster LUT provider.
  virtual const GridPairTable* grid_table() const noexcept;
};
```

`SpatialSettings` owns one `DistanceProvider`. The concrete provider depends on the configured
spatial mode:

| Spatial mode | Provider | Distance calculation |
|---|---|---|
| `location_based` | `DenseDistanceProvider` | Haversine distance between latitude/longitude coordinates |
| `grid_based`, production | `GridLutDistanceProvider` | Euclidean distance indexed by raster-cell offset |
| `grid_based`, verification | `DenseGridDistanceProvider` | Euclidean distance stored as a full matrix |

### Grid Pair Table

Raster locations store their cell `(row, column)` in `Location::coordinate`. For grid distances and
derived movement values, the location IDs themselves do not matter. Only the absolute cell offset
matters:

```text
value(from, to) = lookup[abs(from.row - to.row)][abs(from.col - to.col)]
```

`GridPairTable` therefore owns an offset lookup of approximately `raster_rows × raster_columns`
values instead of an `N × N` location matrix. It does not copy row and column arrays. It reads
coordinates from the existing location database through a non-owning reference, so the location
database must remain alive and must not change shape while the table is in use.

For a single pair, use `at()`:

```cpp
const double distance = table.at(from_location, to_location);
```

For repeated lookups from one source, retain `values_from()`:

```cpp
const auto distances = table.values_from(from_location);
for (std::size_t to_location = 0; to_location < distances.size(); ++to_location) {
  const double distance = distances[to_location];
}
```

`ValuesFromLocation` is a lazy view. It does not allocate or materialize a matrix row. It caches the
source cell, reads each destination cell when indexed, calculates `(delta_row, delta_column)`, and
returns the corresponding lookup value.

`GridPairTable::map()` transforms only the stored offset values and retains the same location
database. Movement models use this to precompute compact grid kernels.

### Raster Backend Selection

`SpatialData::generate_distances()` creates the selected raster provider after locations have been
generated. Backend selection is deliberately source-controlled while exact equivalence is being
verified:

```cpp
enum class GridDistanceBackend { Dense, Lut };
inline constexpr GridDistanceBackend GRID_DISTANCE_BACKEND = GridDistanceBackend::Lut;
```

There is no CMake option or YAML setting for this selection.

## Usage

### Administrative Management

```cpp
// Initialize admin levels
auto admin = new AdminLevelManager();
admin->add_level(0, "Country");
admin->add_level(1, "Province");

// Map locations
admin->map_locations(1);
auto locations = admin->get_locations(1);
```

### Comparing Raster Distances

```cpp
DenseGridDistanceProvider dense(locations, cell_size);
GridLutDistanceProvider lut(locations, cell_size);

for (std::size_t from = 0; from < locations.size(); ++from) {
  for (std::size_t to = 0; to < locations.size(); ++to) {
    assert(dense.distance(from, to) == lut.distance(from, to));
  }
}
```

The unit tests use `EXPECT_DOUBLE_EQ` for every location pair. Exact equality is required because the
LUT is intended to replace the dense raster matrix without changing movement probabilities or random
number streams.

## Distance Behavior by Spatial Mode

### Grid-Based Mode

Raster coordinates store integer row and column indices. Distance is Euclidean and depends only on
`(|delta_row|, |delta_column|)`, so an `N × N` dense matrix can be represented by a much smaller lookup
table. The current hard-coded production selection is `GridDistanceBackend::Lut`; changing the constant
to `Dense` selects the reference implementation on the next build. Both provider classes coexist in
each build, so tests can instantiate and compare them directly.

### Location-Based Mode

The grid LUT does not apply to location-based configurations. These locations contain arbitrary
latitude/longitude coordinates and use `DenseDistanceProvider`, populated with haversine distances
from `Spatial::Coordinate::calculate_distance_in_km`. During circulation, `Population` passes the
provider's existing dense row to the movement model.

### Grid-Based Movement

Grid movement models build a compact movement kernel from the provider's `GridPairTable` during
`prepare()`. During circulation they call `values_from(from_location)` once and reuse that view for
every destination. No dense distance row is generated.

## Key Features

### Administrative Management

- Level hierarchy
- Location mapping
- Region relationships
- Boundary management
- Location queries

### Spatial Analysis

- Distance calculations
- Bit-identical dense/LUT raster distance implementations
- Compact grid movement kernels
- Neighbor finding
- Region containment
- Area calculations
- Coordinate transformations

### Data Processing

- Raster handling
- Grid operations
- Data validation
- Format conversion
- File management

## Dependencies

### Core Components

- `Location`
- `Coordinate`
- `Config`
- `Model`

### Support Systems

- File I/O
- Math utilities
- Data structures
- Error handling

## Notes

### Implementation Guidelines

- Validate spatial data
- Handle edge cases
- Check boundaries
- Verify relationships
- Document assumptions
- Test thoroughly
- Handle errors
- Log operations

### Performance Considerations

- Cache calculations
- Prefer the raster LUT for large grids to avoid quadratic distance storage
- Use `values_from()` when reading many destinations from the same source
- Keep the location database stable while a `GridPairTable` references it
- Keep the dense raster provider as the numerical reference until replacement verification is complete
- Optimize queries
- Manage memory
- Handle large datasets
- Profile operations
- Validate inputs
- Clean up resources
- Document limitations
