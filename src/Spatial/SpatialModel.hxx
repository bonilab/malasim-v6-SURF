/*
 * SpatialModel.hxx
 *
 * Base abstract class for the spatial movement models that are implemented in
 * the simulation.
 */
#ifndef SPATIAL_SPATIALMODEL_H
#define SPATIAL_SPATIALMODEL_H

#include <spdlog/spdlog.h>

#include "Spatial/GIS/AscFile.h"
#include "Utils/TypeDef.h"

namespace Spatial {
class SpatialModel {
public:
  // Disallow copy
  SpatialModel(const SpatialModel&) = delete;
  SpatialModel& operator=(const SpatialModel&) = delete;

  // Disallow move
  SpatialModel(SpatialModel&&) = delete;
  SpatialModel& operator=(SpatialModel&&) = delete;

protected:
  // Prepare the travel raster for the movement model
  static std::vector<double> prepare_surface(const AscFile* travel_raster,
                                             int number_of_locations) {
    // Get the travel times raster
    spdlog::info("Preparing travel surface...");
    if (travel_raster == nullptr) {
      throw std::runtime_error(fmt::format("{} called without travel data loaded", __FUNCTION__));
    }
    // Initialize vector with the correct size
    std::vector<double> travel;
    travel.reserve(number_of_locations);

    // Use the min and max to normalize the raster into a vector
    for (auto row = 0; row < travel_raster->nrows; row++) {
      for (auto col = 0; col < travel_raster->ncols; col++) {
        if (travel_raster->data[row][col] == travel_raster->nodata_value) { continue; }
        travel.push_back(travel_raster->data[row][col]);
      }
    }

    // Return the pointer to the array
    return travel;
  }

public:
  SpatialModel() = default;

  virtual ~SpatialModel() = default;

  // Allow the spatial model to perform any preparation it must do.
  virtual void prepare() {}

  [[nodiscard]] virtual DoubleVector get_v_relative_out_movement_to_destination(
      const int &from_location, const int &number_of_locations,
      const IntVector &v_number_of_residents_by_location) const = 0;
  ;
};
}  // namespace Spatial

#endif
