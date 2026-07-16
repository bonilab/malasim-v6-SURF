#include "WesolowskiSurfaceSM.h"

#include "Simulation/Model.h"

void Spatial::WesolowskiSurfaceSM::prepare() {
  if (spatial_distance_ == nullptr) {
    throw std::runtime_error(fmt::format("{} called without spatial distances", __FUNCTION__));
  }

  const double gamma = gamma_;
  distance_power_ = spatial_distance_->map_with_zero_sentinel(
      [gamma](double distance) { return std::pow(distance, gamma); },
      "WesolowskiSurfaceSM distance powers");

  AscFile* travel_raster =
      Model::get_spatial_data()->get_raster(SpatialData::SpatialFileType::TRAVEL);
  travel = std::move(prepare_surface(travel_raster, number_of_locations_));
}
