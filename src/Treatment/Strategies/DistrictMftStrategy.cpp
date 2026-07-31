/*
 * DistrictMftStrategy.cpp
 *
 * Implement the district MFT strategy.
 */
#include "DistrictMftStrategy.h"

#include "Population/Person/Person.h"
#include "Simulation/Model.h"
#include "Spatial/GIS/SpatialData.h"
#include "Utils/Random.h"

namespace {
// The administrative level this strategy keys off. It is resolved to a numeric
// level id once, in the constructor, so that get_therapy never pays for a
// string comparison.
constexpr auto K_DISTRICT_LEVEL = "district";
}  // namespace

DistrictMftStrategy::DistrictMftStrategy() : IStrategy("DistrictMFT", StrategyType::DistrictMft) {
  const auto* spatial_data = Model::get_spatial_data();
  const auto* boundary = spatial_data->get_boundary(K_DISTRICT_LEVEL);
  if (boundary == nullptr) {
    throw std::runtime_error(
        "The DistrictMFT strategy requires an administrative level named 'district' under "
        "spatial_settings.<mode>.administrative_boundaries.");
  }

  // Cache the numeric level id so the per-treatment lookup in get_therapy can
  // use the integer overload of get_admin_unit.
  district_level_id_ = spatial_data->get_admin_level_id(K_DISTRICT_LEVEL);

  // Size to accommodate either 0-based or 1-based district ids. When the raster
  // is 1-based, slot 0 is simply never referenced.
  district_strategies_.resize(static_cast<std::size_t>(boundary->max_unit_id) + 1);
}

void DistrictMftStrategy::add_therapy(Therapy* therapy) {
  throw std::runtime_error("Invalid function called to add therapy to the District MFT Strategy.");
}

void DistrictMftStrategy::set_district_strategy(int district,
                                                std::unique_ptr<MftStrategy> strategy) {
  // Validate district ID is within bounds
  if (district < 0 || static_cast<std::size_t>(district) >= district_strategies_.size()) {
    throw std::out_of_range(fmt::format("District ID {} is out of valid range [0, {}]", district,
                                        district_strategies_.size() - 1));
  }

  // Check if district already has a strategy assigned
  if (district_strategies_[district] != nullptr) {
    throw std::runtime_error(
        fmt::format("District {} already has an MFT strategy assigned", district));
  }

  district_strategies_[district] = std::move(strategy);
}

Therapy* DistrictMftStrategy::get_therapy(Person* person) {
  // Resolve the MFT for this district. Both lookups here are O(1): the cached
  // level id indexes a vector inside AdminLevelManager, and the district id
  // indexes district_strategies_ directly.
  const auto district = Model::get_spatial_data()->get_admin_unit(district_level_id_,
                                                                  person->get_location());

  const auto* mft = (district >= 0
                     && static_cast<std::size_t>(district) < district_strategies_.size())
                        ? district_strategies_[district].get()
                        : nullptr;
  if (mft == nullptr || mft->therapies.empty()) {
    throw std::runtime_error(
        fmt::format("No MFT assigned to district {} in strategy {}", district, this->name));
  }

  // Select the therapy to give the individual
  auto pr = Model::get_random()->random_flat(0.0, 1.0);
  auto sum = 0.0;
  for (std::size_t ndx = 0; ndx < mft->percentages.size(); ndx++) {
    sum += mft->percentages[ndx];
    if (pr <= sum) { return Model::get_therapy_db()[mft->therapies[ndx]].get(); }
  }

  // The distribution is validated to sum to 1.0 within tolerance at load time,
  // so reaching this point means pr fell inside the floating-point residual at
  // the top of the range. That residual belongs to the final bucket; previously
  // this threw and aborted the run.
  return Model::get_therapy_db()[mft->therapies.back()].get();
}
