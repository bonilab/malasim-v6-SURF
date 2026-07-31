#include "StrategyBuilder.h"

#include <algorithm>
#include <cmath>

#include "AdaptiveCyclingStrategy.h"
#include "CyclingStrategy.h"
#include "DistrictMftStrategy.h"
#include "IStrategy.h"
#include "MFTAgeBasedStrategy.h"
#include "MFTMultiLocationStrategy.h"
#include "MFTRebalancingStrategy.h"
#include "MFTStrategy.h"
#include "NestedMFTMultiLocationStrategy.h"
#include "NestedMFTStrategy.h"
#include "NovelDrugIntroductionStrategy.h"
#include "PublicPrivateMultiLocationStrategy.h"
#include "PublicPrivateStrategy.h"
#include "SFTStrategy.h"
#include "Simulation/Model.h"

StrategyBuilder::StrategyBuilder() = default;

StrategyBuilder::~StrategyBuilder() = default;

std::unique_ptr<IStrategy> StrategyBuilder::build(const YAML::Node &ns, const int &strategy_id) {
  const auto type_name = ns["type"].as<std::string>();
  // NOTE: do not use operator[] here. std::map::operator[] default-inserts a
  // zero-valued entry for an unrecognised name, which silently routes the
  // strategy to StrategyType::SFT (0) instead of reporting the bad config.
  const auto entry = IStrategy::strategy_type_map.find(type_name);
  if (entry == IStrategy::strategy_type_map.end()) {
    throw std::invalid_argument("Unknown strategy type in configuration: " + type_name);
  }
  const auto type = entry->second;
  switch (type) {
    case IStrategy::SFT:
      return build_sft_strategy(ns, strategy_id);
    case IStrategy::Cycling:
      return build_cycling_strategy(ns, strategy_id);
    case IStrategy::AdaptiveCycling:
      return build_adaptive_cycling_strategy(ns, strategy_id);
    case IStrategy::MFT:
      return build_mft_strategy(ns, strategy_id);
    case IStrategy::MFTRebalancing:
      return build_mft_rebalancing_strategy(ns, strategy_id);
    case IStrategy::NestedMFT:
      return build_nested_switching_strategy(ns, strategy_id);
    case IStrategy::MFTMultiLocation:
      return build_mft_multi_location_strategy(ns, strategy_id);
    case IStrategy::NestedMFTMultiLocation:
      return build_nested_mft_different_distribution_by_location_strategy(ns, strategy_id);
    case IStrategy::NovelDrugIntroduction:
      return build_novel_drug_introduction_strategy(ns, strategy_id);
    case IStrategy::DistrictMft:
      return build_district_mft_strategy(ns, strategy_id);
    case IStrategy::MFTAgeBased:
      return build_mft_age_based_strategy(ns, strategy_id);
    case IStrategy::PublicPrivate:
      return build_public_private_strategy(ns, strategy_id);
    case IStrategy::PublicPrivateMultiLocation:
      return build_public_private_multi_location_strategy(ns, strategy_id);
    default:
      throw std::invalid_argument("Strategy type has no builder registered: " + type_name);
  }
}

namespace {

double read_public_share(const YAML::Node &node, const char* field_name) {
  if (!node[field_name]) {
    throw std::invalid_argument(std::string("Missing required strategy field: ") + field_name);
  }
  const auto share = node[field_name].as<double>();
  if (!std::isfinite(share) || share < 0.0 || share > 1.0) {
    throw std::invalid_argument(std::string(field_name) + " must be finite and in [0,1]");
  }
  return share;
}

IStrategy* read_child_strategy(const YAML::Node &node, const char* field_name,
                               const int strategy_id) {
  if (!node[field_name]) {
    throw std::invalid_argument(std::string("Missing required strategy field: ") + field_name);
  }
  const auto child_id = node[field_name].as<int>();
  if (child_id < 0 || child_id >= strategy_id
      || static_cast<std::size_t>(child_id) >= Model::get_strategy_db().size()) {
    throw std::invalid_argument(std::string(field_name)
                                + " must reference a previously constructed strategy");
  }
  return Model::get_strategy_db()[child_id].get();
}

int read_peak_after(const YAML::Node &node) {
  if (!node["peak_after"]) {
    throw std::invalid_argument("Missing required strategy field: peak_after");
  }
  const auto peak_after = node["peak_after"].as<int>();
  if (peak_after < 0) { throw std::invalid_argument("peak_after must be non-negative"); }
  return peak_after;
}

DoubleVector read_location_shares(const YAML::Node &node, const char* field_name,
                                  const int number_of_locations) {
  if (!node[field_name] || !node[field_name].IsSequence()) {
    throw std::invalid_argument(std::string("Missing required strategy sequence: ") + field_name);
  }
  if (node[field_name].size() != static_cast<std::size_t>(number_of_locations)) {
    throw std::invalid_argument(std::string(field_name)
                                + " must contain exactly one value per location");
  }

  DoubleVector shares;
  shares.reserve(node[field_name].size());
  for (const auto &entry : node[field_name]) {
    const auto share = entry.as<double>();
    if (!std::isfinite(share) || share < 0.0 || share > 1.0) {
      throw std::invalid_argument(std::string(field_name) + " values must be finite and in [0,1]");
    }
    shares.push_back(share);
  }
  return shares;
}

void configure_child_strategies(const YAML::Node &node, const int strategy_id,
                                PublicPrivateStrategy* strategy) {
  if (!node["public_strategy_id"] || !node["private_strategy_id"]) {
    throw std::invalid_argument("Public and private strategy IDs are required");
  }
  const auto public_id = node["public_strategy_id"].as<int>();
  const auto private_id = node["private_strategy_id"].as<int>();
  if (public_id == private_id) {
    throw std::invalid_argument("Public and private strategy IDs must be distinct");
  }
  strategy->set_public_strategy(read_child_strategy(node, "public_strategy_id", strategy_id));
  strategy->set_private_strategy(read_child_strategy(node, "private_strategy_id", strategy_id));
}

void configure_child_strategies(const YAML::Node &node, const int strategy_id,
                                PublicPrivateMultiLocationStrategy* strategy) {
  if (!node["public_strategy_id"] || !node["private_strategy_id"]) {
    throw std::invalid_argument("Public and private strategy IDs are required");
  }
  const auto public_id = node["public_strategy_id"].as<int>();
  const auto private_id = node["private_strategy_id"].as<int>();
  if (public_id == private_id) {
    throw std::invalid_argument("Public and private strategy IDs must be distinct");
  }
  strategy->set_public_strategy(read_child_strategy(node, "public_strategy_id", strategy_id));
  strategy->set_private_strategy(read_child_strategy(node, "private_strategy_id", strategy_id));
}

}  // namespace

void StrategyBuilder::add_therapies(const YAML::Node &ns, IStrategy* result) {
  for (auto i = 0; i < ns["therapy_ids"].size(); i++) {
    result->add_therapy(Model::get_therapy_db()[ns["therapy_ids"][i].as<int>()].get());
  }
}

void StrategyBuilder::add_distributions(const YAML::Node &ns, DoubleVector &distribution) {
  for (const auto &node : ns) { distribution.push_back(node.as<double>()); }
}

std::unique_ptr<IStrategy> StrategyBuilder::build_sft_strategy(const YAML::Node &ns,
                                                               const int &strategy_id) {
  auto result = std::make_unique<SFTStrategy>();
  result->id = strategy_id;
  result->name = ns["name"].as<std::string>();
  result->add_therapy(Model::get_therapy_db()[ns["therapy_ids"][0].as<int>()].get());
  return result;
}

std::unique_ptr<IStrategy> StrategyBuilder::build_cycling_strategy(const YAML::Node &ns,
                                                                   const int &strategy_id) {
  auto result = std::make_unique<CyclingStrategy>();
  result->id = strategy_id;
  result->name = ns["name"].as<std::string>();

  result->cycling_time = ns["cycling_time"].as<int>();
  result->next_switching_day = ns["cycling_time"].as<int>();

  add_therapies(ns, result.get());

  return result;
}

std::unique_ptr<IStrategy> StrategyBuilder::build_adaptive_cycling_strategy(
    const YAML::Node &ns, const int &strategy_id) {
  auto result = std::make_unique<AdaptiveCyclingStrategy>();
  result->id = strategy_id;
  result->name = ns["name"].as<std::string>();

  result->trigger_value = ns["trigger_value"].as<double>();
  result->delay_until_actual_trigger = ns["delay_until_actual_trigger"].as<int>();
  result->turn_off_days = ns["turn_off_days"].as<int>();

  add_therapies(ns, result.get());
  return result;
}

std::unique_ptr<IStrategy> StrategyBuilder::build_mft_strategy(const YAML::Node &ns,
                                                               const int &strategy_id) {
  auto result = std::make_unique<MFTStrategy>();
  result->id = strategy_id;
  result->name = ns["name"].as<std::string>();

  // Check for both "distributions" (plural, used in tests) and "distribution" (singular)
  if (ns["distributions"]) {
    add_distributions(ns["distributions"], result->distribution);
  } else if (ns["distribution"]) {
    add_distributions(ns["distribution"], result->distribution);
  }

  add_therapies(ns, result.get());
  return result;
}

std::unique_ptr<IStrategy> StrategyBuilder::build_nested_switching_strategy(
    const YAML::Node &ns, const int &strategy_id) {
  auto result = std::make_unique<NestedMFTStrategy>();
  result->id = strategy_id;
  result->name = ns["name"].as<std::string>();

  add_distributions(ns["start_distribution"], result->start_distribution);
  add_distributions(ns["start_distribution"], result->distribution);
  add_distributions(ns["peak_distribution"], result->peak_distribution);

  result->peak_after = ns["peak_after"].as<int>();

  for (int i = 0; i < ns["strategy_ids"].size(); i++) {
    result->add_strategy(Model::get_strategy_db()[ns["strategy_ids"][i].as<int>()].get());
  }

  return result;
}

std::unique_ptr<IStrategy> StrategyBuilder::build_mft_rebalancing_strategy(const YAML::Node &ns,
                                                                           const int &strategy_id) {
  auto result = std::make_unique<MFTRebalancingStrategy>();
  result->id = strategy_id;
  result->name = ns["name"].as<std::string>();

  add_distributions(ns["distribution"], result->distribution);
  add_distributions(ns["distribution"], result->next_distribution);

  add_therapies(ns, result.get());

  result->update_duration_after_rebalancing = ns["update_duration_after_rebalancing"].as<int>();
  result->delay_until_actual_trigger = ns["delay_until_actual_trigger"].as<int>();
  result->latest_adjust_distribution_time = 0;

  return result;
}

std::unique_ptr<IStrategy> StrategyBuilder::build_mft_multi_location_strategy(
    const YAML::Node &ns, const int &strategy_id) {
  auto result = std::make_unique<MFTMultiLocationStrategy>();
  result->id = strategy_id;
  result->name = ns["name"].as<std::string>();

  result->distribution.clear();
  result->distribution.resize(static_cast<uint64_t>(Model::get_config()->number_of_locations()));

  result->start_distribution.clear();
  result->start_distribution.resize(
      static_cast<uint64_t>(Model::get_config()->number_of_locations()));

  result->peak_distribution.clear();
  result->peak_distribution.resize(
      static_cast<uint64_t>(Model::get_config()->number_of_locations()));

  for (auto loc = 0; loc < Model::get_config()->number_of_locations(); loc++) {
    auto input_loc =
        ns["start_distribution_by_location"].size() < Model::get_config()->number_of_locations()
            ? 0
            : loc;
    add_distributions(ns["start_distribution_by_location"][input_loc], result->distribution[loc]);
  }
  for (auto loc = 0; loc < Model::get_config()->number_of_locations(); loc++) {
    auto input_loc =
        ns["start_distribution_by_location"].size() < Model::get_config()->number_of_locations()
            ? 0
            : loc;
    add_distributions(ns["start_distribution_by_location"][input_loc],
                      result->start_distribution[loc]);
  }

  for (auto loc = 0; loc < Model::get_config()->number_of_locations(); loc++) {
    auto input_loc =
        ns["peak_distribution_by_location"].size() < Model::get_config()->number_of_locations()
            ? 0
            : loc;
    add_distributions(ns["peak_distribution_by_location"][input_loc],
                      result->peak_distribution[loc]);
  }

  add_therapies(ns, result.get());
  result->peak_after = ns["peak_after"].as<int>();
  return result;
}

std::unique_ptr<IStrategy>
StrategyBuilder::build_nested_mft_different_distribution_by_location_strategy(
    const YAML::Node &ns, const int &strategy_id) {
  auto result = std::make_unique<NestedMFTMultiLocationStrategy>();
  result->id = strategy_id;
  result->name = ns["name"].as<std::string>();

  result->distribution.clear();
  result->distribution.resize(static_cast<uint64_t>(Model::get_config()->number_of_locations()));

  result->start_distribution.clear();
  result->start_distribution.resize(
      static_cast<uint64_t>(Model::get_config()->number_of_locations()));

  result->peak_distribution.clear();
  result->peak_distribution.resize(
      static_cast<uint64_t>(Model::get_config()->number_of_locations()));

  for (auto loc = 0; loc < Model::get_config()->number_of_locations(); loc++) {
    auto input_loc =
        ns["start_distribution_by_location"].size() < Model::get_config()->number_of_locations()
            ? 0
            : loc;
    add_distributions(ns["start_distribution_by_location"][input_loc], result->distribution[loc]);
  }
  for (auto loc = 0; loc < Model::get_config()->number_of_locations(); loc++) {
    auto input_loc =
        ns["start_distribution_by_location"].size() < Model::get_config()->number_of_locations()
            ? 0
            : loc;
    add_distributions(ns["start_distribution_by_location"][input_loc],
                      result->start_distribution[loc]);
  }

  for (auto loc = 0; loc < Model::get_config()->number_of_locations(); loc++) {
    auto input_loc =
        ns["peak_distribution_by_location"].size() < Model::get_config()->number_of_locations()
            ? 0
            : loc;
    add_distributions(ns["peak_distribution_by_location"][input_loc],
                      result->peak_distribution[loc]);
  }

  for (auto i = 0; i < ns["strategy_ids"].size(); i++) {
    result->add_strategy(Model::get_strategy_db()[ns["strategy_ids"][i].as<int>()].get());
  }

  result->peak_after = ns["peak_after"].as<int>();
  //    std::cout << result->to_string() << std::endl;

  return result;
}

std::unique_ptr<IStrategy> StrategyBuilder::build_novel_drug_introduction_strategy(
    const YAML::Node &ns, const int strategy_id) {
  auto result = std::make_unique<NovelDrugIntroductionStrategy>();
  result->id = strategy_id;
  result->name = ns["name"].as<std::string>();

  add_distributions(ns["start_distribution"], result->start_distribution);
  add_distributions(ns["start_distribution"], result->distribution);
  add_distributions(ns["peak_distribution"], result->peak_distribution);

  result->peak_after = ns["peak_after"].as<int>();

  for (int i = 0; i < ns["strategy_ids"].size(); i++) {
    result->add_strategy(Model::get_strategy_db()[ns["strategy_ids"][i].as<int>()].get());
  }

  result->newly_introduced_strategy_id = ns["newly_introduced_strategy_id"].as<int>();
  result->tf_threshold = ns["tf_threshold"].as<double>();

  result->replacement_fraction = ns["replacement_fraction"].as<double>();
  result->replacement_duration = ns["replacement_duration"].as<int>();

  return result;
}

// TODO: this one need rework
// Read the YAML and create the DistrictMftStrategy object.
std::unique_ptr<IStrategy> StrategyBuilder::build_district_mft_strategy(const YAML::Node &node,
                                                                        const int &strategy_id) {
  spdlog::warn("This function need to be reworked. Use with caution.");
  // Load the initial parts of the object
  auto strategy = std::make_unique<DistrictMftStrategy>();
  strategy->id = strategy_id;
  strategy->name = node["name"].as<std::string>();

  // Track the districts that have been assigned an MFT, we should see all of
  // them once
  std::vector<int> districts;

  // Get district ID range from SpatialData
  const auto* spatial_data = Model::get_spatial_data();
  const auto min_district_id = spatial_data->get_boundary("district")->min_unit_id;
  const auto max_district_id = spatial_data->get_boundary("district")->max_unit_id;
  const auto expected_district_count = spatial_data->get_boundary("district")->unit_count;

  // Read each of the definitions
  for (auto ndx = 0; ndx < node["definitions"].size(); ndx++) {
    // Read the MFT from the child node
    auto child = node["definitions"][std::to_string(ndx)];
    // Use unique_ptr for exception safety
    auto template_mft = std::make_unique<DistrictMftStrategy::MftStrategy>();

    // Make sure the sizes are valid
    if (child["therapy_ids"].size() != child["distribution"].size()) {
      spdlog::error(
          "The The number of therapies and distributions should be the "
          "same, reading {}",
          ndx);
      throw std::invalid_argument("Matched therapy and distribution array size.");
    }

    // Read the therapy ids and make sure they make sense
    for (auto ndy = 0; ndy < child["therapy_ids"].size(); ndy++) {
      auto id = child["therapy_ids"][ndy].as<int>();
      if (id < 0) {
        spdlog::error("Drug id should not be less than zero, reading {}", ndx);
        throw std::invalid_argument("Drug id should not be less than zero.");
      }
      if (static_cast<std::size_t>(id) >= Model::get_therapy_db().size()) {
        spdlog::error("Drug id {} exceeds count of known drugs ({}), reading {}", id,
                      Model::get_therapy_db().size(), ndx);
        throw std::invalid_argument("Drug id exceeds count of known drugs.");
      }
      template_mft->therapies.push_back(id);
    }

    // Read the distribution percentages for the MFT and make sure they make
    // sense
    auto sum = 0.0;
    for (auto ndy = 0; ndy < child["distribution"].size(); ndy++) {
      auto percent = child["distribution"][ndy].as<double>();
      if (!std::isfinite(percent) || percent < 0.0) {
        spdlog::error("Distribution percentage must be finite and non-negative, reading {}", ndx);
        throw std::invalid_argument("Distribution percentage must be finite and non-negative.");
      }
      if (percent > 1.0) {
        spdlog::error("Distribution percentage cannot be greater than 100%, reading {}", ndx);
        throw std::invalid_argument("Distribution percentage cannot be greater than 100%.");
      }
      sum += percent;
      template_mft->percentages.push_back(percent);
    }
    // NOTE: the sum is compared against 1.0 with a tolerance. Truncating with
    // static_cast<int> rejected every valid distribution, since an accumulation
    // of decimal fractions in binary floating point essentially never lands
    // exactly on 1.0 (e.g. 0.765 + 0.085 + ... accumulates to 0.99999999...).
    constexpr double distribution_tolerance = 1e-6;
    if (std::fabs(sum - 1.0) > distribution_tolerance) {
      spdlog::error("Distribution percentage sum is {}, expected 1.0, reading {}", sum, ndx);
      throw std::invalid_argument("Distribution percentage sum must equal 100%.");
    }

    // Assign the MFT to each of the districts
    for (auto ndy = 0; ndy < child["district_ids"].size(); ndy++) {
      auto id = child["district_ids"][ndy].as<int>();

      // Validate district ID is within valid range
      if (id < min_district_id || id > max_district_id) {
        spdlog::error("Invalid district ID {}, valid range is {} to {}", id, min_district_id,
                      max_district_id);
        throw std::invalid_argument("District ID out of valid range");
      }

      if (std::ranges::find(districts, id) != districts.end()) {
        spdlog::error("District {} encountered a second time, reading {}", id, ndx);
        throw std::invalid_argument("District duplication detected.");
      }

      // Create a new copy of the MFT for this district
      auto district_mft = std::make_unique<DistrictMftStrategy::MftStrategy>(*template_mft);
      strategy->set_district_strategy(id, std::move(district_mft));
      districts.push_back(id);
    }

    // No need to explicitly delete template_mft, unique_ptr handles cleanup
  }

  // All the distributions have been read and assigned, make sure each district
  // has one
  if (districts.size() < expected_district_count) {
    spdlog::error(
        "Number of districts with MFT assigned ({}) is less than total "
        "district count ({})",
        districts.size(), expected_district_count);
    throw std::invalid_argument("Districts missing MFT assignment.");
  }

  // Everything looks good, return the strategy
  return strategy;
}

std::unique_ptr<IStrategy> StrategyBuilder::build_mft_age_based_strategy(const YAML::Node &node,
                                                                         const int &strategy_id) {
  auto result = std::make_unique<MFTAgeBasedStrategy>();
  result->id = strategy_id;
  result->name = node["name"].as<std::string>();
  for (std::size_t i = 0; i < node["therapy_ids"].size(); i++) {
    result->add_therapy(Model::get_therapy_db()[node["therapy_ids"][i].as<int>()].get());
  }
  add_distributions(node["age_boundaries"], result->age_boundaries);

  if (result->age_boundaries.size() != result->therapy_list.size() - 1) {
    spdlog::error(
        "The number of age boundaries should be one less than the "
        "number of therapies.");
    throw std::invalid_argument(
        "The number of age boundaries should be one less than the number of "
        "therapies.");
  }

  return result;
}

std::unique_ptr<IStrategy> StrategyBuilder::build_public_private_strategy(const YAML::Node &node,
                                                                          const int &strategy_id) {
  auto result = std::make_unique<PublicPrivateStrategy>();
  result->id = strategy_id;
  result->name = node["name"].as<std::string>();
  configure_child_strategies(node, strategy_id, result.get());
  result->start_public_share = read_public_share(node, "start_public_share");
  result->peak_public_share = read_public_share(node, "peak_public_share");
  result->public_share = result->start_public_share;
  result->peak_after = read_peak_after(node);
  return result;
}

std::unique_ptr<IStrategy> StrategyBuilder::build_public_private_multi_location_strategy(
    const YAML::Node &node, const int &strategy_id) {
  auto result = std::make_unique<PublicPrivateMultiLocationStrategy>();
  result->id = strategy_id;
  result->name = node["name"].as<std::string>();
  configure_child_strategies(node, strategy_id, result.get());

  const auto number_of_locations = Model::get_config()->number_of_locations();
  result->start_public_share_by_location =
      read_location_shares(node, "start_public_share_by_location", number_of_locations);
  result->peak_public_share_by_location =
      read_location_shares(node, "peak_public_share_by_location", number_of_locations);
  result->public_share_by_location = result->start_public_share_by_location;
  result->peak_after = read_peak_after(node);
  return result;
}
