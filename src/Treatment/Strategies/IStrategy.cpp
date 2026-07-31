#include "IStrategy.h"

std::map<std::string, IStrategy::StrategyType> IStrategy::strategy_type_map{
    {"SFT", StrategyType::SFT},
    {"Cycling", StrategyType::Cycling},
    {"AdaptiveCycling", StrategyType::AdaptiveCycling},
    {"MFT", StrategyType::MFT},
    {"MFTRebalancing", StrategyType::MFTRebalancing},
    {"NestedMFT", StrategyType::NestedMFT},
    {"MFTMultiLocation", StrategyType::MFTMultiLocation},
    {"NestedMFTMultiLocation", StrategyType::NestedMFTMultiLocation},
    {"NovelDrugIntroduction", StrategyType::NovelDrugIntroduction},
    {"DistrictMFT", StrategyType::DistrictMft},
    {"MFTAgeBased", StrategyType::MFTAgeBased},
    {"PublicPrivate", StrategyType::PublicPrivate},
    {"PublicPrivateMultiLocation", StrategyType::PublicPrivateMultiLocation},
};
