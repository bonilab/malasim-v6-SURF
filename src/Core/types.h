#ifndef MALASIM_CORE_TYPES_H
#define MALASIM_CORE_TYPES_H

#include <cstdint>
#include <limits>

namespace core {

// ---------- Basic scalar types ----------

using Age = std::uint8_t;
using AgeClass = std::uint8_t;
using MovingLevel = std::uint8_t;

// max number of people can be 4,294,967,295, which should be sufficient for our simulation
// but cannot be fit into 16 bits, so we use 32 bits for PersonId
using PersonId = std::uint32_t;
// max number of locations can be 4,294,967,295, which should be sufficient for our simulation
// but cannot be fit into 16 bits (32767), so we use 32 bits for LocationId
// 4000 km * 4000 km <-> 800 x 800 grid with 5 km x 5 km cell size = 640,000 location
using LocationId = std::uint32_t;

// could be fit to 16 bits, but we use 32 bits for GenotypeId to be safe
// need to re-check if we can fit all genotypes into 16 bits after we finalize the genotype
// representation
using GenotypeId = std::uint32_t;

// we have small number of drugs and therapies < 256
// uint8_t is enough
using TherapyId = std::uint8_t;
using DrugId = std::uint8_t;

// the max time can be 32767 days, which is about 89 years, which should be sufficient for our
// simulation. Day can be negative.
using SimDay = std::int16_t;
using EventTime = std::int16_t;

// if we simulate for 40 years, and each person can be bitten at most 1 time per day
// then the max number of bites can be 40 * 365 = 14600, which can fit into 16 bits
// the same for trip count
using BiteCount = std::uint16_t;
using TripCount = std::uint16_t;

// ---------- Invalid / sentinel values ----------

inline constexpr Age K_INVALID_AGE = std::numeric_limits<Age>::max();
inline constexpr AgeClass K_INVALID_AGE_CLASS = std::numeric_limits<AgeClass>::max();
inline constexpr MovingLevel K_INVALID_MOVING_LEVEL = std::numeric_limits<MovingLevel>::max();

inline constexpr PersonId K_INVALID_PERSON_ID = std::numeric_limits<PersonId>::max();
inline constexpr LocationId K_INVALID_LOCATION_ID = std::numeric_limits<LocationId>::max();
inline constexpr GenotypeId K_INVALID_GENOTYPE_ID = std::numeric_limits<GenotypeId>::max();
inline constexpr TherapyId K_INVALID_THERAPY_ID = std::numeric_limits<TherapyId>::max();
inline constexpr DrugId K_INVALID_DRUG_ID = std::numeric_limits<DrugId>::max();

inline constexpr SimDay K_INVALID_SIM_DAY = -1;
inline constexpr EventTime K_INVALID_EVENT_TIME = -1;

// ---------- Domain limits ----------

inline constexpr Age K_MAX_HUMAN_AGE = 120;

}  // namespace core
//
#endif  // MALASIM_CORE_TYPES_H
