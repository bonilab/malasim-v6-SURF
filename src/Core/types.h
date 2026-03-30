#ifndef MALASIM_CORE_TYPES_H
#define MALASIM_CORE_TYPES_H

#include <cstdint>
#include <limits>

namespace core {

// ---------- Basic scalar types ----------

using Age = std::uint8_t;
using AgeClass = std::uint8_t;
using MovingLevel = std::uint8_t;

using PersonId = std::uint32_t;
using LocationId = std::uint32_t;
using GenotypeId = std::uint32_t;
using TherapyId = std::uint32_t;
using DrugId = std::uint16_t;

using Day = std::int32_t;
using EventTime = std::int32_t;

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

inline constexpr Day K_INVALID_DAY = -1;
inline constexpr EventTime K_INVALID_EVENT_TIME = -1;

// ---------- Domain limits ----------

inline constexpr Age kMaxHumanAge = 120;

}  // namespace core
//
#endif  // MALASIM_CORE_TYPES_H
