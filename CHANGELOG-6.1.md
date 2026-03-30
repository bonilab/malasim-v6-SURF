# Changelog v6.1

## Bug Fixes

- **Population out-of-bounds access**: Added guards for `K_INVALID_LOCATION_ID` in `Population::add_person()` and `remove_person()` to prevent crashes when Person is created with default constructor (location_ = K_INVALID_LOCATION_ID)
- **MFTAgeBasedStrategyTest**: Reduced simulation time from 100 to 80 years to fit in `core::SimDay` (int16_t max ~89 years)
- **RouletteTest test pollution**: Added `Model::get_instance()->release()` in SetUp/TearDown to properly clean up singleton state between tests; changed from `set_location(i)` to `set_number_of_times_bitten(i)` as temporary identifier to avoid side effects

## Refactor

- **Age and Age Class Types**: Changed `age_` and `age_class_` member variables in `Person` class from `int` to `uint8_t` with symbolic sentinel values
  - Added `Age = uint8_t` and `AgeClass = uint8_t` typedefs in `Person.h`
  - Added `K_INVALID_AGE` (255) and `K_INVALID_AGE_CLASS` (255) constants
  - Updated `Person::age_` and `Person::age_class_` to use new types with appropriate sentinel values

- **Sentinel Comparisons**: Fixed comparisons using `-1` sentinel to use symbolic constants
  - `Person.cpp`: `age_class_ == -1` → `age_class_ == K_INVALID_AGE_CLASS`
  - `Population.h/cpp`: `age_class == -1` → `age_class == K_INVALID_AGE_CLASS`

- **Tests Updated**: `PersonBasicTest.cpp` updated to use `K_INVALID_AGE` constant

- **Test logging**: Changed default log level from `info` to `warn`; set `MALASIM_LOG_LEVEL` env var to control logging (e.g., `MALASIM_LOG_LEVEL=info` for verbose output)

## Files Changed

- `src/Population/Person/Person.h` - Added type aliases and constants
- `src/Population/Person/Person.cpp` - Fixed sentinel comparison
- `src/Population/Population.h` - Updated default parameter
- `src/Population/Population.cpp` - Fixed sentinel comparison
- `tests/Population/Person/PersonBasicTest.cpp` - Updated test assertion
- `tests/Core/Random/RandomTest_random_roulette.cpp` - Fixed test pollution and identifier overflow
- `tests/Treatment/Strategies/MFTAgeBasedStrategyTest.cpp` - Fixed simulation time overflow
- `tests/SpdlogEnvironment.cpp` - Added MALASIM_LOG_LEVEL env var support
