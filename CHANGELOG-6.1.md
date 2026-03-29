# Changelog v6.1

## Refactor

- **Age and Age Class Types**: Changed `age_` and `age_class_` member variables in `Person` class from `int` to `uint8_t` with symbolic sentinel values
  - Added `Age = uint8_t` and `AgeClass = uint8_t` typedefs in `Person.h`
  - Added `K_INVALID_AGE` (255) and `K_INVALID_AGE_CLASS` (255) constants
  - Updated `Person::age_` and `Person::age_class_` to use new types with appropriate sentinel values

- **Sentinel Comparisons**: Fixed comparisons using `-1` sentinel to use symbolic constants
  - `Person.cpp`: `age_class_ == -1` → `age_class_ == K_INVALID_AGE_CLASS`
  - `Population.h/cpp`: `age_class == -1` → `age_class == K_INVALID_AGE_CLASS`

- **Tests Updated**: `PersonBasicTest.cpp` updated to use `K_INVALID_AGE` constant

## Files Changed

- `src/Population/Person/Person.h` - Added type aliases and constants
- `src/Population/Person/Person.cpp` - Fixed sentinel comparison
- `src/Population/Population.h` - Updated default parameter
- `src/Population/Population.cpp` - Fixed sentinel comparison
- `tests/Population/Person/PersonBasicTest.cpp` - Updated test assertion
