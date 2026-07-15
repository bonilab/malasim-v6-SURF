# ImmuneSystem

This module implements the immune system modeling for malaria infection, including age-dependent immunity development and response to parasites.

## Overview

The ImmuneSystem module manages:
- Age-specific immune responses
- Immunity acquisition and decay
- Clinical protection
- Parasite density control
- Treatment interaction effects

## Components

### Base System
- `ImmuneSystem`: Immunity state, mode transitions, and immune response calculations
- `ImmunityClearanceUpdateFunction`: Immunity-based parasite clearance

### Age-Specific Modes
- `ImmuneSystemMode::Infant`: Maternal immunity decay until six months of age
- `ImmuneSystemMode::NonInfant`: Configured acquisition and decay after six months

The six-month event changes the immune-system mode in place after the daily person update,
preserving the current immune value.

## Key Features

### Immunity Types
- Clinical immunity
- Parasite immunity
- Transmission immunity
- Maternal protection
- Drug-enhanced clearance

### Dynamic Responses
- Exposure-based development
- Age-dependent changes
- Decay over time
- Boost from infections
- Treatment interactions

## Implementation

### Immune System Structure
```cpp
class ImmuneSystem {
    Person* person_;
    double latest_value_;
    ImmuneSystemMode mode_;
    bool increase_;
};
```

### Update Functions
- Immunity level calculations
- Clinical protection updates
- Parasite density effects
- Drug interaction modifiers
- Clearance probability

## Usage

```cpp
// Initialize immune system
auto immune_system = new ImmuneSystem(age);

// Update immunity levels
immune_system->update();

// Calculate immune effects
double clinical_immunity = immune_system->get_clinical_immunity();
double parasite_immunity = immune_system->get_parasite_immunity();

// Handle infection interaction
immune_system->handle_infection(parasite_density);
```

## Key Algorithms

### Immunity Development
- Age-based progression
- Exposure-driven enhancement
- Maternal protection decay
- Clinical immunity acquisition
- Parasite immunity buildup

### Clearance Mechanisms
- Base clearance rates
- Immunity modification
- Drug interaction effects
- Age-specific factors
- Exposure history impact

## Dependencies

- Core components:
  - `Person`
  - `Config`
  - `Model`
- Disease components:
  - `ClonalParasitePopulation`
  - `DrugsInBlood`

## Notes

- Age significantly affects immunity
- Exposure history is crucial
- Multiple immunity types interact
- Drug effects are modified by immunity
- Memory management is important
- Update frequency affects accuracy
- Component selection is age-based
