# Events

This module provides the core event system and individual event implementations for the simulation. These events handle various aspects of disease progression, treatment, and individual behavior in the simulation.

## Directory Structure

### Core Components
- `Event.h/cpp`: Base class for all events in the system
  - Provides common functionality for event scheduling and execution
  - Defines the interface that all events must implement

### Subdirectories
- `Population/`: Population-level events (e.g., MDA, mutation events)
- `Environment/`: Environmental events
- `Trials/`: Clinical trial related events

## Event Categories

### Disease Progression Events
- `ProgressToClinicalEvent.h/cpp`: Handles progression to clinical symptoms
- `MatureGametocyteEvent.h/cpp`: Manages gametocyte maturation process
- `MoveParasiteToBloodEvent.h/cpp`: Controls parasite movement to bloodstream
- `EndClinicalEvent.h/cpp`: Manages the end of clinical symptoms

### Treatment Events
- `ReceiveTherapyEvent.h/cpp`: Handles individual treatment administration
- `ReceiveMDATherapyEvent.h/cpp`: Manages Mass Drug Administration treatment
- `TestTreatmentFailureEvent.h/cpp`: Tests for treatment failure
- `ReportTreatmentFailureDeathEvent.h/cpp`: Reports deaths due to treatment failure
- `UpdateWhenDrugIsPresentEvent.h/cpp`: Drug presence monitoring

### Movement Events
- `CirculateToTargetLocationNextDayEvent.h/cpp`: Manages individual movement between locations
- `ReturnToResidenceEvent.h/cpp`: Handles return to residence location

### Lifecycle Events
- `BirthdayEvent.h/cpp`: Manages individual aging
- `SwitchImmuneSystemModeEvent.h/cpp`: Switches from infant to non-infant immunity

### Monitoring Events
- `RaptEvent.h/cpp`: Rapid Assessment of Parasite Treatment event

## Implementation Guidelines

1. **Event Class Structure**
   - All events inherit from `WorldEvent` or `PersonEvent` base classes
   - Follow the standardized format for disallowing copy and move operations
   - Use `[[nodiscard]]` for getters
   - Make member variables private with appropriate getters/setters

2. **Documentation**
   - Use Doxygen-style comments for public interfaces
   - Document parameters and return values
   - Include brief descriptions of event purpose

3. **Memory Management**
   - Use RAII principles
   - Properly handle resource cleanup
   - Use smart pointers where appropriate

4. **Error Handling**
   - Validate input parameters
   - Use appropriate error reporting mechanisms
   - Handle edge cases gracefully

## Dependencies

- Core simulation components:
  - `Model`
  - `Person`
  - `Population`
  - `DrugDatabase`
- Utility classes:
  - `Config`
  - `Random`
  - `Scheduler`

## Usage Example

```cpp
// Creating and scheduling a therapy event
auto therapy_event = new ReceiveTherapyEvent(person, therapy_id);
scheduler->schedule_event(therapy_event);

// Creating a clinical progression event
auto clinical_event = new ProgressToClinicalEvent(person);
scheduler->schedule_event(clinical_event);
```

## Notes

- Events are executed in chronological order
- Events may create other events during execution
- All events must be memory managed by the scheduler
- Follow C++ guidelines for class design and implementation
- Use const-correctness throughout
- Maintain thread safety where applicable
