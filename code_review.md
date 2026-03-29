Code Review: MalaSim-v6-SURF
Executive Summary
A well-structured C++ malaria simulation codebase with generally sound architectural decisions, modern memory management, and good code organization. However, it has several systemic issues — particularly in error handling, ownership clarity, and the overuse of global state — that should be addressed.

1. Architecture & Design
   Strengths
   Clear layering: Core (Scheduler/EventManager), Population, Events, Configuration, Treatment, Reporters are well-separated.
   Event-driven architecture using a templated EventManager.
   Strategy pattern for treatment via IStrategy/Therapy hierarchy.
   Heavy use of std::unique_ptr throughout.
   Issues
   1.1 Global Singleton Overuse (src/Simulation/Model.h:29-31)
   Model::get_instance() is accessed globally from nearly every subsystem. This creates tight coupling, prevents unit testing in isolation, and introduces hidden dependencies. Components like Scheduler, Population, and EventManager all depend on the singleton. Recommendation: Use dependency injection or interface-based design to decouple components.

1.2 God Class: Person (src/Population/Person/Person.h:25-27)
Person inherits from PersonIndexAllHandler, PersonIndexByLocationStateAgeClassHandler, PersonIndexByLocationMovingLevelHandler simultaneously. This violates the Single Responsibility Principle and makes changes to any indexing strategy require modifying Person. Recommendation: Decouple indexing from the entity using observer/visitor patterns.

1.3 RTTI Abuse in EventManager (src/Core/Scheduler/EventManager.h:68,77,92)
has_event<T>() iterates all events and uses dynamic_cast to find events of a specific type — O(n) with expensive casts. Recommendation: Maintain a type-indexed map of event types for O(1) lookups, or use a visitor pattern.

2. Code Quality
   2.1 Raw Pointer Ownership Ambiguity (src/Population/ClonalParasitePopulation.h:69-71, Person.h:308-329)
   Numerous raw pointer data members (Genotype* genotype\_, Population* population*, liver_parasite_type*) lack ownership documentation. It is unclear if they are owned by the class, a registry, or the caller. Recommendation: Use gsl::observer_ptr<T> for non-owning pointers, or add explicit ownership comments.

2.2 Inconsistent Null Pointer Checking
Some locations defensively check for nullptr (e.g., BirthdayEvent.cpp:15-18) while others access pointers without checking (Person.cpp:60). This inconsistency risks crashes at edge cases.

2.3 28+ Incomplete TODOs (src/Population/Person/Person.cpp:84-95 and others)
There are unresolved TODOs including // TODO: remove all events and // TODO: if age exceeds limit, remove person??? — indicating incomplete implementation in core paths.

2.4 Debug Output Left in Production Code (src/Utils/ObjectPool.h:142)

C++
std::cout << "Allocating new chunk..." << '\n'; // Should be spdlog or removed
Mix of std::cout and spdlog — the former should be removed or replaced consistently.

2.5 String Concatenation in Log Calls (src/Simulation/Model.cpp:43,51)

C++
spdlog::info("Loading configuration file: " + utils::Cli::get_instance().get_input_path());
Creates temporary strings unnecessarily. Use fmt-style: spdlog::info("Loading configuration file: {}", path).

2.6 Repetitive Type Aliases (src/Utils/TypeDef.h:26-50)
DoubleVector, DoubleVector2, DoubleVector3, IntVector, IntVector2... etc. defined manually. Consider template aliases or removing unused ones.

3. Memory Management
   3.1 reset(raw_ptr) in Setters (src/Simulation/Model.h:128)

C++
static void set*mdc(ModelDataCollector\* mdc) { get_instance()->mdc*.reset(mdc); }
Accepting a raw pointer in a setter and calling reset() is dangerous — the caller retains a dangling pointer. Fix: Accept std::unique_ptr<ModelDataCollector> and std::move it in.

3.2 ObjectPool Placement New Not Fully Safe (src/Utils/ObjectPool.h:67)
Placement new is used correctly, but there's no guard against callers accidentally calling regular delete on pooled objects. Add static_assert or prominent documentation.

3.3 Dangling Pointer Risk in Population Destructor (src/Population/Population.cpp:34-46)
Population destructs the PersonUniquePtrVector (clearing all owned Person objects), but other subsystems may hold raw Person\* pointers. If destruction order is wrong, accessing those pointers is undefined behavior.

4. Error Handling
   4.1 No Exception Handling Around Core Simulation Loop (src/Core/Scheduler/Scheduler.cpp:21-33)
   Scheduler::run() calls daily_update() in a tight loop with no try-catch. Any uncaught exception (from event execution, configuration access, etc.) will crash the simulation with no meaningful diagnostic.

4.2 Ignored Return Values
population\_->initialize() and similar calls in Model.cpp:92 return void and provide no indication of failure. Failures silently proceed.

4.3 Asserts in Production Logic (src/Utils/Index/PersonIndexByLocationStateAgeClass.cpp:33-36)

C++
assert(vPerson\_.size() > p->get_location() && p->get_location() >= 0);
These are development-only guards. In a release build (with NDEBUG), they vanish. Replace with runtime checks and proper error messages.

4.4 Configuration Validation Unused (src/Configuration/Config.h:45-46)
validate_all_cross_field_validations() is declared but unclear if it's actually called at load time. Configuration should be fully validated at startup.

5. Thread Safety
   5.1 Model Singleton (src/Simulation/Model.h:29-32)
   In C++11, static local variable initialization is thread-safe. However, all access to the singleton's mutable state is unprotected. This is fine while single-threaded, but would need locks if parallelism is ever introduced.

5.2 ObjectPool Thread-Safety Incomplete (src/Utils/ObjectPool.h:48-77)
The lock in the thread-safe variant protects the free list, but placement new (object construction) happens outside the lock. Document that construction must be thread-safe independently.

6. Performance
   6.1 O(n) has_event<T>() with Expensive dynamic_cast (src/Core/Scheduler/EventManager.h:65-71)
   Called frequently. Should use type-indexed storage for O(1) lookup.

6.2 Parallel Vectors Causing Cache Misses (src/Population/Population.h:125-155)
individual*foi_by_location*, individual*relative_biting_by_location*, individual*relative_moving_by_location* are all separate vectors over the same individuals. Consider a Struct-of-Arrays or unified struct.

6.3 No Person Update Batching (src/Population/Population.cpp:104)
All persons are updated daily with no batching, vectorization, or parallelism — potentially a bottleneck at large population sizes.

7. C++ Best Practices
   7.1 Missing [[nodiscard]]
   bool load(...), bool initialize(), std::size_t size(...) should all be [[nodiscard]] to prevent callers from silently ignoring return values.

7.2 Missing noexcept Specifications
Critical paths like Scheduler::run() should document exception behavior with noexcept or have explicit exception handling.

7.3 Inconsistent Modern C++ Usage
Structured bindings and constexpr are used in some files (e.g., EventManager.h:67, ObjectPool.h:114) but older patterns appear elsewhere. No use of std::optional or std::variant where they would improve clarity.

7.4 RAII Not Fully Embraced
Model.cpp has a manual release() function. True RAII would handle cleanup automatically in destructors.

8. Security
   8.1 File Path Not Validated (src/Simulation/Model.cpp:43-44)
   The configuration file path from CLI is passed directly without canonicalization or sandbox validation. A path like ../../etc/passwd would be accepted.

8.2 No Bounds Validation on YAML Configuration
Numeric parameters (ages, location indices, probabilities) are read from YAML without range validation. Out-of-range values can cause buffer overruns, assertion failures, or incorrect simulation results.

8.3 YAML Error Handling Too Generic (src/Events/Population/EventBuilders/IntroduceMutantEventBuilder.cpp:71-123)
Caught exceptions don't identify which field failed parsing. This makes diagnostics very difficult for misconfigured inputs.

9. Recommendations (Priority Order)
   Priority Action
   🔴 Critical Wrap Scheduler::run() / daily_update() in try-catch with diagnostic logging
   🔴 Critical Add range validation for all YAML configuration parameters
   🔴 Critical Document or fix raw pointer ownership semantics
   🟠 High Fix setter reset(raw_ptr) → accept unique_ptr
   🟠 High Replace assert() in production logic with runtime checks
   🟠 High Resolve or remove 28+ TODO items
   🟠 High Replace std::cout with spdlog consistently
   🟡 Medium Optimize EventManager::has_event<T>() to avoid O(n) scan
   🟡 Medium Consolidate parallel vectors in Population
   🟡 Medium Add [[nodiscard]] to important return values
   🟢 Low Refactor toward dependency injection to reduce global singleton reliance
   🟢 Low Consider std::optional/std::variant for optional/discriminated values
   Summary Statistics
   Metric Assessment
   Source files ~289
   Smart pointer coverage ~95% of owned objects
   Error handling coverage ~5% (2–3 try-catch blocks)
   TODO items 28+ unresolved
   Raw pointer clarity Moderate (non-owning use, but undocumented)
   Test coverage Present but limited scope
   Overall Assessment: C+/B− — Good foundational architecture and modern memory management, but robustness, error handling, and testability need significant work before production deployment. The most urgent concerns are lack of error handling around the core simulation loop and insufficient input validation from YAML configuration.
