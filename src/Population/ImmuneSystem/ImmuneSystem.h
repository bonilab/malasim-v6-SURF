#ifndef IMMUNESYSTEM_H
#define IMMUNESYSTEM_H

#include <cstdint>

#include "Core/types.h"

class Model;

class Person;

enum class ImmuneSystemMode : std::uint8_t {
  INFANT,
  NON_INFANT,
};

class ImmuneSystem {
  // OBJECTPOOL(ImmuneSystem)
public:
  // Disallow copy
  ImmuneSystem(const ImmuneSystem &) = delete;
  ImmuneSystem &operator=(const ImmuneSystem &) = delete;
  ImmuneSystem(ImmuneSystem &&) = delete;
  ImmuneSystem &operator=(ImmuneSystem &&) = delete;

  explicit ImmuneSystem(Person* person = nullptr);

  virtual ~ImmuneSystem();

  [[nodiscard]] Person* person() const { return person_; }
  void set_person(Person* person) { person_ = person; }

  [[nodiscard]] ImmuneSystemMode mode() const { return mode_; }
  void initialize_as_infant() { mode_ = ImmuneSystemMode::INFANT; }
  void switch_to_non_infant();

  [[nodiscard]] bool increase() const { return increase_; }
  void set_increase(bool increase) { increase_ = increase; }

  virtual void draw_random_immune();

  virtual void update();

  [[nodiscard]] virtual double get_latest_immune_value() const;

  virtual void set_latest_immune_value(double value);

  [[nodiscard]] virtual double get_current_value() const;

  [[nodiscard]] virtual double get_parasite_size_after_t_days(const int &duration,
                                                              const double &original_size,
                                                              const double &fitness) const;

  [[nodiscard]] virtual double get_clinical_progression_probability() const;

private:
  [[nodiscard]] double get_decay_rate(core::Age age) const;
  [[nodiscard]] double get_acquire_rate(core::Age age) const;
  [[nodiscard]] double get_one_day_decay_factor() const;
  [[nodiscard]] double get_one_day_acquire_factor(core::Age age) const;

  Person* person_{nullptr};
  double latest_value_{0.0};
  ImmuneSystemMode mode_{ImmuneSystemMode::NON_INFANT};
  bool increase_{false};
};

#endif /* IMMUNESYSTEM_H */
