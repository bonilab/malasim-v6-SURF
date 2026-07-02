#ifndef SINGLEHOSTCLONALPARASITEPOPULATIONS_H
#define SINGLEHOSTCLONALPARASITEPOPULATIONS_H

#include <memory>
#include <vector>

#include "Population/ClonalParasitePopulation.h"
#include "Utils/TypeDef.h"

class ClonalParasitePopulation;
class Person;
class DrugType;
class DrugsInBlood;
class ParasiteDensityUpdateFunction;

class SingleHostClonalParasitePopulations {
  // OBJECTPOOL(SingleHostClonalParasitePopulations)
public:
  // Disallow copy
  SingleHostClonalParasitePopulations(const SingleHostClonalParasitePopulations&) = delete;
  SingleHostClonalParasitePopulations& operator=(const SingleHostClonalParasitePopulations&) = delete;

  // Disallow move
  SingleHostClonalParasitePopulations(SingleHostClonalParasitePopulations&&) = delete;
  SingleHostClonalParasitePopulations& operator=(SingleHostClonalParasitePopulations&&) = delete;

  // Use constexpr for compile-time constant
  static constexpr double DEFAULT_LOG_DENSITY = ClonalParasitePopulation::LOG_ZERO_PARASITE_DENSITY;

  // Use explicit for single-parameter constructors
  explicit SingleHostClonalParasitePopulations(Person* person = nullptr);

  // Mark destructor as default if it doesn't need special handling
  virtual ~SingleHostClonalParasitePopulations();

  void init();

  // Iterator type definitions for STL compatibility
  using Iterator = std::vector<std::unique_ptr<ClonalParasitePopulation>>::iterator;
  using ConstIterator = std::vector<std::unique_ptr<ClonalParasitePopulation>>::const_iterator;

  // Iterator methods
  [[nodiscard]] ConstIterator begin() const noexcept { return parasites_.begin(); }
  [[nodiscard]] ConstIterator end() const noexcept { return parasites_.end(); }
  [[nodiscard]] Iterator begin() noexcept { return parasites_.begin(); }
  [[nodiscard]] Iterator end() noexcept { return parasites_.end(); }

  // Access methods
  [[nodiscard]] ClonalParasitePopulation* at(size_t index) const {
    return parasites_.at(index).get();
  }

  [[nodiscard]] ClonalParasitePopulation* operator[](size_t index) const {
    return parasites_[index].get();
  }

  // Mark virtual functions that are meant to be overridden with override in derived classes
  [[nodiscard]] virtual size_t size() const noexcept { return parasites_.size(); }

  [[nodiscard]] bool empty() const noexcept { return parasites_.empty(); }

  virtual void add(std::unique_ptr<ClonalParasitePopulation> blood_parasite);

  virtual void remove(size_t index);

  [[nodiscard]] virtual int latest_update_time() const;

  virtual bool contain(ClonalParasitePopulation* blood_parasite);

  void change_all_parasite_update_function(ParasiteDensityUpdateFunction* from,
                                           ParasiteDensityUpdateFunction* to) const;

  void update();

  void clear_cured_parasites(double cured_threshold);

  void clear();

  void update_by_drugs(DrugsInBlood* drugs_in_blood) const;

  void apply_cnv_reversion(DrugsInBlood* drugs_in_blood) const;

  [[nodiscard]] bool has_detectable_parasite(double detectable_threshold) const;

  [[nodiscard]] bool is_gametocytaemic() const;

  // Use [[nodiscard]] consistently and reference member functions that don't modify state as const
  [[nodiscard]] double log10_total_infectious_density() const noexcept {
    return log10_total_infectious_density_;
  }

  void set_log10_total_infectious_density(double value) noexcept {  // Pass by value for primitives
    log10_total_infectious_density_ = value;
  }

  [[nodiscard]] Person* person() const noexcept { return person_; }

  void set_person(Person* value) noexcept { person_ = value; }

private:
  Person* person_{nullptr};
  std::vector<std::unique_ptr<ClonalParasitePopulation>> parasites_;
  double log10_total_infectious_density_{DEFAULT_LOG_DENSITY};
};

#endif /* SINGLEHOSTCLONALPARASITEPOPULATIONS_H */
