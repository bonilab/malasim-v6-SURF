#ifndef POPULATION_H
#define POPULATION_H

#include <cstddef>
#include <list>
#include <memory>
#include <vector>

#include "Person/Person.h"

using PersonIndexPtrList = std::list<std::unique_ptr<PersonIndex>>;

class Model;
class PersonIndexAll;
class PersonIndexByLocationStateAgeClass;
class PersonIndexByLocationBitingLevel;
class Population {
public:
  Population(Population &&) = delete;
  Population &operator=(Population &&) = delete;
  // Disable copy and assignment
  Population(const Population &) = delete;
  Population &operator=(const Population &) = delete;

  Population();

  virtual ~Population();

  void initialize();
  //
  // void update(int current_time);
  //
  // // Add a person to the population
  void add_person(std::unique_ptr<Person> person);
  //
  // // Remove a person from the population
  void remove_person(Person* person);

  /**
   * This function removes person pointer out of all of the person indexes
   * This will also delete the @person out of memory
   * @param person
   */
  virtual void remove_dead_person(Person* person);

  // /** Return the total number of individuals in the simulation. */
  // virtual std::size_t size();

  /** Return the total number of individuals in the given location. */
  virtual std::size_t size_at(const int &location);

  /**
   * Return the number of individuals in the population
   * If the input location is -1, return total size
   * @param location
   */
  std::size_t size(const int &location = -1, const int &age_class = core::K_INVALID_AGE_CLASS);

  virtual std::size_t size(const int &location, const Person::HostStates &hs, const int &age_class);

  std::size_t size_residents_only(const int &location);

  /**
   * Notify change of a particular person's property to all person indexes
   * @param p
   * @param property
   * @param oldValue
   * @param newValue
   */
  virtual void notify_change(Person* person, const Person::Property &property,
                             const void* old_value, const void* new_value);

  virtual void perform_infection_event();

  void introduce_initial_cases();
  //
  void introduce_parasite(const int &location, Genotype* parasite_type,
                          const int &num_of_infections);

  static void setup_initial_infection(Person* person, Genotype* parasite_type);

  void persist_current_force_of_infection_to_use_n_days_later();

  void perform_birth_event();

  void perform_death_event();

  void generate_individual(int location, int age_class);

  void give_1_birth(const int &location);

  void clear_all_dead_state_individual();

  void perform_circulation_event();

  void perform_circulation_for_1_location(const int &from_location, const int &target_location,
                                          const int &number_of_circulations,
                                          std::vector<Person*> &today_circulations);

  bool has_0_case();

  void initialize_person_indices();

  void update_all_individuals();

  void execute_all_individual_events(int up_to_time);

  void update_current_foi();

  // Notify the population that a person has moved from the source location, to
  // the destination location
  void notify_movement(int source, int destination);

  PersonIndexPtrList* person_index_list() { return person_index_list_.get(); }
  PersonIndexAll* all_persons() { return all_persons_.get(); }

  template <typename T>
  T* get_person_index();

  IntVector get_popsize_by_location() { return popsize_by_location_; }
  void set_popsize_by_location(const IntVector &popsize_by_location) {
    popsize_by_location_ = popsize_by_location;
  }

  [[nodiscard]] std::vector<std::vector<double>> &individual_foi_by_location() {
    return individual_foi_by_location_;
  }

  [[nodiscard]] std::vector<std::vector<double>> &individual_relative_biting_by_location() {
    return individual_relative_biting_by_location_;
  }

  [[nodiscard]] std::vector<std::vector<double>> &individual_relative_moving_by_location() {
    return individual_relative_moving_by_location_;
  }

  [[nodiscard]] std::vector<double> &sum_relative_biting_by_location() {
    return sum_relative_biting_by_location_;
  }

  [[nodiscard]] std::vector<double> &sum_relative_moving_by_location() {
    return sum_relative_moving_by_location_;
  }

  [[nodiscard]] std::vector<double> &current_force_of_infection_by_location() {
    return current_force_of_infection_by_location_;
  }

  [[nodiscard]] std::vector<std::vector<double>> &force_of_infection_for_n_days_by_location() {
    return force_of_infection_for_n_days_by_location_;
  }

  [[nodiscard]] std::vector<std::vector<Person*>> &all_alive_persons_by_location() {
    return all_alive_persons_by_location_;
  }

private:
  std::unique_ptr<PersonIndexAll> all_persons_{nullptr};

  std::unique_ptr<PersonIndexPtrList> person_index_list_{nullptr};
  IntVector popsize_by_location_;

  std::vector<std::vector<double>> individual_foi_by_location_;
  std::vector<std::vector<double>> individual_relative_biting_by_location_;
  std::vector<std::vector<double>> individual_relative_moving_by_location_;

  std::vector<double> sum_relative_biting_by_location_;
  std::vector<double> sum_relative_moving_by_location_;

  std::vector<double> current_force_of_infection_by_location_;
  std::vector<std::vector<double>> force_of_infection_for_n_days_by_location_;
  std::vector<std::vector<Person*>> all_alive_persons_by_location_;
};

template <typename T>
T* Population::get_person_index() {
  for (auto &person_index_ptr : *person_index_list_) {
    T* pi = dynamic_cast<T*>(person_index_ptr.get());
    if (pi != nullptr) { return pi; }
  }
  return nullptr;
}
#endif  // POPULATION_H

