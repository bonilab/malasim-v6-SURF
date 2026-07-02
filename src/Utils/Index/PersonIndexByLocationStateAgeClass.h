#ifndef PERSONINDEXBYLOCATIONSTATEAGECLASS_H
#define    PERSONINDEXBYLOCATIONSTATEAGECLASS_H

#include "Utils/TypeDef.h"
#include "Population/Person/Person.h"
#include "PersonIndex.h"

class PersonIndexByLocationStateAgeClass : public PersonIndex {
public:
  //disable copy and assign
  PersonIndexByLocationStateAgeClass(const PersonIndexByLocationStateAgeClass &) = delete;
  void operator=(const PersonIndexByLocationStateAgeClass &) = delete;

  PersonPtrVector4 &vPerson() {
    return vPerson_;
  }
  void vPerson(const PersonPtrVector4 &v) {
    vPerson_ = v;
  }

private:
  PersonPtrVector4 vPerson_;

 public:
  //    PersonIndexByLocationStateAgeClass();
  PersonIndexByLocationStateAgeClass(const int &no_location = 1, const int &no_host_state = 1,
                                     const int &no_age_class = 1);

  //    PersonIndexByLocationStateAgeClass(const PersonIndexByLocationStateAgeClass& orig);
  virtual ~PersonIndexByLocationStateAgeClass();

  void Initialize(const int &no_location = 1, const int &no_host_state = 1, const int &no_age_class = 1);

  virtual void add(Person *p);

  virtual void remove(Person *p);

  virtual std::size_t size() const;

  virtual void update();

  virtual void notify_change(Person *p, const Person::Property &property, const void *oldValue, const void *newValue);

 private:
  void remove_without_set_index(Person *p);

  void add(Person *p, core::LocationId location, const Person::HostStates &host_state, core::AgeClass age_class);

  void change_property(Person *p, core::LocationId location, const Person::HostStates &host_state, core::AgeClass age_class);
};

#endif    /* PERSONINDEXBYLOCATIONSTATEAGECLASS_H */

