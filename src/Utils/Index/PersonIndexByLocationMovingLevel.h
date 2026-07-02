#ifndef PERSONINDEXBYLOCATIONMOVINGLEVEL_H
#define    PERSONINDEXBYLOCATIONMOVINGLEVEL_H

#include "Utils/TypeDef.h"
#include "Population/Person/Person.h"
#include "PersonIndex.h"

class PersonIndexByLocationMovingLevel : public PersonIndex {
public:
  //disable copy and assign
  PersonIndexByLocationMovingLevel(const PersonIndexByLocationMovingLevel &) = delete;
    void operator=(const PersonIndexByLocationMovingLevel &) = delete;

public:
  PersonPtrVector3 &vPerson() {
    return vPerson_;
  }
  void vPerson(const PersonPtrVector3 &v) {
    vPerson_ = v;
  }

private:
  PersonPtrVector3 vPerson_;
 public:
  PersonIndexByLocationMovingLevel(const int &no_location = 1, const int &no_level = 1);

  //    PersonIndexByLocationMovingLevel(const PersonIndexByLocationMovingLevel& orig);
  virtual ~PersonIndexByLocationMovingLevel();

  void Initialize(const int &no_location = 1, const int &no_level = 1);

  virtual void add(Person *p);

  virtual void remove(Person *p);

  virtual std::size_t size() const;

  virtual void update();

  virtual void notify_change(Person *p, const Person::Property &property, const void *oldValue, const void *newValue);

 private:
  void remove_without_set_index(Person *p);

  void add(Person *p, core::LocationId location, core::MovingLevel moving_level);

  void change_property(Person *p, core::LocationId location, core::MovingLevel moving_level);

};

#endif    /* PERSONINDEXBYLOCATIONMOVINGLEVEL_H */

