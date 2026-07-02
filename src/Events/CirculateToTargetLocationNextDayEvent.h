/*
 * CirculateToTargetLocationNextDayEvent.h
 *
 * Define the event to move the individual to the next location.
 */
#ifndef CIRCULATETOTARGETLOCATIONNEXTDAYEVENT_H
#define CIRCULATETOTARGETLOCATIONNEXTDAYEVENT_H

// #include "Core/ObjectPool.h"
#include "Core/types.h"
#include "Event.h"

class Person;
class Scheduler;

class CirculateToTargetLocationNextDayEvent : public PersonEvent {
  //  OBJECTPOOL(CirculateToTargetLocationNextDayEvent)
public:
  // disallow copy and move
  CirculateToTargetLocationNextDayEvent(const CirculateToTargetLocationNextDayEvent &) = delete;
  CirculateToTargetLocationNextDayEvent &operator=(const CirculateToTargetLocationNextDayEvent &) =
      delete;
  CirculateToTargetLocationNextDayEvent(CirculateToTargetLocationNextDayEvent &&) = delete;
  CirculateToTargetLocationNextDayEvent &operator=(CirculateToTargetLocationNextDayEvent &&) =
      delete;

  explicit CirculateToTargetLocationNextDayEvent(Person* person) : PersonEvent(person) {}
  ~CirculateToTargetLocationNextDayEvent() override = default;

  [[nodiscard]] core::LocationId target_location() const { return target_location_; }
  void set_target_location(core::LocationId value) { target_location_ = value; }

  [[nodiscard]] const std::string name() const override {
    return "CirculateToTargetLocationNextDayEvent";
  }

private:
  core::LocationId target_location_{0};
  void do_execute() override;
};

#endif
