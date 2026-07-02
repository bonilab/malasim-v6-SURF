/*
 * ReportTreatmentFailureDeathEvent.h
 *
 * Defines the event that reports that an individual died of malaria following
 * treatment.
 */
#ifndef REPORTTREATMENTFAILUREDEATHEVENT_H
#define REPORTTREATMENTFAILUREDEATHEVENT_H

// #include "Core/ObjectPool.h"
#include "Core/types.h"
#include "Event.h"

class Person;
class Scheduler;

class ReportTreatmentFailureDeathEvent : public PersonEvent {
  // OBJECTPOOL(ReportTreatmentFailureDeathEvent)
public:
  ReportTreatmentFailureDeathEvent &operator=(const ReportTreatmentFailureDeathEvent &) = delete;
  ReportTreatmentFailureDeathEvent &operator=(ReportTreatmentFailureDeathEvent &&) = delete;
  // disallow copy and move
  ReportTreatmentFailureDeathEvent(const ReportTreatmentFailureDeathEvent &) = delete;
  ReportTreatmentFailureDeathEvent(ReportTreatmentFailureDeathEvent &&) = delete;

  explicit ReportTreatmentFailureDeathEvent(Person* person)
      : PersonEvent(person), age_class_(0), location_id_(0), therapy_id_(0) {}
  ~ReportTreatmentFailureDeathEvent() override = default;

  [[nodiscard]] const std::string name() const override {
    return "ReportTreatmentFailureDeathEvent";
  }

  [[nodiscard]] core::AgeClass age_class() const { return age_class_; }
  void set_age_class(core::AgeClass value) { age_class_ = value; }
  [[nodiscard]] int location_id() const { return location_id_; }
  void set_location_id(int value) { location_id_ = value; }
  [[nodiscard]] int therapy_id() const { return therapy_id_; }
  void set_therapy_id(int value) { therapy_id_ = value; }

private:
  core::AgeClass age_class_{core::K_INVALID_AGE_CLASS};
  int location_id_{-1};
  int therapy_id_{-1};
  void do_execute() override;
};

#endif
