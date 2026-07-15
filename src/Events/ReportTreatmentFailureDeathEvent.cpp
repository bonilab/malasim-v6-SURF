/*
 * ReportTreatmentFailureDeathEvent.cpp
 *
 * Implement the event to report that an individual died of malaria following
 * treatment.
 */
#include "ReportTreatmentFailureDeathEvent.h"

#include "Core/Scheduler/Scheduler.h"
#include "MDC/ModelDataCollector.h"
#include "Population/Person/Person.h"
#include "Simulation/Model.h"

void ReportTreatmentFailureDeathEvent::do_execute() {
  auto* person = get_person();
  if (person == nullptr) { throw std::runtime_error("Person is nullptr"); }
  Model::get_mdc()->record_1_treatment_failure_by_therapy(person->get_location(),
                                                          person->get_age_class(), therapy_id());
}
