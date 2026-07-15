#ifndef SWITCH_IMMUNE_SYSTEM_MODE_EVENT_H
#define SWITCH_IMMUNE_SYSTEM_MODE_EVENT_H

#include "Event.h"

class Person;

class SwitchImmuneSystemModeEvent : public PersonEvent {
public:
  SwitchImmuneSystemModeEvent(const SwitchImmuneSystemModeEvent &) = delete;
  SwitchImmuneSystemModeEvent(SwitchImmuneSystemModeEvent &&) = delete;
  SwitchImmuneSystemModeEvent &operator=(const SwitchImmuneSystemModeEvent &) = delete;
  SwitchImmuneSystemModeEvent &operator=(SwitchImmuneSystemModeEvent &&) = delete;

  explicit SwitchImmuneSystemModeEvent(Person* person);
  ~SwitchImmuneSystemModeEvent() override;

  [[nodiscard]] const std::string name() const override { return "SwitchImmuneSystemModeEvent"; }

protected:
  void do_execute() override;
};

#endif  // SWITCH_IMMUNE_SYSTEM_MODE_EVENT_H
