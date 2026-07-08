#ifndef CHANGEMUTATIONMASKEVENT_H
#define CHANGEMUTATIONMASKEVENT_H


#include "Events/Event.h"
#include "Core/Scheduler/Scheduler.h"
#include <string>
#include <vector>
#include "Simulation/Model.h"

class ChangeMutationMaskEvent  : public WorldEvent {
public:
  //disallow copy, assign and move
  ChangeMutationMaskEvent(const ChangeMutationMaskEvent&) = delete;
  void operator=(const ChangeMutationMaskEvent&) = delete;
  ChangeMutationMaskEvent(ChangeMutationMaskEvent&&) = delete;
  void operator=(ChangeMutationMaskEvent&&) = delete;

public:
  std::vector<bool> mask_;
public:
  explicit ChangeMutationMaskEvent(const std::vector<bool> &mask, const int &at_time = -1);

  ~ChangeMutationMaskEvent() override = default;

  const std::string name() const override {
    return "ChangeMutationMaskEvent";
  }

private:
  void do_execute() override;
};

#endif //CHANGEMUTATIONMASKEVENT_H
