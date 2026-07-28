#ifndef ADAPTIVECYCLINGSTRATEGY_H
#define ADAPTIVECYCLINGSTRATEGY_H

#include <vector>

#include "IStrategy.h"

class AdaptiveCyclingStrategy : public IStrategy {
public:
  // disallow copy and assign and move
  AdaptiveCyclingStrategy(const AdaptiveCyclingStrategy &) = delete;
  void operator=(const AdaptiveCyclingStrategy &) = delete;
  AdaptiveCyclingStrategy(AdaptiveCyclingStrategy &&) = delete;
  AdaptiveCyclingStrategy &operator=(AdaptiveCyclingStrategy &&) = delete;

  std::vector<Therapy*> therapy_list;
  int index{0};
  double trigger_value{0.1};
  int delay_until_actual_trigger{365};
  int turn_off_days{365};
  int latest_switch_time{0};

  AdaptiveCyclingStrategy();

  ~AdaptiveCyclingStrategy() override;

  void add_therapy(Therapy* therapy) override;

  virtual void switch_therapy();

  Therapy* get_therapy(Person* person) override;

  [[nodiscard]] std::string to_string() const override;

  void update_end_of_time_step() override;

  void adjust_started_time_point(const int &current_time) override;

  void monthly_update() override;
};

#endif /* ADAPTIVECYCLINGSTRATEGY_H */
