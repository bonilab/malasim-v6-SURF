#ifndef SMARTMFTSTRATEGY_H
#define SMARTMFTSTRATEGY_H

#include "MFTStrategy.h"

class MFTRebalancingStrategy : public MFTStrategy {
public:
  // disallow copy and assign and move
  MFTRebalancingStrategy(const MFTRebalancingStrategy &) = delete;
  void operator=(const MFTRebalancingStrategy &) = delete;
  MFTRebalancingStrategy(MFTRebalancingStrategy &&) = delete;
  MFTRebalancingStrategy &operator=(MFTRebalancingStrategy &&) = delete;

  int update_duration_after_rebalancing{365};
  int latest_adjust_distribution_time{0};
  int delay_until_actual_trigger{365};
  int next_update_time{0};
  std::vector<double> next_distribution{0};

  MFTRebalancingStrategy();

  //        SmartMFTStrategy(const SmartMFTStrategy & orig);
  ~MFTRebalancingStrategy() override;

  [[nodiscard]] std::string to_string() const override;

  void update_end_of_time_step() override;

  void adjust_started_time_point(const int &current_time) override;
};

#endif /* SMARTMFTSTRATEGY_H */
