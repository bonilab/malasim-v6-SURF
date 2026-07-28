#ifndef POMS_NESTEDSWITCHINGDIFFERENTDISTRIBUTIONBYLOCATION_H
#define POMS_NESTEDSWITCHINGDIFFERENTDISTRIBUTIONBYLOCATION_H

#include "IStrategy.h"
#include "Utils/TypeDef.h"

class Config;

class NestedMFTMultiLocationStrategy : public IStrategy {
public:
  // disallow copy and assign and move
  NestedMFTMultiLocationStrategy(const NestedMFTMultiLocationStrategy &) = delete;
  void operator=(const NestedMFTMultiLocationStrategy &) = delete;
  NestedMFTMultiLocationStrategy(NestedMFTMultiLocationStrategy &&) = delete;
  NestedMFTMultiLocationStrategy &operator=(NestedMFTMultiLocationStrategy &&) = delete;

  NestedMFTMultiLocationStrategy();

  //    NestedSwitchingStrategy(const NestedSwitchingStrategy& orig);
  ~NestedMFTMultiLocationStrategy() override;

  virtual void add_strategy(IStrategy* strategy);

  void add_therapy(Therapy* therapy) override;

  Therapy* get_therapy(Person* person) override;

  [[nodiscard]] std::string to_string() const override;

  /**
   * This function will be executed at end of time step, to check and switch therapy if needed
   */
  void update_end_of_time_step() override;

  void adjust_distribution(const int &time);

  void adjust_started_time_point(const int &current_time) override;

  void monthly_update() override;

  std::vector<IStrategy*> strategy_list;
  DoubleVector2 distribution;
  DoubleVector2 start_distribution;
  DoubleVector2 peak_distribution;
  int starting_time{0};
  int peak_after{0};
};

#endif  // POMS_NESTEDSWITCHINGDIFFERENTDISTRIBUTIONBYLOCATION_H
