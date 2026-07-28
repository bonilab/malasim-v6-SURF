#ifndef MFTSTRATEGY_H
#define MFTSTRATEGY_H

#include <vector>

#include "IStrategy.h"

class Random;

class Therapy;

class MFTStrategy : public IStrategy {
public:
  // Disallow copy
  MFTStrategy(const MFTStrategy &) = delete;
  MFTStrategy &operator=(const MFTStrategy &) = delete;

  // Disallow move
  MFTStrategy(MFTStrategy &&) = delete;
  MFTStrategy &operator=(MFTStrategy &&) = delete;

  std::vector<Therapy*> therapy_list;
  std::vector<double> distribution;

  MFTStrategy();

  //    MFTStrategy(const MFTStrategy& orig);
  ~MFTStrategy() override;

  void add_therapy(Therapy* therapy) override;

  Therapy* get_therapy(Person* person) override;

  void update_end_of_time_step() override;

  [[nodiscard]] std::string to_string() const override;

  void adjust_started_time_point(const int &current_time) override;

  void monthly_update() override;
};

#endif /* MFTSTRATEGY_H */
