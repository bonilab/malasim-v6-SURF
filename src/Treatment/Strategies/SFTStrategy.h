#ifndef SFTSTRATEGY_H
#define SFTSTRATEGY_H

#include <vector>

#include "IStrategy.h"

class SFTStrategy : public IStrategy {
public:
  // disallow copy and assign and move
  SFTStrategy(const SFTStrategy &) = delete;
  void operator=(const SFTStrategy &) = delete;
  SFTStrategy(SFTStrategy &&) = delete;
  SFTStrategy &operator=(SFTStrategy &&) = delete;

  [[nodiscard]] std::vector<Therapy*> get_therapy_list() const { return therapy_list_; }
  void set_therapy_list(const std::vector<Therapy*> &therapy_list) { therapy_list_ = therapy_list; }

  SFTStrategy();

  //    SFTStrategy(const SFTStrategy& orig);
  ~SFTStrategy() override;

  virtual std::vector<Therapy*> &get_therapy_list();

  void add_therapy(Therapy* therapy) override;

  Therapy* get_therapy(Person* person) override;

  [[nodiscard]] std::string to_string() const override;

  void update_end_of_time_step() override;

  void adjust_started_time_point(const int &current_time) override;

  void monthly_update() override;

private:
  std::vector<Therapy*> therapy_list_;
};

#endif /* SFTSTRATEGY_H */
