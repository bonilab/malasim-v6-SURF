#ifndef NONINFANTIMMUNECOMPONENT
#define    NONINFANTIMMUNECOMPONENT

#include "ImmuneComponent.h"

class NonInfantImmuneComponent : public ImmuneComponent {
  //disallow copy and assign
  NonInfantImmuneComponent(const NonInfantImmuneComponent&) = delete;
  void operator=(const NonInfantImmuneComponent&) = delete;

public:
  NonInfantImmuneComponent(ImmuneSystem *immune_system = nullptr);

  // NonInfantImmuneComponent(const NonInfantImmuneComponent& orig);
  virtual ~NonInfantImmuneComponent();

  virtual double get_decay_rate(const int &age = 0) const override;

  virtual double get_acquire_rate(const int &age = 0) const override;

private:

};

#endif    /* NONINFANTIMMUNECOMPONENT */

