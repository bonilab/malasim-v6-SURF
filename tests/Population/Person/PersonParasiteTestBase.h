#ifndef PERSON_PARASITE_TEST_BASE_H
#define PERSON_PARASITE_TEST_BASE_H

#include "PersonTestBase.h"

class PersonParasiteTest : public PersonTestBase {
protected:
  int current_time_ = 30;

  void SetUp() override {
    PersonTestBase::SetUp();
    mock_scheduler_->set_current_time(current_time_);
  }
};

#endif  // PERSON_PARASITE_TEST_BASE_H
