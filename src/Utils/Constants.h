//
// Created by Kien Tran on 1/8/25.
//

#ifndef CONSTANTS_H
#define CONSTANTS_H

#include "Core/types.h"
class Constants {
public:
  static constexpr double PI = 3.14159265358979323846;
  static constexpr core::SimDay DAYS_IN_YEAR = 365;
  // 100 parasites total, equivalent to 0.00002 parasites per microliter.
  static constexpr double DEFAULT_LOG10_PARASITE_DENSITY_CURED = -4.699;
};

#endif  // CONSTANTS_H
