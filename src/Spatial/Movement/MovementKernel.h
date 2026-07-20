#ifndef SPATIAL_MOVEMENT_MOVEMENTKERNEL_H
#define SPATIAL_MOVEMENT_MOVEMENTKERNEL_H

#include <cmath>
#include <cstddef>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>

#include "Spatial/GIS/GridPairTable.h"

namespace Spatial {
template <typename Func>
[[nodiscard]] GridPairTable make_grid_movement_kernel(const GridPairTable &distances,
                                                      Func func,
                                                      const std::string &context) {
  constexpr double K_EPSILON = std::numeric_limits<double>::epsilon();
  const auto reads_as_zero = [](double value) { return std::fabs(value) < K_EPSILON; };

  auto kernel = distances.map(
      [&](double distance) { return reads_as_zero(distance) ? 0.0 : func(distance); });

  const auto &distance_values = distances.stored_values();
  const auto &kernel_values = kernel.stored_values();
  for (size_t i = 0; i < distance_values.size(); ++i) {
    if (reads_as_zero(distance_values[i]) || !reads_as_zero(kernel_values[i])) { continue; }

    std::ostringstream message;
    message << context << ": the zero sentinel is ambiguous. A pair at distance "
            << distance_values[i] << " maps to " << kernel_values[i]
            << ", which is_zero() cannot distinguish from the zero-distance sentinel (epsilon "
            << K_EPSILON << "). Check the model parameters against the raster extent.";
    throw std::runtime_error(message.str());
  }
  return kernel;
}
}  // namespace Spatial

#endif  // SPATIAL_MOVEMENT_MOVEMENTKERNEL_H
