#ifndef CONVERSION_UTIL_H
#define CONVERSION_UTIL_H

#include "BoostGeometryTypes.h"

#include <Eigen/Dense>

class ConversionUtil {
public:
  static Eigen::Vector3d toEigen(point_t &p) {
    return Eigen::Vector3d(p.get<0>(), p.get<1>(), 0.0);
  }
};
#endif // CONVERSION_UTIL_H