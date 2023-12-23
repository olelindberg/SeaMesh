#ifndef BOOST_GEOMETRY_UTIL_H
#define BOOST_GEOMETRY_UTIL_H

#include "BoostGeometryTypes.h"

#include <Eigen/Dense>

class BoostGeomtryUtil {
public:
  static void makeConvexHull(const std::vector<Eigen::Vector3d> &points,
                             polygon_t &hull) {
    boost::geometry::model::multi_point<point_t> mpt;
    for (auto &point : points)
      boost::geometry::append(mpt, point_t(point(0), point(1)));
    boost::geometry::convex_hull(mpt, hull);
  }

  static void makePolygon(const std::vector<Eigen::Vector3d> &points,
                          polygon_t &polygon) {
    for (auto &point : points)
      boost::geometry::append(polygon, point_t(point(0), point(1)));
  }
};

#endif // BOOST_GEOMETRY_UTIL_H