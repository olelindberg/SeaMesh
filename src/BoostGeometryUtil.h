#ifndef BOOST_GEOMETRY_UTIL_H
#define BOOST_GEOMETRY_UTIL_H

#include "BoostGeometryTypes.h"

#include <Eigen/Dense>

class BoostGeomtryUtil
{
public:
  static void segment_rtree(const multi_polygon_t &multi_polygon, segment_rtree_t &segment_rtree)
  {

    for (size_t p = 0; p < multi_polygon.size(); ++p)
    {
      auto const &poly = multi_polygon[p];

      // Outer ring
      auto const &outer = poly.outer();
      for (size_t i = 0; i + 1 < outer.size(); ++i)
        segment_rtree.insert({segment_t(outer[i], outer[i + 1]), {p, 0, i}});

      // Inner rings
      for (size_t r = 0; r < poly.inners().size(); ++r)
      {
        auto const &inner = poly.inners()[r];
        for (size_t i = 0; i + 1 < inner.size(); ++i)
          segment_rtree.insert({segment_t(inner[i], inner[i + 1]), {p, r + 1, i}});
      }
    }
  }

  static void makeConvexHull(const std::vector<Eigen::Vector3d> &points, polygon_t &hull)
  {
    boost::geometry::model::multi_point<point_t> mpt;
    for (auto &point : points)
      boost::geometry::append(mpt, point_t(point(0), point(1)));
    boost::geometry::convex_hull(mpt, hull);
  }

  static void makePolygon(const std::vector<Eigen::Vector3d> &points, polygon_t &polygon)
  {
    for (auto &point : points)
      boost::geometry::append(polygon, point_t(point(0), point(1)));
  }

  static void removeSmallPolygons(double area_min, multi_polygon_t &multi_polygon)
  {
    for (auto it = multi_polygon.begin(); it != multi_polygon.end();)
    {
      if (boost::geometry::area(it->outer()) < area_min)
      {
        it = multi_polygon.erase(it);
      }
      else
      {
        ++it;
      }
    }
  }

  static void printMultiPolygon(const std::shared_ptr<multi_polygon_t> &multi_polygon);
};

#endif // BOOST_GEOMETRY_UTIL_H