#ifndef BOOST_GEOMETRY_TYPES_H
#define BOOST_GEOMETRY_TYPES_H

#include <boost/geometry.hpp>

#include <boost/geometry/algorithms/within.hpp>
#include <boost/geometry/geometries/adapted/boost_tuple.hpp>
#include <boost/geometry/geometries/box.hpp>
#include <boost/geometry/geometries/geometries.hpp>
#include <boost/geometry/geometries/multi_polygon.hpp>
#include <boost/geometry/geometries/point.hpp>
#include <boost/geometry/geometries/polygon.hpp>
#include <boost/geometry/index/rtree.hpp>

typedef boost::geometry::model::point<double, 2, boost::geometry::cs::cartesian>
    point_t;

typedef boost::geometry::model::segment<point_t> segment_t;
typedef boost::geometry::model::box<point_t> box_t;

typedef std::pair<box_t, unsigned> value_t;

// Template parmeters:
// 1: Point type.
// 2: Clockwise or counter clockwise?
// 3: Open or closed polygons?
typedef boost::geometry::model::polygon<point_t, false, true>
    polygon_ccw_closed_t;
typedef boost::geometry::model::polygon<point_t, false, false>
    polygon_ccw_open_t;
typedef boost::geometry::model::polygon<point_t, true> polygon_t;

// typedef boost::geometry::model::polygon<point_t, false, false> polygon_t;
typedef boost::geometry::model::multi_polygon<polygon_t> multi_polygon_t;

#endif // BOOST_GEOMETRY_TYPES_H