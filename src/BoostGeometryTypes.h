#ifndef BOOST_GEOMETRY_TYPES_H
#define BOOST_GEOMETRY_TYPES_H

#include <boost/geometry.hpp>

#include <boost/geometry/algorithms/within.hpp>
#include <boost/geometry/geometries/adapted/boost_tuple.hpp>
#include <boost/geometry/geometries/box.hpp>
#include <boost/geometry/geometries/geometries.hpp>
#include <boost/geometry/geometries/linestring.hpp>
#include <boost/geometry/geometries/multi_linestring.hpp>
#include <boost/geometry/geometries/multi_polygon.hpp>
#include <boost/geometry/geometries/point.hpp>
#include <boost/geometry/geometries/polygon.hpp>
#include <boost/geometry/index/rtree.hpp>

namespace bg  = boost::geometry;
namespace bgi = bg::index;

typedef bg::model::point<double, 2, bg::cs::cartesian> point_t;
typedef bg::model::segment<point_t>                    segment_t;
typedef bg::model::box<point_t>                        box_t;
typedef bg::model::ring<point_t>                       ring_t;
typedef bg::model::linestring<point_t>                 linestring_t;

// Template parmeters:
// 1: Point type.
// 2: Clockwise or counter clockwise?
// 3: Open or closed polygons?
typedef bg::model::polygon<point_t, false, true>  polygon_ccw_closed_t;
typedef bg::model::polygon<point_t, false, false> polygon_ccw_open_t;
typedef bg::model::polygon<point_t, true>         polygon_t;

typedef bg::model::multi_polygon<polygon_t>       multi_polygon_t;
typedef bg::model::multi_linestring<linestring_t> multi_linestring_t;

struct SegmentInfo
{
  size_t polygon_index;
  size_t ring_index; // 0 = outer ring, 1..N = inner rings
  size_t segment_index;
};

typedef std::pair<box_t, int>             polygon_value_t;
typedef std::pair<box_t, unsigned>        value_t;
typedef std::pair<segment_t, SegmentInfo> segment_value_t;

typedef bgi::rtree<segment_value_t, bgi::rstar<16>> segment_rtree_t;

#endif // BOOST_GEOMETRY_TYPES_H