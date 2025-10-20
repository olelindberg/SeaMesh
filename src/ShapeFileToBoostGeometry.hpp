#ifndef SHAPE_FILE_TO_BOOST_GEOMETRY_HPP
#define SHAPE_FILE_TO_BOOST_GEOMETRY_HPP

#include <shapefil.h>

#include <cmath>
#include <limits>

class ShapeFileToBoostGeometry
{
public:
  static void create(const std::vector<SHPObject> &shape_polygons_closed, multi_polygon_t &multi_polygon)
  {
    multi_polygon.resize(shape_polygons_closed.size());
    for (int i = 0; i < shape_polygons_closed.size(); ++i)
    {
      multi_polygon[i].outer().reserve(shape_polygons_closed[i].nVertices);
      for (int j = 0; j < shape_polygons_closed[i].nVertices; ++j)
      {
        boost::geometry::append(multi_polygon[i].outer(), point_t(shape_polygons_closed[i].padfX[j], shape_polygons_closed[i].padfY[j]));
      }
    }
  }
};

#endif // SHAPE_FILE_TO_BOOST_GEOMETRY_HPP
