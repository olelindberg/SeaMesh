#ifndef SHAPE_FILE_UTIL_HPP
#define SHAPE_FILE_UTIL_HPP

#include <shapefil.h>

#include <cmath>
#include <limits>

class ShapeFileUtil
{
public:
  ///
  /// Check if polygon is clockwise orientated.
  ///
  static bool isClockwisePolygon(SHPObject &shpobj)
  {

    double angle_sum = 0.0;

    int  i  = shpobj.nVertices - 2;
    auto x1 = shpobj.padfX[i + 1] - shpobj.padfX[i];
    auto y1 = shpobj.padfY[i + 1] - shpobj.padfY[i];
    for (int i = 0; i < shpobj.nVertices - 1; ++i)
    {

      auto x2 = shpobj.padfX[i + 1] - shpobj.padfX[i];
      auto y2 = shpobj.padfY[i + 1] - shpobj.padfY[i];

      auto dot = x1 * x2 + y1 * y2; // Dot product between[x1, y1] and [x2, y2]
      auto det = x1 * y2 - y1 * x2; // Determinant
      angle_sum += atan2(det, dot); // atan2(y, x) or atan2(sin, cos)

      x1 = x2;
      y1 = y2;
    }

    bool isclockwise = true;
    if (angle_sum > 0.0)
      isclockwise = false;

    return isclockwise;
  }

  ///
  /// Merge open polygons
  ///
  static void mergeOpenPolygons(const std::vector<SHPObject> &shape_polygons, std::vector<SHPObject> &shape_polygons_closed)
  {

    //-------------------------------------------------------------------------//
    // Find open polygons:
    //-------------------------------------------------------------------------//
    std::vector<SHPObject> shape_polygons_open;
    for (auto &shape_polygon : shape_polygons)
    {

      if (!ShapeFileUtil::isClosedPolygon(shape_polygon))
        shape_polygons_open.push_back(shape_polygon);
      else
        shape_polygons_closed.push_back(shape_polygon);
    }

    //-------------------------------------------------------------------------//
    // Merge open polygons:
    //-------------------------------------------------------------------------//
    for (auto polygon1 = shape_polygons_open.begin(); polygon1 != shape_polygons_open.end();)
    {
      for (auto polygon2 = shape_polygons_open.begin(); polygon2 != shape_polygons_open.end();)
      {

        auto iterator_updated = false;
        if (polygon1 != polygon2)
        {

          if (distance(polygon1->padfX[0], polygon1->padfY[0], polygon2->padfX[0], polygon2->padfY[0]) < std::numeric_limits<double>::epsilon())
          {
            std::cout << "We have a match1\n";
            std::cout << polygon1->padfX[0] << " " << polygon1->padfY[0] << std::endl;
            std::cout << polygon2->padfX[0] << " " << polygon2->padfY[0] << std::endl;
          }

          else if (distance(polygon1->padfX[0], polygon1->padfY[0], polygon2->padfX[polygon2->nVertices - 1], polygon2->padfY[polygon2->nVertices - 1]) <
                   std::numeric_limits<double>::epsilon())
          {

            // Create new polygon:
            SHPObject polygon_new;
            ShapeFileUtil::merge(*polygon2, *polygon1, polygon_new);

            // Assign the new polygon and erase the old:
            *polygon1        = polygon_new;
            polygon2         = shape_polygons_open.erase(polygon2);
            iterator_updated = true;
          }

          else if (distance(polygon1->padfX[polygon1->nVertices - 1], polygon1->padfY[polygon1->nVertices - 1], polygon2->padfX[polygon2->nVertices - 1],
                            polygon2->padfY[polygon2->nVertices - 1]) < std::numeric_limits<double>::epsilon())
          {
            std::cout << "We have a match3\n";
            std::cout << polygon1->padfX[polygon1->nVertices - 1] << " " << polygon1->padfY[polygon1->nVertices - 1] << std::endl;
            std::cout << polygon2->padfX[polygon2->nVertices - 1] << " " << polygon2->padfY[polygon2->nVertices - 1] << std::endl;
          }

          else if (distance(polygon1->padfX[polygon1->nVertices - 1], polygon1->padfY[polygon1->nVertices - 1], polygon2->padfX[0], polygon2->padfY[0]) <
                   std::numeric_limits<double>::epsilon())
          {

            // Merge and create new polygon:
            SHPObject polygon_new;
            ShapeFileUtil::merge(*polygon1, *polygon2, polygon_new);

            *polygon1        = polygon_new;
            polygon2         = shape_polygons_open.erase(polygon2);
            iterator_updated = true;
          }
        }
        if (!iterator_updated)
          ++polygon2;
      }
      ++polygon1;
    }

    //-------------------------------------------------------------------------//
    // Add merged open polygons, which are now closed, to closed polygons:
    //-------------------------------------------------------------------------//
    for (auto &shape_polygon : shape_polygons_open)
    {
      if (ShapeFileUtil::isClosedPolygon(shape_polygon))
      {
        shape_polygons_closed.push_back(shape_polygon);
      }
    }
  }

  static std::string shape_type_to_string(int shapeType)
  {
    switch (shapeType)
    {
    case 0:
      return "Null Shape";
    case 1:
      return "Point";
    case 3:
      return "PolyLine";
    case 5:
      return "Polygon";
    case 8:
      return "MultiPoint";
    case 11:
      return "PointZ";
    case 13:
      return "PolyLineZ";
    case 15:
      return "PolygonZ";
    case 18:
      return "MultiPointZ";
    case 21:
      return "PointM";
    case 23:
      return "PolyLineM";
    case 25:
      return "PolygonM";
    case 28:
      return "MultiPointM";
    case 31:
      return "MultiPatch";
    default:
      return "Unknown";
    }
  }

  static std::string read_prj(const std::string &filename)
  {
    std::ifstream file(filename);
    if (!file.is_open())
    {
      std::cerr << "Could not open .prj file: " << filename << std::endl;
      return "";
    }
    return std::string((std::istreambuf_iterator<char>(file)), std::istreambuf_iterator<char>());
  }

private:
  static double distance(double x1, double y1, double x2, double y2)
  {
    auto dx = x2 - x1;
    auto dy = y2 - y1;
    return std::sqrt(dx * dx + dy * dy);
  }

  static void merge(const SHPObject &polygon1, const SHPObject &polygon2, SHPObject &polygon_new)
  {

    polygon_new.nVertices = polygon1.nVertices + polygon2.nVertices - 1;
    polygon_new.padfX     = new double[polygon_new.nVertices];
    polygon_new.padfY     = new double[polygon_new.nVertices];

    int id = 0;
    for (int i = 0; i < polygon1.nVertices; ++i)
    {
      polygon_new.padfX[id] = polygon1.padfX[i];
      polygon_new.padfY[id] = polygon1.padfY[i];
      ++id;
    }
    for (int i = 1; i < polygon2.nVertices; ++i)
    {
      polygon_new.padfX[id] = polygon2.padfX[i];
      polygon_new.padfY[id] = polygon2.padfY[i];
      ++id;
    }
  }

  ///
  /// Check if the polygon is closed.
  ///
  static bool isClosedPolygon(const SHPObject &shpobj)
  {

    auto dx   = shpobj.padfX[shpobj.nVertices - 1] - shpobj.padfX[0];
    auto dy   = shpobj.padfY[shpobj.nVertices - 1] - shpobj.padfY[0];
    auto dist = std::sqrt(dx * dx + dy * dy);

    bool isclosed = true;
    if (dist > std::numeric_limits<double>::epsilon())
      isclosed = false;

    return isclosed;
  }
};

#endif // SHAPE_FILE_UTIL_HPP
