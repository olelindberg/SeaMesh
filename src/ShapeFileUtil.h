#ifndef SHAPE_FILE_UTIL_H
#define SHAPE_FILE_UTIL_H

#include <shapefil.h>

#include <cmath>
#include <limits>

class ShapeFileUtil {
public:
  ///
  /// Check if the polygon is closed.
  ///
  static bool isClosedPolygon(SHPObject &shpobj) {

    auto dx = shpobj.padfX[shpobj.nVertices - 1] - shpobj.padfX[0];
    auto dy = shpobj.padfY[shpobj.nVertices - 1] - shpobj.padfY[0];
    auto dist = std::sqrt(dx * dx + dy * dy);

    bool isclosed = true;
    if (dist > std::numeric_limits<double>::epsilon())
      isclosed = false;

    return isclosed;
  }

  ///
  /// Check if polygon is clockwise orientated.
  ///
  static bool isClockwisePolygon(SHPObject &shpobj) {

    double angle_sum = 0.0;

    int i = shpobj.nVertices - 2;
    auto x1 = shpobj.padfX[i + 1] - shpobj.padfX[i];
    auto y1 = shpobj.padfY[i + 1] - shpobj.padfY[i];
    for (int i = 0; i < shpobj.nVertices - 1; ++i) {

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

  static void merge(const SHPObject &polygon1, const SHPObject &polygon2,
                    SHPObject &polygon_new) {

    polygon_new.nVertices = polygon1.nVertices + polygon2.nVertices - 1;
    polygon_new.padfX = new double[polygon_new.nVertices];
    polygon_new.padfY = new double[polygon_new.nVertices];

    int id = 0;
    for (int i = 0; i < polygon1.nVertices; ++i) {
      polygon_new.padfX[id] = polygon1.padfX[i];
      polygon_new.padfY[id] = polygon1.padfY[i];
      ++id;
    }
    for (int i = 1; i < polygon2.nVertices; ++i) {
      polygon_new.padfX[id] = polygon2.padfX[i];
      polygon_new.padfY[id] = polygon2.padfY[i];
      ++id;
    }
  }
};

#endif // SHAPE_FILE_UTIL_H
