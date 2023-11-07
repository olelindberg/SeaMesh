#include <shapefil.h>

#include <iostream>
#include <memory>
#include <vector>

class Point {
public:
  double x = 0.0;
  double y = 0.0;
};

class AABB {
public:
  Point lower;
  Point upper;
};

class GeometricObject {
public:
  AABB aabb;
};

class Polyline : public GeometricObject {
public:
  std::vector<Point> vertices;
};

class PolylineSet : public GeometricObject {
public:
  std::vector<Polyline> polylines;
};

class ShapeFileReader {
public:
  std::shared_ptr<PolylineSet> polylineset;

  bool read(std::string filename) {

    SHPHandle shphandle = SHPOpen(filename.c_str(), "rb");

    int pnEntities;
    int pnShapeType;
    double padfMinBound[4];
    double padfMaxBound[4];
    SHPGetInfo(shphandle, &pnEntities, &pnShapeType, padfMinBound,
               padfMaxBound);

    if (pnShapeType == SHPT_ARC) {

      polylineset = std::make_shared<PolylineSet>();

      polylineset->aabb.lower.x = padfMinBound[0];
      polylineset->aabb.lower.y = padfMinBound[1];
      polylineset->aabb.upper.x = padfMaxBound[0];
      polylineset->aabb.upper.y = padfMaxBound[1];

      polylineset->polylines.reserve(pnEntities);

      for (int i = 0; i < pnEntities; ++i) {
        auto shpobj = SHPReadObject(shphandle, i);

        if (shpobj->nSHPType == SHPT_ARC && shpobj->nParts == 1) {

          Polyline polyline;
          polyline.aabb.lower.x = shpobj->dfXMin;
          polyline.aabb.lower.y = shpobj->dfYMin;
          polyline.aabb.upper.x = shpobj->dfXMax;
          polyline.aabb.upper.y = shpobj->dfYMax;
          polyline.vertices.reserve(shpobj->nVertices);

          for (int j = 0; j < shpobj->nVertices; ++j) {
            Point point;
            point.x = shpobj->padfX[j];
            point.y = shpobj->padfY[j];
            polyline.vertices.push_back(point);
          }
          polylineset->polylines.push_back(polyline);
        }
      }
    }
    return true;
  }
};

int main() {

  ShapeFileReader shapefilereader;
  shapefilereader.read("../data/Topo-bathy/shape_format/data/Kystlinie.shp");
  auto polylineset = shapefilereader.polylineset;
  std::cout << polylineset->polylines.size() << std::endl;

  return 1;
}
