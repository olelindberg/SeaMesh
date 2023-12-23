#ifndef SHAPE_FILE_READER_H
#define SHAPE_FILE_READER_H

#include "BoostGeometryUtil.h"
#include "ShapeFileUtil.h"

#include <shapefil.h>

class ShapeFileReader {
public:
  bool read(std::string filename) {

    SHPHandle shphandle = SHPOpen(filename.c_str(), "rb");

    int pnEntities;
    int pnShapeType;
    double padfMinBound[4];
    double padfMaxBound[4];
    SHPGetInfo(shphandle, &pnEntities, &pnShapeType, padfMinBound,
               padfMaxBound);

    if (pnShapeType == SHPT_ARC) {

      _polygons = std::make_shared<std::vector<SHPObject>>();

      for (int i = 0; i < pnEntities; ++i) {

        auto shpobj = SHPReadObject(shphandle, i);

        if (shpobj->nSHPType == SHPT_ARC && shpobj->nParts == 1)
          _polygons->push_back(*shpobj);
      }
    }

    return true;
  }

  const std::shared_ptr<std::vector<SHPObject>> &getPolygons() {
    return _polygons;
  }

private:
  std::shared_ptr<std::vector<SHPObject>> _polygons;
};

#endif // SHAPE_FILE_READER_H