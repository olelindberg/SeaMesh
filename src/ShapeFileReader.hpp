#ifndef SHAPE_FILE_READER_HPP
#define SHAPE_FILE_READER_HPP

#include "BoostGeometryUtil.h"
#include "ShapeFileUtil.hpp"
#include "geospatial_util.h"
#include "logger.h"

#include <shapefil.h>

#include <filesystem>

class ShapeFileReader
{
public:
  bool read(std::string filename)
  {

    SHPHandle shphandle = SHPOpen(filename.c_str(), "rb");

    if (shphandle)
    {
      int    pnEntities;
      int    pnShapeType;
      double padfMinBound[4];
      double padfMaxBound[4];
      SHPGetInfo(shphandle, &pnEntities, &pnShapeType, padfMinBound, padfMaxBound);

      Logger::info("Shapefile reader: Opened shapefile: " + filename);
      Logger::info("Shapefile reader: Number of entities: " + std::to_string(pnEntities));
      Logger::info("Shapefile reader: Shape type: " + ShapeFileUtil::shape_type_to_string(pnShapeType));
      Logger::info("Shapefile reader: Bounding box: [" + std::to_string(padfMinBound[0]) + ", " + std::to_string(padfMinBound[1]) + "] - [" +
                   std::to_string(padfMaxBound[0]) + ", " + std::to_string(padfMaxBound[1]) + "]");

      auto projection_file_path = std::filesystem::path(filename).replace_extension(".prj");
      auto projection           = ShapeFileUtil::read_prj(projection_file_path.string());
      GeospatialUtil::print_projection("Shapefile reader: ", projection);

      if (pnShapeType == SHPT_ARC)
      {

        _polygons = std::make_shared<std::vector<SHPObject>>();

        for (int i = 0; i < pnEntities; ++i)
        {

          auto shpobj = SHPReadObject(shphandle, i);

          if (shpobj->nSHPType == SHPT_ARC && shpobj->nParts == 1)
            _polygons->push_back(*shpobj);
        }
      }
    }
    else
    {
      return false;
    }

    return true;
  }

  const std::shared_ptr<std::vector<SHPObject>> &getPolygons() { return _polygons; }

private:
  std::shared_ptr<std::vector<SHPObject>> _polygons;
};

#endif // SHAPE_FILE_READER_HPP