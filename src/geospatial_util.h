#ifndef GEOSPATIAL_UTIL_H
#define GEOSPATIAL_UTIL_H

#include "logger.h"

#include "ogr_spatialref.h"

#include <string>

class GeospatialUtil
{
public:
  static void print_projection(const std::string &prefix, const std::string &projection_info)
  {

    if (projection_info.size() > 0)
    {
      OGRSpatialReference oSRS;
      if (oSRS.importFromWkt(projection_info.c_str()) == OGRERR_NONE)
      {
        char *pszPrettyWkt = nullptr;
        oSRS.exportToPrettyWkt(&pszPrettyWkt, TRUE);
        Logger::info(prefix + std::string("Projection information:"));
        std::cout << std::endl << pszPrettyWkt << std::endl << std::endl;
        CPLFree(pszPrettyWkt);
      }
    }
    else
      Logger::info("Geotiff reader: No projection found in dataset.");

    // Example function body
  }
};

#endif // GEOSPATIAL_UTIL_H