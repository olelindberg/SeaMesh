#ifndef GEOSPATIAL_UTIL_H
#define GEOSPATIAL_UTIL_H

#include "BoostGeometryTypes.h"
#include "logger.h"

#include "ogr_spatialref.h"
#include <proj.h>

#include <iostream>

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

  static bool transform_multi_polygon(multi_polygon_t &mp, const std::string &src_epsg = "EPSG:25832", const std::string &dst_epsg = "EPSG:3034")
  {
    PJ_CONTEXT *ctx = proj_context_create();
    if (!ctx)
      return false;

    PJ *crs_trans = proj_create_crs_to_crs(ctx, src_epsg.c_str(), dst_epsg.c_str(), nullptr);
    if (!crs_trans)
    {
      proj_context_destroy(ctx);
      return false;
    }

    // Iterate over all polygons, rings, and points
    for (auto &poly : mp)
    {
      for (auto &ring : poly.inners())
      {
        // interior rings
        for (auto &pt : ring)
        {
          PJ_COORD c = proj_coord(bg::get<0>(pt), bg::get<1>(pt), 0, 0);
          c          = proj_trans(crs_trans, PJ_FWD, c);
          bg::set<0>(pt, c.xy.x);
          bg::set<1>(pt, c.xy.y);
        }
      }

      auto &outer_ring = poly.outer();
      for (auto &pt : outer_ring)
      {
        PJ_COORD c = proj_coord(bg::get<0>(pt), bg::get<1>(pt), 0, 0);
        c          = proj_trans(crs_trans, PJ_FWD, c);
        bg::set<0>(pt, c.xy.x);
        bg::set<1>(pt, c.xy.y);
      }
    }

    proj_destroy(crs_trans);
    proj_context_destroy(ctx);
    return true;
  }
};

#endif // GEOSPATIAL_UTIL_H