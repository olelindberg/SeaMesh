
#include "BoostGeometryCreator.h"
#include "boost_multipolygon_vtk_writer.h"
#include "geospatial_util.h"
#include "land_polygon_refinement_mask.h"
#include "test_directional_adaptive_octree.h"
#include "test_gpkg_to_boost_multipolygon.h"
//#include "BoostGeometryRtree.h"
//#include "BoostGeometryUtil.h"
//#include "ConversionUtil.h"
#include "GeoTiffReader.h"
#include "logger.h"
//#include "MortonTest.h"
//#include "OcTree.h"
//#include "OctreeNode.h"
//#include "QuadTree.h"
//#include "QuadTreeBuilder.h"
//#include "QuadTreePointsInserter.h"
//#include "QuadTreeSeaBuilder.h"
//#include "QuadTreeSearch.h"
#include "ShapeFileReader.hpp"
#include "ShapeFileToBoostGeometry.hpp"
//#include "Timer.h"
//#include "TreesUtil.h"
//
#include "octree_builder_land_polygon.h"
#include "octree_vtk_writer.h"
//
//#include "test_octree.h"
//
//#include <Eigen/Dense>
//
//#include <boost/geometry.hpp>
//
//#include <iostream>
//#include <memory>
//#include <vector>

// namespace bg = boost::geometry;

int main()
{

  Logger::info("Running SeaMesh octree builder ...");

  std::string projection = "EPSG:3034";
  double      area_min   = 1.0e6;

  // int quadtree_depth = 2;
  double quadtree_leaf_length   = 0.000001;
  double lengthmin              = 10;
  int    levelmax               = 8;
  bool   make_quadtree          = true;
  bool   show_quadtree          = true;
  bool   show_direction_vectors = true;
  //[INFO] min : (3076896.422822, 3882017.010298)
  //[INFO] max : (3420324.370211, 4318581.557322)

  double xmin = 3000000.0;
  double ymin = 3800000.0;
  double zmin = 0.0;

  double xmax = 3500000.0;
  double ymax = 4400000.0;
  double zmax = 100000.0;

  double lengthx = xmax - xmin;
  double lengthy = ymax - ymin;
  double lengthz = zmax - zmin;

  Logger::info("Domain bounds:");
  Logger::info("  xmin: " + std::to_string(xmin) + " xmax: " + std::to_string(xmax));
  Logger::info("  ymin: " + std::to_string(ymin) + " ymax: " + std::to_string(ymax));
  Logger::info("  zmin: " + std::to_string(zmin) + " zmax: " + std::to_string(zmax));
  Logger::info("  lengthx: " + std::to_string(lengthx));
  Logger::info("  lengthy: " + std::to_string(lengthy));

  GeoTiffReader geotiffreader;
  geotiffreader.read();

  //-------------------------------------------------------------------------//
  // Create boost geometry multi polygon:
  //-------------------------------------------------------------------------//
  Logger::info("Creating boost geometry multi polygon ...");
  auto multi_polygon = std::make_shared<multi_polygon_t>();

  test_gpkg_to_boost_multipolygon("../../data/Klimadatastyrelsen/TopografiskLandpolygon/landpolygon.gpkg", *multi_polygon, "landpolygon_2500", projection);

  //-------------------------------------------------------------------------//
  // Delete islands too small:
  //-------------------------------------------------------------------------//
  Logger::info("Removing small polygons ...");
  BoostGeomtryUtil::removeSmallPolygons(area_min, *multi_polygon);

  //-------------------------------------------------------------------------//
  // Quick sanity: print bbox and polygon count
  //-------------------------------------------------------------------------//
  box_t bbox;
  bg::envelope(*multi_polygon, bbox);
  auto min = bbox.min_corner();
  auto max = bbox.max_corner();
  Logger::info("Coastline bounding box:");
  Logger::info("  min: (" + std::to_string(min.get<0>()) + ", " + std::to_string(min.get<1>()) + ")");
  Logger::info("  max: (" + std::to_string(max.get<0>()) + ", " + std::to_string(max.get<1>()) + ")");
  Logger::info("Number of coastline polygons: " + std::to_string(multi_polygon->size()));

  //-------------------------------------------------------------------------//
  // Build segment R-tree for fast intersection tests
  //------------------------------------------------------------------------//
  Logger::info("Building segment R-tree ...");
  segment_rtree_t segment_rtree;
  BoostGeomtryUtil::segment_rtree(*multi_polygon, segment_rtree);

  //-------------------------------------------------------------------------//
  // Directional (anisotropic) octree with Morton codes and balancing.
  // Root domain: [0,1]^3, morton=0, level=0
  //-------------------------------------------------------------------------//
  Logger::info("Building directional adaptive octree ...");
  LandPolygonRefineMask                      choose_refinement_pattern_example(levelmax, segment_rtree);
  std::unique_ptr<DirectionalAdaptiveOctree> root =
      std::make_unique<DirectionalAdaptiveOctree>(0ULL, 0, 0, 0, xmin, xmax, ymin, ymax, zmin, zmax, RefineMask::REFINE_NONE);
  DirectionalAdaptiveOctreeUtil::build_directional_recursive(root.get(), levelmax, choose_refinement_pattern_example);
  DirectionalAdaptiveOctreeUtil::balance_octree(root.get());

  //-------------------------------------------------------------------------//
  // Write outputs:
  //-------------------------------------------------------------------------//
  Logger::info("Writing coastline VTP ...");
  auto pd = BoostMultiPolygonVtkWriter::multi_polygon_to_lines(*multi_polygon);
  BoostMultiPolygonVtkWriter::write_vtp_lines(pd, "/home/ole/tmp/land_boundaries.vtp");

  Logger::info("Writing octree VTU ...");
  DirectionalAdaptiveOctreeUtil::write_vtu(root.get(), "/home/ole/tmp/seamesh_octree.vtu");

  return 1;
}
