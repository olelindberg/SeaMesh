
#include "BoostGeometryCreator.h"
#include "land_polygon_refinement_mask.h"
#include "test_directional_adaptive_octree.h"
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

namespace bg = boost::geometry;

int main()
{

  Logger::info("Running SeaMesh octree builder ...");

  double area_min = 1.0e6;

  // int quadtree_depth = 2;
  double quadtree_leaf_length   = 0.000001;
  double lengthmin              = 10;
  int    levelmax               = 8;
  bool   make_quadtree          = true;
  bool   show_quadtree          = true;
  bool   show_direction_vectors = true;

  double xmin = 350000.0;
  double ymin = 5900000.0;
  double zmin = 0.0;

  double xmax = 950000.0;
  double ymax = 6500000.0;
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
  // Reading shape file polygons:
  //-------------------------------------------------------------------------//
  Logger::info("Reading coastline polygons ...");
  ShapeFileReader shapefilereader;
  shapefilereader.read("../../data/seamesh/Kystlinie_fixed.shp");
  auto shape_polygons = shapefilereader.getPolygons();
  Logger::info("Number of shape file polygons: " + std::to_string(shape_polygons->size()));

  //-------------------------------------------------------------------------//
  // Create boost geometry multi polygon:
  //-------------------------------------------------------------------------//
  Logger::info("Creating boost geometry multi polygon ...");
  auto multi_polygon = std::make_shared<multi_polygon_t>();
  ShapeFileToBoostGeometry::create(*shape_polygons, *multi_polygon);

  //-------------------------------------------------------------------------//
  // Delete islands too small:
  //-------------------------------------------------------------------------//
  Logger::info("Removing small polygons ...");
  BoostGeomtryUtil::removeSmallPolygons(area_min, *multi_polygon);

  Logger::info("Number of coastline polygons: " + std::to_string(multi_polygon->size()));

  //-------------------------------------------------------------------------//
  // Directional (anisotropic) octree with Morton codes and balancing.
  // Root domain: [0,1]^3, morton=0, level=0
  //-------------------------------------------------------------------------//
  Logger::info("Building directional adaptive octree ...");
  segment_rtree_t segment_rtree;
  BoostGeomtryUtil::segment_rtree(*multi_polygon, segment_rtree);

  LandPolygonRefineMask choose_refinement_pattern_example(levelmax, segment_rtree);

  std::unique_ptr<DirectionalAdaptiveOctree> root =
      std::make_unique<DirectionalAdaptiveOctree>(0ULL, 0, 0, 0, xmin, xmax, ymin, ymax, zmin, zmax, RefineMask::REFINE_NONE);

  Logger::info("Initial directional build ...");
  DirectionalAdaptiveOctreeUtil::build_directional_recursive(root.get(), levelmax, choose_refinement_pattern_example);

  DirectionalAdaptiveOctreeUtil::balance_octree(root.get());

  Logger::info("Writing octree VTU ...");
  DirectionalAdaptiveOctreeUtil::write_vtu(root.get(), "/home/ole/tmp/seamesh_octree.vtu");

  return 1;
}
