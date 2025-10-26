
#include "test_directional_adaptive_octree.h"

#include "BoostGeometryCreator.h"
//#include "BoostGeometryRtree.h"
//#include "BoostGeometryUtil.h"
//#include "ConversionUtil.h"
#include "GeoTiffReader.h"
#include "Logger.h"
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

enum class COASTLINE_ID
{
  TOPO_BATHY_KYSTLINIE,
  CIRCULAR
};

class LandPolygonRefineMask : public RefineMaskSelector
{
public:
  LandPolygonRefineMask(int levelmax, segment_rtree_t &segment_rtree) : _levelmax(levelmax), _segment_rtree(segment_rtree) {}

  virtual RefineMask operator()(const DirectionalAdaptiveOctree *node, double cx, double cy, double cz, int level_x, int level_y, int level_z) const override
  {
    std::cout << "LandPolygonRefineMask: levelx=" << level_x << " levely=" << level_y << " levelz=" << level_z << " box=[(" << node->xmin << "," << node->ymin
              << ")-(" << node->xmax << "," << node->ymax << ")]\n";
    if (level_x >= _levelmax || level_y >= _levelmax || level_z >= _levelmax)
      return RefineMask::REFINE_NONE;

    box_t box(point_t(node->xmin, node->ymin), point_t(node->xmax, node->ymax));

    double center_x = 0.5 * (node->xmin + node->xmax);

    std::vector<segment_value_t> candidates;
    _segment_rtree.query(bgi::intersects(box), std::back_inserter(candidates));

    if (!candidates.empty())
    {
      std::cout << "  found " << candidates.size() << " intersections.\n";
      return RefineMask::REFINE_X | RefineMask::REFINE_Y;
    }
    std::cout << "  no intersections found.\n";
    return RefineMask::REFINE_NONE;
  }

private:
  int              _levelmax;
  segment_rtree_t &_segment_rtree;
};

int main()
{

  //  test_directional_adaptive_octree();

  //  return 0;
  // test_octree("/home/ole/tmp");

  //  MortonTest::run();

  // return 0;

  int polynomial_degree = 2;

  Logger::info("Starting ...");

  // Coordinate transforms from geographical to unit cube:

  // 1) Horizontal Geographical to unit cube:
  // 4) Vertical z to unite cube:

  // Vertical sigma transform:
  // Earth surface projection transform:

  Logger::info("Creating mesh ...");

  // Create a simple oc tree:
  Logger::info("Sea mesh (octree) ...");
  //  auto sea = std::make_shared<OcTree>(Eigen::Vector3d(0.0, 0.0, 0.0), Eigen::Vector3d(1.0, 1.0, 1.0), 0, 0);

  Logger::info("Surface mesh (quadtree) ...");
  //  QuadTree::ptr surface;
  //  TreesUtil::createSurfaceMesh(sea, surface);

  Logger::info("Bottom mesh (quadtree) ...");
  //  QuadTree::ptr seabed;
  //  TreesUtil::createBottomMesh(sea, seabed);

  // Creating static fields:
  Logger::info("Creating Static fields ...");

  // Create bathymetry:
  Logger::info("Bathymetry");

  // Creating all prognostic fields:
  Logger::info("Creating prognostic fields ...");

  // Create surface elevation:
  Logger::info("Surface elevation");

  // Create horizontal velocity field:
  Logger::info("Horizontal velocity");

  // Create salinity:
  Logger::info("Salinity");

  // Create temperature:
  Logger::info("Temperature");

  // Create diagnostic fields:
  Logger::info("Creating diagnostic fields ...");

  // Create vertical velocity:
  Logger::info("Vertical velocity");

  // Surface vind stress:
  Logger::info("Surface wind stress");

  // Surface pressure:
  Logger::info("Surface pressure");

  // Create coriolis field:
  Logger::info("Coriolis field");

  // Create density field:
  Logger::info("Density field");

  // Create turbulent kininetic energy:
  Logger::info("Turbulent kinetic energy");

  // Create turbulent dissipation rate:
  Logger::info("Turbulent dissipation rate");

  //  return 0;

  GeoTiffReader geotiffreader;
  //  vtkSmartPointer<vtkImageActor> actorr;
  // geotiffreader.read(actorr);

  double area_min = 1.0e6;

  // int quadtree_depth = 2;
  double       quadtree_leaf_length   = 0.000001;
  double       lengthmin              = 10;
  int          levelmax               = 12;
  bool         make_quadtree          = true;
  bool         show_quadtree          = true;
  bool         show_direction_vectors = true;
  COASTLINE_ID coastline_id           = COASTLINE_ID::TOPO_BATHY_KYSTLINIE;

  double xmin = 350000.0;
  double ymin = 5900000.0;
  double zmin = 0.0;

  double xmax = 950000.0;
  double ymax = 6500000.0;
  double zmax = 100000.0;

  double lengthx = xmax - xmin;
  double lengthy = ymax - ymin;
  double lengthz = zmax - zmin;
  Logger::info("Domain size x: " + std::to_string(lengthx) + " y: " + std::to_string(lengthy) + " z: " + std::to_string(lengthz));

  Logger::info("Creating coastline polygons ...");
  auto multi_polygon = std::make_shared<multi_polygon_t>();
  if (coastline_id == COASTLINE_ID::TOPO_BATHY_KYSTLINIE)
  {
    //-------------------------------------------------------------------------//
    // Reading shape file polygons:
    //-------------------------------------------------------------------------//
    Logger::info("Reading shape file polygons ...");
    ShapeFileReader shapefilereader;
    // shapefilereader.read("../../data/Topo-bathy/shape_format/data/Kystlinie.shp");
    shapefilereader.read("../../data/seamesh/Kystlinie_fixed.shp");
    auto shape_polygons = shapefilereader.getPolygons();

    //-------------------------------------------------------------------------//
    // Merge open connected polygons to closed polygons:
    //-------------------------------------------------------------------------//
    Logger::info("Merging open polygons ...");
    std::vector<SHPObject> shape_polygons_closed;
    //    ShapeFileUtil::mergeOpenPolygons(*shape_polygons, shape_polygons_closed);

    //-------------------------------------------------------------------------//
    // Create boost geometry multi polygon:
    //-------------------------------------------------------------------------//
    Logger::info("Creating boost geometry multi polygon ...");
    ShapeFileToBoostGeometry::create(*shape_polygons, *multi_polygon);

    //-------------------------------------------------------------------------//
    // Delete islands too small:
    //-------------------------------------------------------------------------//
    Logger::info("Removing small polygons ...");
    BoostGeomtryUtil::removeSmallPolygons(area_min, *multi_polygon);
  }
  else if (coastline_id == COASTLINE_ID::CIRCULAR)
  {
    auto polygon = BoostGeomtryCreator::circularPolygon(650000.0, 6200000.0, 100000.0, 100, false);
    multi_polygon->push_back(*polygon);
  }

  // //  auto polygon = BoostGeomtryCreator::squarePolygon(650000.0, 6200000.0, 300000.0);
  // //  multi_polygon->push_back(*polygon);

  // for (auto &polygon : (*multi_polygon)[0].inners())
  // {
  //   std::cout << boost::geometry::area(polygon) << std::endl;
  // }

  // Logger::info("Building octree ...");
  // OctreeNode root_node(0, 0, xmin, xmax, ymin, ymax, 0.0, dx);

  // OctreeBuilderLandPolygon octree_builder_land_polygon(levelmax, *multi_polygon);
  // octree_builder_land_polygon.build_tree(&root_node);

  // export_to_vtu(&root_node, "/home/ole/tmp/seamesh_octree.vtu");

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

  Logger::info("Writing octree VTU ...");
  DirectionalAdaptiveOctreeUtil::write_vtu(root.get(), "/home/ole/tmp/directional_octree_before_balance.vtu");

  Logger::info("Balancing octree ...");
  DirectionalAdaptiveOctreeUtil::balance_octree(root.get());

  Logger::info("Writing octree VTU ...");
  DirectionalAdaptiveOctreeUtil::write_vtu(root.get(), "/home/ole/tmp/directional_octree_after_balance.vtu");

  return 1;
}
