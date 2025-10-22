
#include "test_directional_adaptive_octree.h"

#include "BoostGeometryCreator.h"
#include "BoostGeometryRtree.h"
#include "BoostGeometryUtil.h"
#include "ConversionUtil.h"
#include "GeoTiffReader.h"
#include "Logger.h"
#include "MortonTest.h"
#include "OcTree.h"
#include "OctreeNode.h"
#include "QuadTree.h"
#include "QuadTreeBuilder.h"
#include "QuadTreePointsInserter.h"
#include "QuadTreeSeaBuilder.h"
#include "QuadTreeSearch.h"
#include "ShapeFileReader.hpp"
#include "ShapeFileToBoostGeometry.hpp"
#include "Timer.h"
#include "TreesUtil.h"
#include "VtkQuadTreeUtil.h"
#include "VtkRunner.hpp"
#include "VtkUtil.h"
#include "octree_builder_land_polygon.h"
#include "test_octree.h"

#include <Eigen/Dense>

#include <boost/geometry.hpp>

#include <iostream>
#include <memory>
#include <vector>

namespace bg = boost::geometry;

enum class COASTLINE_ID
{
  TOPO_BATHY_KYSTLINIE,
  CIRCULAR
};

int main()
{

  test_directional_adaptive_octree();

  return 0;
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
  auto sea = std::make_shared<OcTree>(Eigen::Vector3d(0.0, 0.0, 0.0), Eigen::Vector3d(1.0, 1.0, 1.0), 0, 0);

  Logger::info("Surface mesh (quadtree) ...");
  QuadTree::ptr surface;
  TreesUtil::createSurfaceMesh(sea, surface);

  Logger::info("Bottom mesh (quadtree) ...");
  QuadTree::ptr seabed;
  TreesUtil::createBottomMesh(sea, seabed);

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

  GeoTiffReader                  geotiffreader;
  vtkSmartPointer<vtkImageActor> actorr;
  // geotiffreader.read(actorr);

  double area_min = 1.0e6;

  // int quadtree_depth = 2;
  double       quadtree_leaf_length   = 0.000001;
  double       lengthmin              = 10;
  int          levelmax               = 10;
  bool         make_quadtree          = true;
  bool         show_quadtree          = true;
  bool         show_direction_vectors = true;
  COASTLINE_ID coastline_id           = COASTLINE_ID::TOPO_BATHY_KYSTLINIE;
  double       xmin                   = 350000.0;
  double       ymin                   = 5900000.0;
  double       xmax                   = 950000.0;
  double       ymax                   = 6500000.0;

  double lengthx = xmax - xmin;
  double lengthy = ymax - ymin;
  Logger::info("Domain size x: " + std::to_string(lengthx) + " y: " + std::to_string(lengthy));
  double dx = lengthx / std::pow(2, levelmax);
  double dy = lengthy / std::pow(2, levelmax);
  Logger::info("Minimum cell size dx: " + std::to_string(dx) + " dy: " + std::to_string(dy));

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

  //  auto polygon = BoostGeomtryCreator::squarePolygon(650000.0, 6200000.0, 300000.0);
  //  multi_polygon->push_back(*polygon);

  for (auto &polygon : (*multi_polygon)[0].inners())
  {
    std::cout << boost::geometry::area(polygon) << std::endl;
  }

  // Build quadtree:
  Logger::info("Building quadtree ...");
  QuadTree::ptr quadtree = nullptr;
  quadtree               = std::make_shared<QuadTree>(Eigen::Vector3d(350000.0, 5900000.0, 0.0), Eigen::Vector3d(950000.0, 6500000.0, 0.0), 0, 0);

  Logger::info("Creating quadtree sea builder ...");
  //  QuadTreeSeaBuilder qt_sea_builder(levelmax, lengthmin, multi_polygon);

  Logger::info("Creating quadtree ...");
  //  qt_sea_builder.create(quadtree);

  Logger::info("Balancing quadtree ...");
  //  qt_sea_builder.balance(quadtree);

  OctreeNode root_node(0, 0, xmin, xmax, ymin, ymax, 0.0, dx);

  OctreeBuilderLandPolygon octree_builder_land_polygon(levelmax, *multi_polygon);
  octree_builder_land_polygon.build_tree(&root_node);

  export_to_vtu(&root_node, "/home/ole/tmp/seamesh_octree.vtu");

  return 0;

  //-------------------------------------------------------------------------
  // VTK stuff beyond this point!
  //-------------------------------------------------------------------------
  Logger::info("VTK stuff ...");

  //-------------------------------------------------------------------------
  // 1) Create all actors and followers:
  //-------------------------------------------------------------------------
  std::vector<vtkSmartPointer<vtkActor>>    actors;
  std::vector<vtkSmartPointer<vtkFollower>> followers;

  Logger::info("Creating vtk polylines ...");
  ;
  vtkNew<vtkNamedColors> colors;
  vtkNew<vtkStringArray> colorArray;
  colors->GetColorNames(colorArray);

  {
    int polyline_id = 0;
    for (auto &polygon : (*multi_polygon))
    {
      actors.push_back(VtkUtil::polyline(polygon.outer(), colorArray->GetValue(polyline_id % colorArray->GetMaxId())));
      ++polyline_id;
      // break;
    }
  }
  Logger::info("Creating vtk quadtree ...");
  ;
  if (quadtree && show_quadtree)
  {
    VtkQuadTreeUtil::show_quad_tree(quadtree, actors, followers);
  }
  std::cout << "Number of vtk actors: " << actors.size() << std::endl;

  //  actors.push_back(actorr);
  VtkRunner vtk_runner(actorr, actors, followers);

  return 1;
}
