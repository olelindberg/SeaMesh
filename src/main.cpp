#include "BoostGeometryCreator.h"
#include "BoostGeometryRtree.h"
#include "BoostGeometryUtil.h"
#include "ConversionUtil.h"
#include "QuadTree.h"
#include "QuadTreeBuilder.h"
#include "QuadTreePointsInserter.h"
#include "QuadTreeSearch.h"
#include "ShapeFileReader.hpp"
#include "ShapeFileToBoostGeometry.hpp"
#include "Timer.h"
#include "VtkQuadTreeUtil.h"
#include "VtkRunner.hpp"
#include "VtkUtil.h"

#include <Eigen/Dense>

#include <boost/geometry.hpp>

#include <iostream>
#include <memory>
#include <vector>

namespace bg = boost::geometry;

enum class COASTLINE_ID { TOPO_BATHY_KYSTLINIE, CIRCULAR };

int main() {

  double area_min = 1.0e6;

  // int quadtree_depth = 2;
  double       quadtree_leaf_length   = 0.000001;
  double       lmin                   = 200;
  int          quadtree_subdivisions  = 3;
  bool         make_quadtree          = true;
  bool         show_quadtree          = true;
  bool         show_direction_vectors = true;
  COASTLINE_ID coastline_id           = COASTLINE_ID::TOPO_BATHY_KYSTLINIE;

  auto multi_polygon = std::make_shared<multi_polygon_t>();
  if (coastline_id == COASTLINE_ID::TOPO_BATHY_KYSTLINIE) {

    //-------------------------------------------------------------------------//
    // Reading shape file polygons:
    //-------------------------------------------------------------------------//
    ShapeFileReader shapefilereader;
    shapefilereader.read("../../data/Topo-bathy/shape_format/data/Kystlinie.shp");
    auto shape_polygons = shapefilereader.getPolygons();

    //-------------------------------------------------------------------------//
    // Merge open connected polygons to closed polygons:
    //-------------------------------------------------------------------------//
    std::vector<SHPObject> shape_polygons_closed;
    ShapeFileUtil::mergeOpenPolygons(*shape_polygons, shape_polygons_closed);

    //-------------------------------------------------------------------------//
    // Create boost geometry multi polygon:
    //-------------------------------------------------------------------------//
    ShapeFileToBoostGeometry::create(shape_polygons_closed, *multi_polygon);

    //-------------------------------------------------------------------------//
    // Delete islands too small:
    //-------------------------------------------------------------------------//
    BoostGeomtryUtil::removeSmallPolygons(area_min, *multi_polygon);

  } else if (coastline_id == COASTLINE_ID::CIRCULAR) {

    auto polygon = BoostGeomtryCreator::circularPolygon(650000.0, 6200000.0, 100000.0, 100, false);
    multi_polygon->push_back(*polygon);
  }

  auto polygon = BoostGeomtryCreator::squarePolygon(650000.0, 6200000.0, 300000.0);
  multi_polygon->push_back(*polygon);

  // int poly_idx = 0;
  // for (auto const &poly : *multi_polygon) {
  //   std::cout << "Polygon " << poly_idx++ << ":\n";
  //
  //   // Exterior ring
  //   std::cout << "  Exterior ring:\n";
  //   for (auto const &p : poly.outer()) {
  //     //  std::cout << "    (" << bg::get<0>(p) << ", " << bg::get<1>(p) <<
  //     //  ")\n";
  //   }
  //
  //   // Interior rings (holes)
  //   int hole_idx = 0;
  //   for (auto const &hole : poly.inners()) {
  //     std::cout << "  Hole " << hole_idx++ << ":\n";
  //     for (auto const &p : hole) {
  //       // std::cout << "    (" << bg::get<0>(p) << ", " << bg::get<1>(p) <<
  //       // ")\n";
  //     }
  //   }
  // }

  //  return 0;

  for (auto &polygon : (*multi_polygon)[0].inners()) {
    std::cout << boost::geometry::area(polygon) << std::endl;
  }

  BoostRtreeSearch<polygon_t> search_tree((*multi_polygon));

  // Build quadtree:
  QuadTree::ptr quadtree = nullptr;
  quadtree               = std::make_shared<QuadTree>(Eigen::Vector3d(350000.0, 5900000.0, 0.0), Eigen::Vector3d(950000.0, 6500000.0, 0.0), 0, 0);

  QuadTreePointsInserter qt_point_inserter(20, lmin);

  int poly_idx = 0;
  for (auto const &poly : *multi_polygon) {
    // Exterior ring
    for (auto const &p : poly.outer()) {
      qt_point_inserter.create(quadtree, p);
    }
  }

  qt_point_inserter.balance(quadtree);

  //-------------------------------------------------------------------------
  // VTK stuff beyond this point!
  //-------------------------------------------------------------------------
  std::cout << "VTK stuff ...\n";

  //-------------------------------------------------------------------------
  // 1) Create all actors and followers:
  //-------------------------------------------------------------------------
  std::vector<vtkSmartPointer<vtkActor>>    actors;
  std::vector<vtkSmartPointer<vtkFollower>> followers;

  std::cout << "Creating vtk polylines ...\n";
  vtkNew<vtkNamedColors> colors;
  vtkNew<vtkStringArray> colorArray;
  colors->GetColorNames(colorArray);

  {
    int polyline_id = 0;
    for (auto &polygon : (*multi_polygon)) {
      actors.push_back(VtkUtil::polyline(polygon.outer(), colorArray->GetValue(polyline_id % colorArray->GetMaxId())));
      ++polyline_id;
      // break;
    }
  }
  std::cout << "Creating vtk quadtree ...\n";
  if (quadtree && show_quadtree) {
    VtkQuadTreeUtil::show_quad_tree(quadtree, actors, followers);
  }
  std::cout << "Number of vtk actors: " << actors.size() << std::endl;

  VtkRunner vtk_runner(actors, followers);

  return 1;
}
