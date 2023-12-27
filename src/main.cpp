#include "BoostGeometryCreator.h"
#include "BoostGeometryRtree.h"
#include "BoostGeometryUtil.h"
#include "ConversionUtil.h"
#include "QuadTree.h"
#include "QuadTreeBuilder.h"
#include "QuadTreeSearch.h"
#include "ShapeFileReader.h"
#include "Timer.h"
#include "VtkQuadTreeUtil.h"

#include "vtkAutoInit.h"
VTK_MODULE_INIT(vtkRenderingOpenGL2); // VTK was built with vtkRenderingOpenGL2
VTK_MODULE_INIT(vtkInteractionStyle);

#include "VtkUtil.h"

#include "vtkRenderWindow.h"
#include "vtkRenderer.h"
#include <vtkAxesActor.h>
#include <vtkFollower.h>
#include <vtkInteractorStyleTrackballCamera.h>
#include <vtkLineSource.h>
#include <vtkOrientationMarkerWidget.h>
#include <vtkRenderWindowInteractor.h>

#include <Eigen/Dense>

#include <boost/geometry.hpp>

#include <iostream>
#include <memory>
#include <vector>

double distance(double x1, double y1, double x2, double y2) {
  auto dx = x2 - x1;
  auto dy = y2 - y1;
  return std::sqrt(dx * dx + dy * dy);
}

int main() {

  double area_min = 1.0e9;

  // int quadtree_depth = 2;
  double quadtree_leaf_length = 0.000001;
  int quadtree_subdivisions = 5;
  bool make_quadtree = true;
  bool show_quadtree = true;
  bool show_direction_vectors = true;

  auto multi_polygon = std::make_shared<multi_polygon_t>();
  if (true) {

    //-------------------------------------------------------------------------//
    // Reading shape file polygons:
    //-------------------------------------------------------------------------//
    ShapeFileReader shapefilereader;
    shapefilereader.read(
        "../../data/Topo-bathy/shape_format/data/Kystlinie.shp");
    auto shape_polygons = shapefilereader.getPolygons();

    //-------------------------------------------------------------------------//
    // Find open polygons:
    //-------------------------------------------------------------------------//
    std::vector<SHPObject> shape_polygons_closed;
    std::vector<SHPObject> shape_polygons_open;
    for (auto &shape_polygon : *shape_polygons) {

      if (!ShapeFileUtil::isClosedPolygon(shape_polygon))
        shape_polygons_open.push_back(shape_polygon);
      else
        shape_polygons_closed.push_back(shape_polygon);
    }

    //-------------------------------------------------------------------------//
    // Merge open polygons:
    //-------------------------------------------------------------------------//
    for (auto polygon1 = shape_polygons_open.begin();
         polygon1 != shape_polygons_open.end();) {
      for (auto polygon2 = shape_polygons_open.begin();
           polygon2 != shape_polygons_open.end();) {

        auto iterator_updated = false;
        if (polygon1 != polygon2) {

          if (distance(polygon1->padfX[0], polygon1->padfY[0],
                       polygon2->padfX[0], polygon2->padfY[0]) <
              std::numeric_limits<double>::epsilon()) {
            std::cout << "We have a match1\n";
            std::cout << polygon1->padfX[0] << " " << polygon1->padfY[0]
                      << std::endl;
            std::cout << polygon2->padfX[0] << " " << polygon2->padfY[0]
                      << std::endl;
          }

          else if (distance(polygon1->padfX[0], polygon1->padfY[0],
                            polygon2->padfX[polygon2->nVertices - 1],
                            polygon2->padfY[polygon2->nVertices - 1]) <
                   std::numeric_limits<double>::epsilon()) {

            // Create new polygon:
            SHPObject polygon_new;
            ShapeFileUtil::merge(*polygon2, *polygon1, polygon_new);

            // Assign the new polygon and erase the old:
            *polygon1 = polygon_new;
            polygon2 = shape_polygons_open.erase(polygon2);
            iterator_updated = true;

          }

          else if (distance(polygon1->padfX[polygon1->nVertices - 1],
                            polygon1->padfY[polygon1->nVertices - 1],
                            polygon2->padfX[polygon2->nVertices - 1],
                            polygon2->padfY[polygon2->nVertices - 1]) <
                   std::numeric_limits<double>::epsilon()) {
            std::cout << "We have a match3\n";
            std::cout << polygon1->padfX[polygon1->nVertices - 1] << " "
                      << polygon1->padfY[polygon1->nVertices - 1] << std::endl;
            std::cout << polygon2->padfX[polygon2->nVertices - 1] << " "
                      << polygon2->padfY[polygon2->nVertices - 1] << std::endl;
          }

          else if (distance(polygon1->padfX[polygon1->nVertices - 1],
                            polygon1->padfY[polygon1->nVertices - 1],
                            polygon2->padfX[0], polygon2->padfY[0]) <
                   std::numeric_limits<double>::epsilon()) {

            // Merge and create new polygon:
            SHPObject polygon_new;
            ShapeFileUtil::merge(*polygon1, *polygon2, polygon_new);

            *polygon1 = polygon_new;
            polygon2 = shape_polygons_open.erase(polygon2);
            iterator_updated = true;
          }
        }
        if (!iterator_updated)
          ++polygon2;
      }
      ++polygon1;
    }

    //-------------------------------------------------------------------------//
    // Add merged open polygons, which are now closed, to closed polygons:
    //-------------------------------------------------------------------------//
    for (auto &shape_polygon : shape_polygons_open) {
      if (ShapeFileUtil::isClosedPolygon(shape_polygon)) {
        shape_polygons_closed.push_back(shape_polygon);
      }
    }

    //-------------------------------------------------------------------------//
    // Create boost geometry multi polygon:
    //-------------------------------------------------------------------------//
    multi_polygon->resize(shape_polygons_closed.size());
    for (int i = 0; i < shape_polygons_closed.size(); ++i) {
      (*multi_polygon)[i].outer().reserve(shape_polygons_closed[i].nVertices);
      for (int j = 0; j < shape_polygons_closed[i].nVertices; ++j) {
        boost::geometry::append((*multi_polygon)[i].outer(),
                                point_t(shape_polygons_closed[i].padfX[j],
                                        shape_polygons_closed[i].padfY[j]));
      }
    }

    //-------------------------------------------------------------------------//
    // Delete islands too small:
    //-------------------------------------------------------------------------//

    for (auto it = multi_polygon->begin(); it != multi_polygon->end();) {
      if (boost::geometry::area(it->outer()) < area_min) {
        it = multi_polygon->erase(it);
      } else {
        ++it;
      }
    }

  } else {

    auto polygon = BoostGeomtryCreator::circularPolygon(650000.0, 6200000.0,
                                                        100000.0, 100, false);
    multi_polygon->push_back(*polygon);
  }

  auto polygon =
      BoostGeomtryCreator::squarePolygon(650000.0, 6200000.0, 300000.0);
  multi_polygon->push_back(*polygon);

  for (auto &polygon : (*multi_polygon)[0].inners()) {
    std::cout << boost::geometry::area(polygon) << std::endl;
  }

  BoostRtreeSearch<polygon_t> search_tree((*multi_polygon));

  // Build quadtree:
  QuadTree::ptr quadtree = nullptr;
  quadtree = std::make_shared<QuadTree>(
      Eigen::Vector3d(350000.0, 5900000.0, 0.0),
      Eigen::Vector3d(950000.0, 6500000.0, 0.0), 0, 0);
  QuadTreeBuilder treebuilder(1, quadtree_leaf_length);
  std::vector<segment_t> segments;
  if (make_quadtree) {
    Timer timer("Quadtree");
    treebuilder.create(quadtree);

    //    treebuilder.balance(quadtree);

    // Assign quad index and count leaf nodes:
    QuadTreeSearch search;
    std::shared_ptr<QuadTreeIndexAssigner> quadtree_index_assigner =
        std::make_shared<QuadTreeIndexAssigner>();
    search.set_search_event(quadtree_index_assigner);
    search.find_leafs(quadtree);
    std::cout << "#leaf nodes    " << quadtree_index_assigner->get_id()
              << std::endl;

    std::cout << "Quadtree bounding box \n";
    std::cout << quadtree->get_xmin() << " " << quadtree->get_ymin()
              << std::endl;
    std::cout << quadtree->get_xmax() << " " << quadtree->get_ymax()
              << std::endl;

    std::cout << "Area covered by quadtree is ";
    std::cout << (quadtree->get_xmax() - quadtree->get_xmin()) / 1000.0
              << "km times ";
    std::cout << (quadtree->get_ymax() - quadtree->get_ymin()) / 1000.0
              << "km\n";
    for (int i = 0; i < quadtree_subdivisions; ++i) {

      segments.clear();

      auto leafevent = std::make_shared<QuadTreeLeafSearch>();
      QuadTreeSearch leafsearch(leafevent);

      auto neighbourevent = std::make_shared<QuadTreeNeighbourSearch>();
      QuadTreeSearch neighboursearch(neighbourevent);

      std::cout << "Finding leafs ..." << std::endl;
      leafevent->clear();
      leafsearch.find_leafs(quadtree);
      auto leafs = leafevent->get_leafs();

      // Remove leafs inside polygons:
      //      for (auto it = leafs.begin(); it != leafs.end(); ++it) {
      //      }

      std::cout << "Finding closest points ..." << std::endl;
      std::vector<Eigen::Vector3d> directions;
      std::vector<Eigen::Vector3d> origins;
      for (auto &leaf : leafs) {

        Eigen::Vector3d tmp;
        leaf->center(tmp);
        point_t center(tmp(0), tmp(1));

        segment_t segment;
        boost::geometry::closest_points(center, *multi_polygon, segment);
        segments.push_back(segment);

        origins.push_back(ConversionUtil::toEigen(segment.first));
        directions.push_back(
            (origins.back() - ConversionUtil::toEigen(segment.second))
                .normalized());
      }

      std::cout << "Levelset uphill refinement ..." << std::endl;
      for (auto &leaf : leafs) {
        Eigen::Vector3d tmp;
        leaf->center(tmp);
        point_t center(tmp(0), tmp(1));
        if (!boost::geometry::within(center, *multi_polygon)) {
          auto direction = directions[leaf->get_id()];

          auto limit = Eigen::Vector3d(1.0, 1.0, 0.0)
                           .normalized()
                           .dot(Eigen::Vector3d::UnitX());

          auto projx = direction.dot(Eigen::Vector3d::UnitX());
          auto projy = direction.dot(Eigen::Vector3d::UnitY());

          neighbourevent->clear();
          neighboursearch.find_neighbours(leaf);
          auto edge_neighbours = neighbourevent->get_edge_neighbours();
          if (projy <= -limit) // South
          {
            if (edge_neighbours[0].empty()) {
              leaf->flag = true;
            } else {
              for (auto &neighbour : edge_neighbours[0]) {
                if (direction.dot(directions[neighbour->get_id()]) <= 0.0) {
                  leaf->flag = true;
                  neighbour->flag = true;
                }
              }
            }
          } else if (projx >= limit) // East
          {
            if (edge_neighbours[1].empty()) {
              leaf->flag = true;
            } else {
              for (auto &neighbour : edge_neighbours[1]) {
                if (direction.dot(directions[neighbour->get_id()]) <= 0.0) {
                  leaf->flag = true;
                  neighbour->flag = true;
                }
              }
            }
          } else if (projy > limit) // North
          {
            if (edge_neighbours[2].empty()) {
              leaf->flag = true;
            } else {
              for (auto &neighbour : edge_neighbours[2]) {
                if (direction.dot(directions[neighbour->get_id()]) <= 0.0) {
                  leaf->flag = true;
                  neighbour->flag = true;
                }
              }
            }
          } else if (projx < limit) // West
          {
            if (edge_neighbours[3].empty()) {
              leaf->flag = true;
            } else {
              for (auto &neighbour : edge_neighbours[3]) {
                if (direction.dot(directions[neighbour->get_id()]) <= 0.0) {
                  leaf->flag = true;
                  neighbour->flag = true;
                }
              }
            }
          }
        } // if not inside
      }   // for leafs

      std::cout << "Subdividing quad tree ..." << std::endl;
      for (auto &leaf : leafs)
        if (leaf->flag) {
          treebuilder.subdivide(leaf);
        }
      std::cout << "Balancing quad tree ..." << std::endl;
      treebuilder.balance(quadtree);

      // Assign quad index and count leaf nodes:
      QuadTreeSearch search;
      std::shared_ptr<QuadTreeIndexAssigner> quadtree_index_assigner =
          std::make_shared<QuadTreeIndexAssigner>();
      search.set_search_event(quadtree_index_assigner);
      search.find_leafs(quadtree);
      std::cout << "#leaf nodes    " << quadtree_index_assigner->get_id()
                << std::endl;
    }
  }

  //-------------------------------------------------------------------------
  // VTK stuff beyond this point!
  //-------------------------------------------------------------------------
  std::cout << "VTK stuff ...\n";

  //-------------------------------------------------------------------------
  // 1) Create all actors and followers:
  //-------------------------------------------------------------------------
  std::vector<vtkSmartPointer<vtkActor>> actors;
  std::vector<vtkSmartPointer<vtkFollower>> followers;

  if (show_direction_vectors)
    for (auto &segment : segments) {

      auto origin = ConversionUtil::toEigen(segment.first);
      auto origin2 = ConversionUtil::toEigen(segment.second);
      auto direction = origin - ConversionUtil::toEigen(segment.second);

      // Create two points, P0 and P1
      //    double p0[3] = {origin(0), origin(1), 0.0};
      //    double p1[3] = {origin(0) + direction(0), origin(1) +
      //    direction(1), 0.0};
      double p0[3] = {origin(0), origin(1), 0.0};
      double p1[3] = {origin2(0), origin2(1), 0.0};

      vtkNew<vtkLineSource> lineSource;
      lineSource->SetPoint1(p0);
      lineSource->SetPoint2(p1);

      // Visualize
      vtkNew<vtkNamedColors> colors;

      vtkNew<vtkPolyDataMapper> mapper;
      mapper->SetInputConnection(lineSource->GetOutputPort());
      vtkNew<vtkActor> actor;
      actor->SetMapper(mapper);
      actor->GetProperty()->SetLineWidth(2);
      actor->GetProperty()->SetColor(colors->GetColor3d("Peacock").GetData());

      actors.push_back(actor);
    }

  std::cout << "Creating vtk points ...\n";
  //  for (auto polyline : polylineset->polylines) {
  //  actors.push_back(VtkUtil::points(polyline.vertices, 2.0));
  // break;
  // }

  std::cout << "Creating vtk polylines ...\n";
  vtkNew<vtkNamedColors> colors;
  vtkNew<vtkStringArray> colorArray;
  colors->GetColorNames(colorArray);
  // std::cout << colorArray << std::endl;
  // int polyline_id = 0;
  // for (auto polyline : polylineset->polylines) {
  //   actors.push_back(VtkUtil::polyline(
  //       polyline.vertices,
  //       colorArray->GetValue(polyline_id % colorArray->GetMaxId())));
  //   ++polyline_id;
  //   // break;
  // }

  {
    int polyline_id = 0;
    for (auto &polygon : (*multi_polygon)) {
      actors.push_back(VtkUtil::polyline(
          polygon.outer(),
          colorArray->GetValue(polyline_id % colorArray->GetMaxId())));
      ++polyline_id;
      // break;
    }
  }
  std::cout << "Creating vtk quadtree ...\n";
  if (quadtree && show_quadtree) {
    VtkQuadTreeUtil::show_quad_tree(quadtree, actors, followers);
  }
  std::cout << "Number of vtk actors: " << actors.size() << std::endl;

  //  if (quadtree && show_quadtree)
  //    VtkQuadTreeUtil::make_actors(quadtree, actors, followers);

  //-------------------------------------------------------------------------
  // 2) Renderer and window stuff:
  //-------------------------------------------------------------------------
  std::cout << "Creating vtk renderer ...\n";
  vtkSmartPointer<vtkRenderer> renderer = vtkRenderer::New();

  std::cout << "Adding vtk actors to renderer ...\n";
  for (auto &actor : actors)
    renderer->AddActor(actor);

  std::cout << "Setting backgound color ...\n";
  renderer->SetBackground(0.2, 0.2, 0.2);

  std::cout << "Resetting camara ...\n";
  //  renderer->ResetCamera();

  std::cout << "Adding followers ...\n";
  for (auto &follower : followers) {
    renderer->AddActor(follower);
    follower->SetCamera(renderer->GetActiveCamera());
  }

  std::cout << "Creating vtk render window ...\n";
  vtkSmartPointer<vtkRenderWindow> renderWindow = vtkRenderWindow::New();
  renderWindow->AddRenderer(renderer);
  renderWindow->SetSize(1200, 900);

  std::cout << "Creating vtk render window ...\n";
  vtkNew<vtkRenderWindowInteractor> renderWindowInteractor;
  renderWindowInteractor->SetRenderWindow(renderWindow);

  std::cout << "Creating interactor style ...\n";
  vtkSmartPointer<vtkInteractorStyleTrackballCamera> interactorStyle =
      vtkInteractorStyleTrackballCamera::New();
  renderWindowInteractor->SetInteractorStyle(interactorStyle);
  renderWindowInteractor->GetInteractorStyle()->SetCurrentRenderer(renderer);

  //    vtkSmartPointer<vtkAxesActor> axes =
  //    vtkSmartPointer<vtkAxesActor>::New();
  //    vtkSmartPointer<vtkOrientationMarkerWidget> widget =
  //    vtkSmartPointer<vtkOrientationMarkerWidget>::New();
  //    widget->SetOutlineColor(0.9300, 0.5700, 0.1300);
  //    widget->SetOrientationMarker(axes);
  //    widget->SetViewport(0.0, 0.0, 0.4, 0.4);
  //    widget->SetInteractor(renderWindowInteractor);
  //    widget->SetEnabled(1);
  //    widget->InteractiveOn();
  //
  std::cout << "Resetting camara ...\n";
  renderer->ResetCamera();

  std::cout << "Rendering ...\n";
  renderWindow->Render();

  std::cout << "Starting vtk window interactor ...\n";
  renderWindowInteractor->Start();

  return 1;
}
