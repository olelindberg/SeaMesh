#include <iostream>

#include "../../BoostPolygonTypes.hpp"
#include "../../VtkRunner.hpp"

int main() {

  std::cout << "//----------------------------------------------------//\n";
  std::cout << "// test02_boost_polygon_voronoi                       //\n";
  std::cout << "//----------------------------------------------------//\n";

  std::vector<Point> points;
  std::vector<Segment> segments;
  int testid = 2;
  switch (testid) {
  case 0:
    std::cout << "test 0\n";
    segments.push_back(Segment(Point(0.0, 0.0), Point(1.0, 0.0)));
    segments.push_back(Segment(Point(1.0, 0.0), Point(1.0, 1.0)));
    segments.push_back(Segment(Point(1.0, 1.0), Point(0.0, 1.0)));
    segments.push_back(Segment(Point(0.0, 1.0), Point(0.0, 0.0)));
    break;
  case 1:
    std::cout << "test 1\n";
    segments.push_back(Segment(Point(0.0, 0.0), Point(0.0, 1.0)));
    segments.push_back(Segment(Point(0.0, 1.0), Point(1.0, 1.0)));
    segments.push_back(Segment(Point(1.0, 1.0), Point(1.0, 0.0)));
    segments.push_back(Segment(Point(1.0, 0.0), Point(0.0, 0.0)));
    break;
  case 2:
    std::cout << "test 2\n";
    segments.push_back(Segment(Point(0.0, 0.0), Point(0.0, 1.0)));
    segments.push_back(Segment(Point(0.0, 1.0), Point(1.0, 1.0)));
    segments.push_back(Segment(Point(1.0, 1.0), Point(1.0, 0.0)));
    segments.push_back(Segment(Point(1.0, 0.0), Point(0.0, 0.0)));

    segments.push_back(Segment(Point(2.0, 0.0), Point(2.0, 1.0)));
    segments.push_back(Segment(Point(2.0, 1.0), Point(3.0, 1.0)));
    segments.push_back(Segment(Point(3.0, 1.0), Point(3.0, 0.0)));
    segments.push_back(Segment(Point(3.0, 0.0), Point(2.0, 0.0)));

    segments.push_back(Segment(Point(2.0, 2.0), Point(2.0, 3.0)));
    segments.push_back(Segment(Point(2.0, 3.0), Point(3.0, 3.0)));
    segments.push_back(Segment(Point(3.0, 3.0), Point(3.0, 2.0)));
    segments.push_back(Segment(Point(3.0, 2.0), Point(2.0, 2.0)));

    segments.push_back(Segment(Point(-1.0, -1.0), Point(-1.0, 4.0)));
    segments.push_back(Segment(Point(-1.0, 4.0), Point(4.0, 4.0)));
    segments.push_back(Segment(Point(4.0, 4.0), Point(4.0, -1.0)));
    segments.push_back(Segment(Point(4.0, -1.0), Point(-1.0, -1.0)));
    break;
  default:
    break;
  }

  boost::polygon::voronoi_diagram<double> vd;
  boost::polygon::construct_voronoi(segments.begin(), segments.end(), &vd);

  std::vector<vtkSmartPointer<vtkActor>> actors;
  std::vector<vtkSmartPointer<vtkFollower>> followers;

  for (auto it = vd.edges().begin(); it != vd.edges().end(); ++it) {
    if (it->is_primary()) {

      if (it->vertex0() && it->vertex1()) {
        double p0[3] = {it->vertex0()->x(), it->vertex0()->y(), 0.0};
        double p1[3] = {it->vertex1()->x(), it->vertex1()->y(), 0.0};

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
    }

    if (false && it->is_secondary()) {

      const boost::polygon::voronoi_vertex<double> *vertex0 = nullptr;
      const boost::polygon::voronoi_vertex<double> *vertex1 = nullptr;

      if (it->vertex0()) {
        vertex0 = it->vertex0();
        vertex1 = it->next()->vertex1();
      } else if (it->vertex1()) {
        vertex0 = it->prev()->vertex0();
        vertex1 = it->vertex1();
      }

      if (vertex0 && vertex1) {

        double p0[3] = {vertex0->x(), vertex0->y(), 0.0};
        double p1[3] = {vertex1->x(), vertex1->y(), 0.0};

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
        actor->GetProperty()->SetColor(colors->GetColor3d("Green").GetData());

        actors.push_back(actor);
      }
    }
  }

  std::cout << "Number of vtk actors: " << actors.size() << std::endl;

  VtkRunner vtk_runner(actors, followers);

  std::cout << "//----------------------------------------------------//\n";
  std::cout << "// test02_boost_polygon_voronoi, done                 //\n";
  std::cout << "//----------------------------------------------------//\n";

  return 0;
}