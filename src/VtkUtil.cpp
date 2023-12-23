#include "VtkUtil.h"

#include <vtkAxesActor.h>
#include <vtkColorTransferFunction.h>
#include <vtkFloatArray.h>
#include <vtkGlyph3D.h>
#include <vtkMaskPoints.h>
#include <vtkNamedColors.h>
#include <vtkOrientationMarkerWidget.h>
#include <vtkPointData.h>
#include <vtkPointSource.h>
#include <vtkPolyData.h>
#include <vtkPolyDataMapper.h>
#include <vtkPolyDataNormals.h>
#include <vtkPolyLine.h>
#include <vtkProperty.h>
#include <vtkSphereSource.h>
#include <vtkVertexGlyphFilter.h>

#include <cmath>
#include <limits>
#include <string>

vtkSmartPointer<vtkActor>
VtkUtil::pointsWithColorScale(const std::vector<Eigen::Vector3d> &pointset,
                              double zscale) {
  vtkSmartPointer<vtkPoints> points = vtkPoints::New();
  for (auto &point : pointset)
    points->InsertNextPoint(point(0), point(1), zscale * point(2));
  return pointsWithColorScale(points, zscale);
}

vtkSmartPointer<vtkActor>
VtkUtil::points(const std::vector<Eigen::Vector3d> &pointset, double radius,
                std::string color) {
  vtkSmartPointer<vtkPoints> vtk_pointset = vtkPoints::New();
  vtkSmartPointer<vtkCellArray> vertices = vtkCellArray::New();
  for (auto &point : pointset) {
    vtkIdType pid[1];
    pid[0] = vtk_pointset->InsertNextPoint(point(0), point(1), point(2));
    vertices->InsertNextCell(1, pid);
  }
  return points(vtk_pointset, vertices, radius, color);
}

vtkSmartPointer<vtkActor>
VtkUtil::points(vtkSmartPointer<vtkPoints> &points,
                vtkSmartPointer<vtkCellArray> &vertices, double radius,
                std::string color) {

  vtkNew<vtkPolyData> point;
  point->SetPoints(points);
  point->SetVerts(vertices);

  vtkNew<vtkPolyDataMapper> mapper;
  mapper->SetInputData(point);

  vtkNew<vtkActor> actor;
  actor->SetMapper(mapper);

  vtkNew<vtkNamedColors> colors;
  actor->GetProperty()->SetColor(colors->GetColor3d(color).GetData());
  actor->GetProperty()->SetPointSize(radius);

  return actor;
}

vtkSmartPointer<vtkActor>
VtkUtil::pointsWithColorScale(vtkSmartPointer<vtkPoints> &points,
                              double zscale) {

  vtkSmartPointer<vtkPolyData> polydata = vtkSmartPointer<vtkPolyData>::New();
  polydata->SetPoints(points);

  // Add distances to each point
  vtkSmartPointer<vtkFloatArray> zcoordinate =
      vtkSmartPointer<vtkFloatArray>::New();
  zcoordinate->SetNumberOfComponents(1);
  zcoordinate->SetName("zcoordinate");

  // Evaluate the signed distance function at all of the grid points
  for (vtkIdType pointId = 0; pointId < points->GetNumberOfPoints();
       ++pointId) {
    double p[3];
    points->GetPoint(pointId, p);
    zcoordinate->InsertNextValue(p[2]);
  }
  polydata->GetPointData()->SetScalars(zcoordinate);

  vtkNew<vtkVertexGlyphFilter> glyph;
  glyph->SetInputData(polydata);
  glyph->Update();

  double colors[5][3] = {{0, 0, 1}, {0, 1, 1}, {0, 1, 0}, {1, 1, 0}, {1, 0, 0}};
  vtkSmartPointer<vtkColorTransferFunction> colorFunction =
      vtkSmartPointer<vtkColorTransferFunction>::New();
  colorFunction->SetUseBelowRangeColor(true);
  colorFunction->SetBelowRangeColor(colors[0]);
  colorFunction->SetNanColor(colors[0]);
  colorFunction->AddRGBPoint(-18.0 * std::fabs(zscale), colors[0][0],
                             colors[0][1], colors[0][2]);
  colorFunction->AddRGBPoint(-14.0 * std::fabs(zscale), colors[1][0],
                             colors[1][1], colors[1][2]);
  colorFunction->AddRGBPoint(-10.0 * std::fabs(zscale), colors[2][0],
                             colors[2][1], colors[2][2]);
  colorFunction->AddRGBPoint(-6.0 * std::fabs(zscale), colors[3][0],
                             colors[3][1], colors[3][2]);
  colorFunction->AddRGBPoint(-2.0 * std::fabs(zscale), colors[4][0],
                             colors[4][1], colors[4][2]);
  colorFunction->Build();

  vtkNew<vtkPolyDataMapper> mapper;
  mapper->SetInputConnection(glyph->GetOutputPort());
  mapper->ScalarVisibilityOn();
  mapper->SetLookupTable(colorFunction);
  mapper->UseLookupTableScalarRangeOn();

  vtkNew<vtkActor> actor;
  actor->SetMapper(mapper);
  actor->GetProperty()->SetPointSize(5);

  return actor;
}

vtkSmartPointer<vtkActor>
VtkUtil::polyline(const std::vector<Eigen::Vector3d> &pointset,
                  std::string color) {
  vtkSmartPointer<vtkPoints> points = vtkPoints::New();
  for (auto &point : pointset)
    points->InsertNextPoint(point(0), point(1), point(2));
  return polyline(points, color);
}

vtkSmartPointer<vtkActor>
VtkUtil::polyline(const std::vector<point_t> &pointset, std::string color) {
  vtkSmartPointer<vtkPoints> points = vtkPoints::New();
  for (auto &point : pointset) {
    points->InsertNextPoint(point.get<0>(), point.get<1>(), 0.0);
  }
  return polyline(points, color);
}

vtkSmartPointer<vtkActor> VtkUtil::polyline(vtkSmartPointer<vtkPoints> &points,
                                            std::string color) {
  vtkNew<vtkPolyLine> polyLine;
  polyLine->GetPointIds()->SetNumberOfIds(points->GetNumberOfPoints());
  for (unsigned int i = 0; i < points->GetNumberOfPoints(); i++)
    polyLine->GetPointIds()->SetId(i, i);

  vtkNew<vtkCellArray> cells;
  cells->InsertNextCell(polyLine);

  vtkNew<vtkPolyData> polyData;
  polyData->SetPoints(points);
  polyData->SetLines(cells);

  // Setup actor and mapper
  vtkNew<vtkPolyDataMapper> mapper;
  mapper->SetInputData(polyData);

  vtkNew<vtkActor> actor;
  actor->SetMapper(mapper);
  vtkNew<vtkNamedColors> colors;
  actor->GetProperty()->SetColor(colors->GetColor3d(color).GetData());
  return actor;
}

vtkSmartPointer<vtkInteractorObserver> VtkUtil::orientationMarker() {
  vtkSmartPointer<vtkAxesActor> axes = vtkSmartPointer<vtkAxesActor>::New();
  vtkSmartPointer<vtkOrientationMarkerWidget> widget =
      vtkSmartPointer<vtkOrientationMarkerWidget>::New();
  widget->SetOutlineColor(0.9300, 0.5700, 0.1300);
  widget->SetOrientationMarker(axes);
  widget->SetViewport(0.0, 0.0, 0.4, 0.4);

  return widget;
}
