#ifndef SMSH_VTK_UTIL_H
#define SMSH_VTK_UTIL_H

#include "BoostGeometryUtil.h"

#include <vtkActor.h>
#include <vtkArrowSource.h>
#include <vtkConeSource.h>
#include <vtkInteractorObserver.h>
#include <vtkMatrix4x4.h>
#include <vtkMinimalStandardRandomSequence.h>
#include <vtkNamedColors.h>
#include <vtkPoints.h>
#include <vtkPolyDataMapper.h>
#include <vtkPolyDataNormals.h>
#include <vtkPolyLine.h>
#include <vtkRotationalExtrusionFilter.h>
#include <vtkSmartPointer.h>
#include <vtkSplineFilter.h>
#include <vtkTransform.h>
#include <vtkTransformPolyDataFilter.h>

#include <Eigen/Dense>

class VtkUtil {
public:
  /// @brief Creates a vtk actor for visualization of a point set
  /// @param pointset The vector of points
  /// @param radius Radius of the sphere
  /// @param color Color of the shere
  /// @return Pointer to a vtk actor
  static vtkSmartPointer<vtkActor> points(const std::vector<Eigen::Vector3d> &pointset, double radius = 30.0, std::string color = "red");

  /// @brief Creates a vtk actor for visualization of a point set
  /// @param points The points
  /// @param radius Radius of the sphere
  /// @param color Color of the shere
  /// @return Pointer to a vtk actor
  static vtkSmartPointer<vtkActor> points(vtkSmartPointer<vtkPoints> &points, vtkSmartPointer<vtkCellArray> &vertices, double radius = 30.0, std::string color = "red");

  /// @brief Creates a pixel visualization of a point set with color scale
  /// @param pointset
  /// @param zscale
  /// @return
  static vtkSmartPointer<vtkActor> pointsWithColorScale(const std::vector<Eigen::Vector3d> &pointset, double zscale);
  static vtkSmartPointer<vtkActor> pointsWithColorScale(vtkSmartPointer<vtkPoints> &points, double zscale = 1.0);

  static vtkSmartPointer<vtkInteractorObserver> orientationMarker();

  static vtkSmartPointer<vtkActor> polyline(const std::vector<Eigen::Vector3d> &pointset, std::string color = "Tomato");

  static vtkSmartPointer<vtkActor> polyline(const std::vector<point_t> &pointset, std::string color = "Tomato") {
    vtkSmartPointer<vtkPoints> points = vtkPoints::New();
    for (auto &point : pointset) {
      points->InsertNextPoint(point.get<0>(), point.get<1>(), 0.0);
    }
    return polyline(points, color);
  }

  static vtkSmartPointer<vtkActor> polyline(vtkSmartPointer<vtkPoints> &points, std::string color = "Tomato") {
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
    //    actor->GetProperty()->SetColor(colors->GetColor3d(color).GetData());
    return actor;
  }

  static std::vector<vtkSmartPointer<vtkActor>> momentArrowActor(double rotations, double radius = 0.2) {
    int    N      = 16;
    double pi     = std::acos(-1.0);
    double angmax = 2.0 * pi * rotations;
    double x      = std::cos(0.65 * angmax);
    double y      = std::sin(0.65 * angmax);
    double z      = 0.0 * std::fabs(3.0 * radius * angmax / (2.0 * pi));

    vtkNew<vtkConeSource> cone;
    double                cone_height = 0.35 * std::fabs(angmax);
    cone->SetResolution(N);
    cone->SetRadius(2.0 * radius);
    cone->SetHeight(cone_height);
    if (angmax > 0) {
      cone->SetCenter(x - 0.5 * cone_height * y, y + 0.5 * cone_height * x, z);
      cone->SetDirection(-y, x, 0.0);
    } else {
      cone->SetCenter(x + 0.5 * cone_height * y, y - 0.5 * cone_height * x, z);
      cone->SetDirection(y, -x, 0.0);
    }
    cone->Update();

    vtkNew<vtkPolyDataNormals> cone_normals;
    cone_normals->SetInputConnection(cone->GetOutputPort());
    cone_normals->SetFeatureAngle(60);

    vtkNew<vtkPolyDataMapper> cone_mapper;
    cone_mapper->SetInputConnection(cone_normals->GetOutputPort());

    vtkNew<vtkActor> cone_actor;
    cone_actor->SetMapper(cone_mapper);

    vtkNew<vtkPoints> points;
    for (int i = 0; i < N; ++i) {
      double pi    = std::acos(-1.0);
      double theta = 2.0 * pi * double(i) / double(N - 1);
      double x     = radius * std::cos(theta);
      double y     = radius * std::sin(theta);
      points->InsertPoint(i, 1.0 + x, 0.0, y);
    }

    vtkNew<vtkCellArray> poly;
    poly->InsertNextCell(N);
    for (int i = 0; i < N; ++i)
      poly->InsertCellPoint(i);

    vtkNew<vtkPolyData> profile;
    profile->SetPoints(points);
    profile->SetPolys(poly);

    vtkNew<vtkRotationalExtrusionFilter> extrude;
    extrude->SetInputData(profile);
    extrude->SetResolution(360);
    extrude->SetTranslation(z);
    extrude->SetDeltaRadius(0.0);
    extrude->SetAngle(0.65 * angmax * 180.0 / pi);

    vtkNew<vtkPolyDataNormals> normals;
    normals->SetInputConnection(extrude->GetOutputPort());
    normals->SetFeatureAngle(60);

    vtkNew<vtkPolyDataMapper> map;
    map->SetInputConnection(normals->GetOutputPort());

    vtkNew<vtkActor> spring;
    spring->SetMapper(map);

    std::vector<vtkSmartPointer<vtkActor>> actors = {cone_actor, spring};

    return actors;
  }

  static vtkSmartPointer<vtkActor> forceArrowActor(const std::vector<double> &force, double force_scale, double radius = 0.1) {
    vtkNew<vtkNamedColors> colors;

    // Create an arrow.
    vtkNew<vtkArrowSource> arrowSource;
    arrowSource->SetTipResolution(32);
    arrowSource->SetShaftResolution(32);
    arrowSource->SetTipRadius(2.0 * radius);
    arrowSource->SetShaftRadius(radius);

    // Generate a random start and end point
    double startPoint[3];
    double endPoint[3];
    startPoint[0] = force[0];
    startPoint[1] = force[1];
    startPoint[2] = force[2];

    endPoint[0] = force[0];
    endPoint[1] = force[1] + force[4] / force_scale;
    endPoint[2] = force[2];

    // Compute a basis
    double normalizedX[3];
    double normalizedY[3];
    double normalizedZ[3];

    // The X axis is a vector from start to end
    vtkMath::Subtract(endPoint, startPoint, normalizedX);
    double length = vtkMath::Norm(normalizedX);
    vtkMath::Normalize(normalizedX);

    // The Z axis is an arbitrary vector cross X
    vtkNew<vtkMinimalStandardRandomSequence> rng;
    rng->SetSeed(8775070); // For testing.
    double arbitrary[3];
    for (auto i = 0; i < 3; ++i) {
      rng->Next();
      arbitrary[i] = rng->GetRangeValue(-10, 10);
    }
    vtkMath::Cross(normalizedX, arbitrary, normalizedZ);
    vtkMath::Normalize(normalizedZ);

    // The Y axis is Z cross X
    vtkMath::Cross(normalizedZ, normalizedX, normalizedY);
    vtkNew<vtkMatrix4x4> matrix;

    // Create the direction cosine matrix
    matrix->Identity();
    for (auto i = 0; i < 3; i++) {
      matrix->SetElement(i, 0, normalizedX[i]);
      matrix->SetElement(i, 1, normalizedY[i]);
      matrix->SetElement(i, 2, normalizedZ[i]);
    }

    // Apply the transforms
    vtkNew<vtkTransform> transform;
    transform->Translate(startPoint);
    transform->Concatenate(matrix);
    transform->Scale(length, 0.25 * 230.0, 0.25 * 230.0);

    // Transform the polydata
    vtkNew<vtkTransformPolyDataFilter> transformPD;
    transformPD->SetTransform(transform);
    transformPD->SetInputConnection(arrowSource->GetOutputPort());

    // Create a mapper and actor for the arrow

    vtkNew<vtkPolyDataNormals> normals;
    normals->SetInputConnection(transformPD->GetOutputPort());
    normals->SetFeatureAngle(60);

    vtkNew<vtkPolyDataMapper> mapper;
    mapper->SetInputConnection(normals->GetOutputPort());

    vtkNew<vtkActor> actor;
    actor->SetMapper(mapper);

    return actor;
  }
};

#endif // SMSH_VTK_UTIL_H
