#ifndef VTK_QUAD_TREE_UTIL_H
#define VTK_QUAD_TREE_UTIL_H

#include "QuadTree.h"
#include "VtkUtil.h"

#include "vtkVectorText.h"
#include <vtkActor.h>
#include <vtkFollower.h>
#include <vtkNamedColors.h>
#include <vtkPlaneSource.h>
#include <vtkPolyDataMapper.h>
#include <vtkPolyLine.h>
#include <vtkProperty.h>

#include <memory>
#include <vector>

class VtkQuadTreeUtil {
public:
  template <typename point_t>
  static vtkSmartPointer<vtkFollower> indexFollower(int cnt, const point_t &pos,
                                                    double scale) {
    vtkSmartPointer<vtkVectorText> textSource =
        vtkSmartPointer<vtkVectorText>::New();
    textSource->SetText(std::to_string(cnt).c_str());
    textSource->Update();

    // Create a mapper and actor
    vtkSmartPointer<vtkPolyDataMapper> mapper =
        vtkSmartPointer<vtkPolyDataMapper>::New();
    mapper->SetInputConnection(textSource->GetOutputPort());

    vtkSmartPointer<vtkFollower> follower = vtkSmartPointer<vtkFollower>::New();
    follower->SetMapper(mapper);
    follower->GetProperty()->SetColor(0.0, 1.0, 0.0);
    follower->SetPosition(pos(0), pos(1), pos(2));
    follower->SetScale(scale);
    return follower;
  }

  static void
  show_quad_tree(std::shared_ptr<QuadTree> &q,
                 std::vector<vtkSmartPointer<vtkActor>> &actors,
                 std::vector<vtkSmartPointer<vtkFollower>> &followers) {

    vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
    vtkSmartPointer<vtkCellArray> cells = vtkSmartPointer<vtkCellArray>::New();

    make_actors(q, points, cells);

    // Create a polydata to store everything in
    vtkNew<vtkPolyData> polyData;

    // Add the points to the dataset
    polyData->SetPoints(points);

    // Add the lines to the dataset
    polyData->SetLines(cells);

    // Setup actor and mapper
    vtkNew<vtkPolyDataMapper> mapper;
    mapper->SetInputData(polyData);

    vtkNew<vtkNamedColors> colors;

    vtkNew<vtkActor> actor;
    actor->SetMapper(mapper);
    actor->GetProperty()->SetColor(colors->GetColor3d("Tomato").GetData());

    actors.push_back(actor);
  }

  static void make_actors(std::shared_ptr<QuadTree> &q,
                          vtkSmartPointer<vtkPoints> &points,
                          vtkSmartPointer<vtkCellArray> &cells) {
    if (q->trees.empty()) {

      // Create five points.
      double p0[3] = {q->lower(0), q->lower(1), 0.0};
      double p1[3] = {q->upper(0), q->lower(1), 0.0};
      double p2[3] = {q->upper(0), q->upper(1), 0.0};
      double p3[3] = {q->lower(0), q->upper(1), 0.0};

      // Create a vtkPoints object and store the points in it
      int num_points = points->GetNumberOfPoints();
      points->InsertNextPoint(p0);
      points->InsertNextPoint(p1);
      points->InsertNextPoint(p2);
      points->InsertNextPoint(p3);

      //      Eigen::Vector3d centroid;
      //      q->center(centroid);
      //      std::vector<Eigen::Vector3d> centroids;
      //      centroids.push_back(centroid);

      //    followers.push_back(indexFollower(q->get_id(), centroid, 50));

      vtkNew<vtkPolyLine> polyLine;
      polyLine->GetPointIds()->SetNumberOfIds(5);
      for (unsigned int i = 0; i < 4; i++)
        polyLine->GetPointIds()->SetId(i, num_points + i);
      polyLine->GetPointIds()->SetId(4, num_points + 0);

      // Create a cell array to store the lines in and add the lines to it
      cells->InsertNextCell(polyLine);

      //   actors.push_back(VtkUtil::points(centroids, 20.0, "Magenta"));
    } else {
      for (auto &subtree : q->trees)
        if (subtree)
          make_actors(subtree, points, cells);
    }
  }

  static void make_actors(const std::shared_ptr<QuadTree> &q,
                          const std::vector<int> &face_ids,
                          std::vector<vtkSmartPointer<vtkActor>> &actors) {
    if (q) {
      if (std::find(face_ids.begin(), face_ids.end(), q->get_id()) !=
          face_ids.end()) {

        vtkNew<vtkPlaneSource> planeSource;
        planeSource->SetOrigin(q->lower(0), q->lower(1), 0.0);
        planeSource->SetPoint1(q->upper(0), q->lower(1), 0.0);
        planeSource->SetPoint2(q->lower(0), q->upper(1), 0.0);
        planeSource->Update();

        vtkNew<vtkPolyDataMapper> planeMapper;
        planeMapper->SetInputConnection(planeSource->GetOutputPort());

        vtkNew<vtkActor> planeActor;
        planeActor->SetMapper(planeMapper);
        planeActor->GetProperty()->SetOpacity(0.7);

        actors.push_back(planeActor);
      }

      for (auto &qi : q->trees)
        make_actors(qi, face_ids, actors);
    }
  }
};

#endif // VTK_QUAD_TREE_UTIL_H