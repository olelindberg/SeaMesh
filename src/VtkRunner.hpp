#ifndef VTK_UTIL_HPP
#define VTK_UTIL_HPP

#include "vtkAutoInit.h"
VTK_MODULE_INIT(vtkRenderingOpenGL2); // VTK was built with vtkRenderingOpenGL2
VTK_MODULE_INIT(vtkInteractionStyle);

#include "vtkRenderWindow.h"
#include "vtkRenderer.h"
#include <vtkAxesActor.h>
#include <vtkFollower.h>
#include <vtkInteractorStyleTrackballCamera.h>
#include <vtkLineSource.h>
#include <vtkNamedColors.h>
#include <vtkOrientationMarkerWidget.h>
#include <vtkPolyDataMapper.h>
#include <vtkProperty.h>
#include <vtkRenderWindowInteractor.h>

#include <vector>

//-------------------------------------------------------------------------
// 2) Renderer and window stuff:
//-------------------------------------------------------------------------

class VtkRunner {
public:
  VtkRunner(std::vector<vtkSmartPointer<vtkActor>> actors,
            std::vector<vtkSmartPointer<vtkFollower>> followers) {
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
  }
};

#endif // VTK_UTIL_HPP
