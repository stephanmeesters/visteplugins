#ifndef bmia_CameraTools_h
#define bmia_CameraTools_h

#include "DTITool.h"

#include <vtkPolyDataMapper.h>
#include <vtkObjectFactory.h>
#include <vtkActor.h>
#include <vtkSmartPointer.h>
#include <vtkRenderWindow.h>
#include <vtkRenderer.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkPolyData.h>
#include <vtkSphereSource.h>
#include <vtkInteractorStyleTrackballCamera.h>

#include <vtkCamera.h>

#include "gui/MetaCanvas/vtkSubCanvas.h"

#include <vtkTransform.h>

namespace bmia
{

class CameraTools :  public plugin::AdvancedPlugin
{
    Q_OBJECT
    Q_INTERFACES(bmia::plugin::Plugin)
    Q_INTERFACES(bmia::plugin::AdvancedPlugin)

public:

    QString getPluginVersion()
    {
        return "1.0.0";
    }

    CameraTools();
    ~CameraTools();

    virtual void init();

private:


};

}

#endif  // bmia_CameraTools_h
