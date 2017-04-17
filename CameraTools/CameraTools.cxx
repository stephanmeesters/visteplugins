#include "CameraTools.h"

#define VTK_CREATE(type, name) \
  vtkSmartPointer<type> name = vtkSmartPointer<type>::New()

namespace bmia
{

typedef struct
{
	bool used;
	double position[3];
	double focalPoint[3];
	double viewAngle;
	double viewUp[3];
	double clippingRange[2];
	double parallelScale;
	double roll;

	vtkTransform* transform;

} SavedCameraSettings;

// Define interaction style
class KeyPressInteractorStyle : public vtkInteractorStyleTrackballCamera
{
private:

	bmia::Core* core;
	bool shiftDown;
	QList<SavedCameraSettings*> savedSettingsList;

  public:
    static KeyPressInteractorStyle* New();
    vtkTypeMacro(KeyPressInteractorStyle, vtkInteractorStyleTrackballCamera);

	KeyPressInteractorStyle()
	{
		shiftDown = false;

		for(int i = 0; i<10; i++)
		{
			SavedCameraSettings* savedCameraSettings = new SavedCameraSettings;
			savedCameraSettings->used = false;
			savedSettingsList.append(savedCameraSettings);
		}
	}

	void SetCore(bmia::Core* core)
	{
		this->core = core;
	}
 
    virtual void OnKeyPress() 
    {
		// Get the keypress
		vtkRenderWindowInteractor *rwi = this->Interactor;
		std::string key = rwi->GetKeySym();
		shiftDown = rwi->GetShiftKey();
 
		// Output the key that was pressed
		//std::cout << "Pressed " << key << std::endl;
	
		vtkCamera* camera = this->core->canvas()->GetRenderer3D()->GetActiveCamera();

		// Handle an arrow key
		if(key == "Up")
		{
			camera->Elevation(-1);
			camera->OrthogonalizeViewUp();
		}

		if(key == "Down")
		{
			camera->Elevation(1);
			camera->OrthogonalizeViewUp();
		}

		if(key == "Left")
		{
			camera->Azimuth(1);
			camera->OrthogonalizeViewUp();
		}

		if(key == "Right")
		{
			camera->Azimuth(-1);
			camera->OrthogonalizeViewUp();
		}

		if(key == "s" || key == "S")
		{
			double dist = camera->GetDistance();
			double origin[3]; 
			camera->GetFocalPoint(origin);

			camera->SetPosition(dist,origin[1],origin[2]);
			this->core->canvas()->GetRenderer3D()->ResetCameraClippingRange();
			camera->SetRoll(-90.0);
			camera->OrthogonalizeViewUp();
		}

		if(key == "c" || key == "C")
		{
			double dist = camera->GetDistance();
			double origin[3]; 
			camera->GetFocalPoint(origin);

			camera->SetPosition(origin[0],dist,origin[2]);
			this->core->canvas()->GetRenderer3D()->ResetCameraClippingRange();
			camera->SetRoll(-180.0);
			camera->OrthogonalizeViewUp();
		}

		if(key == "a" || key == "A")
		{
			double dist = camera->GetDistance();
			double origin[3]; 
			camera->GetFocalPoint(origin);

			camera->SetPosition(origin[0],origin[1],dist);
			this->core->canvas()->GetRenderer3D()->ResetCameraClippingRange();
			camera->SetRoll(0.0);
			camera->OrthogonalizeViewUp();
		}
 
		// Handle a "normal" key
		if(key == "r" || key == "R")
		{
			std::cout << "Reset camera." << std::endl;
			this->core->canvas()->GetRenderer3D()->ResetCamera();
		}

		// Save camera positions
		for(int i = 0; i<1; i++)
		{
			//if(atoi(key.c_str()) == i)
			
			
				SavedCameraSettings* savedCameraSettings = savedSettingsList.at(i);

				std::cout << shiftDown;

				// save
				if(key == "1")
				{
					camera->GetPosition(savedCameraSettings->position);
					camera->GetFocalPoint(savedCameraSettings->focalPoint);
					savedCameraSettings->viewAngle = camera->GetViewAngle();
					camera->GetViewUp(savedCameraSettings->viewUp);
					camera->GetClippingRange(savedCameraSettings->clippingRange);
					savedCameraSettings->parallelScale = camera->GetParallelScale();
					savedCameraSettings->roll = camera->GetRoll();
					//savedCameraSettings->transform = camera->GetViewTransformObject();
					savedCameraSettings->used = true;
					std::cout << 1;
				}

				// apply
				else if(key == "2" && savedCameraSettings->used)
				{
					camera->SetPosition(savedCameraSettings->position);
					camera->SetFocalPoint(savedCameraSettings->focalPoint);
					camera->SetViewAngle(savedCameraSettings->viewAngle);
					camera->SetViewUp(savedCameraSettings->viewUp);
					camera->SetClippingRange(savedCameraSettings->clippingRange);
					camera->SetParallelScale(savedCameraSettings->parallelScale);
					camera->SetRoll(savedCameraSettings->roll);
					camera->OrthogonalizeViewUp();
					//camera->ApplyTransform(savedCameraSettings->transform);
					//camera->ComputeViewPlaneNormal();
					//camera->OrthogonalizeViewUp();
					std::cout << 2;
				}
			
		}
		

		this->core->render();
 
		// Forward events
		vtkInteractorStyleTrackballCamera::OnKeyPress();
    }

	/*virtual void OnKeyRelease() 
    {
		// Get the keypress
		vtkRenderWindowInteractor *rwi = this->Interactor;
		std::string key = rwi->GetKeySym();
	}*/
 
};
vtkStandardNewMacro(KeyPressInteractorStyle);

//------------------------[ Plugin constructor ]-----------------------\\

CameraTools::CameraTools() : plugin::AdvancedPlugin("CameraTools")
{

}

//------------------------[ Plugin destructor ]-----------------------\\

CameraTools::~CameraTools()
{

}

//------------------------[ Initialization ]-----------------------\\

void CameraTools::init()
{
	vtkSmartPointer<KeyPressInteractorStyle> style = 
		vtkSmartPointer<KeyPressInteractorStyle>::New();

	vtkSubCanvas * subcanvas = this->fullCore()->canvas()->GetSubCanvas3D();
	subcanvas->SetInteractorStyle(style);
	style->SetCurrentRenderer(this->fullCore()->canvas()->GetRenderer3D());

	style->SetCore(this->fullCore());
}

}

Q_EXPORT_PLUGIN2(libCameraTools, bmia::CameraTools)