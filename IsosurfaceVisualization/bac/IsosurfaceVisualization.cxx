#include "IsosurfaceVisualization.h"

#define VTK_CREATE(type, name) \
  vtkSmartPointer<type> name = vtkSmartPointer<type>::New()

// init values
#define DEFAULT_SMOOTHING 1.0
#define DEFAULT_ALPHA 1.0
#define DEFAULT_SPECULAR 0.00

namespace bmia
{

///
///      INITIALIZATION
///

//------------------------[ Plugin constructor ]-----------------------\\

IsosurfaceVisualization::IsosurfaceVisualization() : plugin::AdvancedPlugin("IsosurfaceVisualization")
{
    this->widget = NULL;
    this->form   = NULL;
}

//------------------------[ Plugin destructor ]-----------------------\\

IsosurfaceVisualization::~IsosurfaceVisualization()
{
    delete this->widget;
    delete this->form;
}

//------------------------[ Initialization ]-----------------------\\

void IsosurfaceVisualization::init()
{
    this->widget = new QWidget();
    this->form = new Ui::IsosurfaceVisualizationForm();
    this->form->setupUi(this->widget);

    // Link events in the GUI to function calls
    this->connectAll();
    this->assembly = vtkPropAssembly::New();

    this->current_modelInfo = NULL;
    this->overlay_modelInfo = NULL;

    // Update the combo box
    this->comboBoxDataChanged();

    // setup picker
    this->setupClippingPlanesPicker();

    // Get the canvas
	vtkMedicalCanvas * canvas = this->fullCore()->canvas();

	QString sliceActorNames[3];
    sliceActorNames[0] = "X Plane";
    sliceActorNames[1] = "Y Plane";
    sliceActorNames[2] = "Z Plane";

	// Loop through all three axes
	for(int axis = 0; axis < 3; ++axis)
	{
        vtkSmartPointer<vtkSphereSource> diskSource =
            vtkSmartPointer<vtkSphereSource>::New();
        diskSource->SetRadius(4);

        vtkSmartPointer<vtkPolyDataMapper> diskMapper =
            vtkSmartPointer<vtkPolyDataMapper>::New();
        diskMapper->SetInputConnection(diskSource->GetOutputPort());

        vtkActor* diskActor =
            vtkActor::New();
        diskActor->SetMapper(diskMapper);
        diskActor->SetPosition(180,216,0);
        diskActor->GetProperty()->SetColor(1,0,0);

        canvas->GetSubCanvas2D(axis)->GetRenderer()->AddActor(diskActor);

        this->pointer2DList.append(diskActor);
	}

	// Add measured points
	for(int i = 0; i < 2; i++)
	{
        MeasuredPoint* point = new MeasuredPoint;

        point->set = false;

        vtkSmartPointer<vtkSphereSource> diskSource =
            vtkSmartPointer<vtkSphereSource>::New();
        diskSource->SetRadius(1);

        vtkSmartPointer<vtkPolyDataMapper> diskMapper =
            vtkSmartPointer<vtkPolyDataMapper>::New();
        diskMapper->SetInputConnection(diskSource->GetOutputPort());

        vtkActor* diskActor =
            vtkActor::New();
        diskActor->SetMapper(diskMapper);

        if(i == 0)
            diskActor->GetProperty()->SetColor(1,0,0);
        else
            diskActor->GetProperty()->SetColor(0,1,0);


        point->sphere = diskActor;

        this->measuredPointList.append(point);
	}

	this->measuredLine = NULL;

    this->scalarBar = NULL; //temp


}

//------------[ Setup pointer interactor for clipping planes ]----------------\\

void IsosurfaceVisualization::setupClippingPlanesPicker()
{
    // Create a new styleTrackballPP
    this->styleTrackballPP = vtkInteractorStyleTrackballPositionPicker::New();

    // Set up the options for styleTrackballPP
    this->styleTrackballPP->SetRenderProcess(this->fullCore()->canvas()->GetRenderer3D());
    this->styleTrackballPP->SetParentClass(this);

    // Set this styleTrackballPP to subcanvas
    vtkSubCanvas * subcanvas = this->fullCore()->canvas()->GetSubCanvas3D();
    subcanvas->SetInteractorStyle(this->styleTrackballPP);

    // Create a cellPicker and set up its options
    vtkCellPicker * cellPicker = vtkCellPicker::New();

    // Set this cellPicker to subcanvas
    subcanvas->GetInteractor()->SetPicker(cellPicker);
    cellPicker->Delete();
}

///
///      DATA I/O
///

//------------------------[ Dataset added ]-----------------------\\

void IsosurfaceVisualization::dataSetAdded(data::DataSet * d)
{
    // Assert the data set pointer (should never be NULL)
    Q_ASSERT(d);

	// Get the kind of the data set
    QString kind = d->getKind();

    // Load scalar volume
    if(kind == "scalar volume")
    {
        // Keep track of the datasets used by this plugin
        this->datasets.append(d);

        // Add to UI combobox for selection of data
        this->form->comboBoxDataset->addItem(d->getName());

        // Add to UI combobox for selection of overlay
        this->form->comboBoxOverlay->addItem(d->getName());

        // Create model info
        createModelInfo(d);
    }

    // Transfer Functions (LUTs)
	else if (kind == "transfer function")
	{
	    // Keep track of transfer function datasets
	    this->tf_datasets.append(d);

        // Add to UI combobox for selection of overlay plane LUT
		this->form->comboBoxOverlayLUT->addItem(d->getName());

		// Add to UI combobox for selection of base layer plane LUT
		this->form->comboBoxBaseLayerLUT->addItem(d->getName());

		// Add to UI combobox for selection of curvature LUT
		this->form->comboBoxCurvatureLUT->addItem(d->getName());

        // Create the lookup table from the transfer function
        createLookupTable(d);
	}

	// Add fibers
	else if (kind == "fibers")
	{
	    // Check if fiber has polydata
	    if (d->getVtkPolyData() == NULL)
			return;

        // Create new fiber struct
        SortedFibers* sortedFibers = new SortedFibers;

        // Initialize struct
        sortedFibers->ds = d;
		sortedFibers->userSelectedLine = 0;

        // Add the new data set to the list of currently available fiber sets
        this->sortedFibersList.append(sortedFibers);

        // Add to UI combobox for distance measurements to fibers
        this->form->comboBoxFiberData->addItem(d->getName());
	}
}

//------------------------[ Dataset changed ]-----------------------\\

void IsosurfaceVisualization::dataSetChanged(data::DataSet * d)
{
    // Assert the data set pointer (should never be NULL)
    Q_ASSERT(d);

	// Get the kind of the data set
    QString kind = d->getKind();

    // Scalar volume
    if(kind == "scalar volume")
    {
        // Get index
        int dsIndex = this->datasets.indexOf(d);

        // Change names in comboboxes
        this->form->comboBoxDataset->setItemText(dsIndex,d->getName());
        this->form->comboBoxOverlay->setItemText(dsIndex,d->getName());
    }

    // Transfer Functions (LUTs)
    else if (kind == "transfer function")
	{
	    // Check if the data set has been added to this plugin
        if (!(this->tf_datasets.contains(d)))
            return;

        // Change the name of the data set in the GUI
        int index = this->tf_datasets.indexOf(d);

        // Create new transfer function LUTs
        createLookupTable(d,index);

        // Update the LUTs
        comboBoxOverlayLUTChanged();
        comboBoxBaseLayerLUTChanged();
		comboBoxCurvatureLUTChanged();
	}
}

//------------------------[ Dataset removed ]-----------------------\\

void IsosurfaceVisualization::dataSetRemoved(data::DataSet * d)
{
    // Assert the data set pointer (should never be NULL)
    Q_ASSERT(d);

	// Get the kind of the data set
    QString kind = d->getKind();

    // Scalar volume
    if(kind == "scalar volume")
    {
        // Check if the data set has been added to this plugin
        if (!(this->datasets.contains(d)))
            return;

        // Get index
        int dsIndex = this->datasets.indexOf(d);

        // Remove from UI combobox for selection of data
        this->form->comboBoxDataset->removeItem(dsIndex);

        // Remove from UI combobox for selection of overlay
        this->form->comboBoxOverlay->removeItem(dsIndex);

        // Remove model info
        this->modelInfoList.removeAt(dsIndex);
    }

    // Transfer Functions (LUTs)
    else if (kind == "transfer function")
	{

	}

	// Fibers
	else if (kind == "fibers")
	{
	    // Check if the data set exists
		int dsIndex = this->findInputDataSet(d);

        // Does not exist, return
		if (dsIndex == -1)
			return;

        // Remove from UI combobox for selection of overlay
        this->form->comboBoxFiberData->removeItem(dsIndex);

        // Remove from collection
        this->sortedFibersList.removeAt(dsIndex);
	}
}

int IsosurfaceVisualization::findInputDataSet(data::DataSet * ds)
{
	int index = 0;

	// Loop through all input fiber data sets
	for (QList<SortedFibers*>::iterator i = this->sortedFibersList.begin(); i != this->sortedFibersList.end(); ++i, ++index)
	{
		// Return the index if we've found the target data set
		if ((*i)->ds == ds)
			return index;
	}

	return -1;
}

///
///      PROCESSING
///

//------------------------[ Create model info ]-----------------------\\

void IsosurfaceVisualization::createModelInfo(data::DataSet * d)
{
    // Create model info struct for new dataset
    ModelInfo* modelInfo = new ModelInfo;
    modelInfo->ds = d;
    modelInfo->alpha = DEFAULT_ALPHA;
    modelInfo->prop = NULL; // null pointer
    modelInfo->polydata = NULL;
    modelInfo->colorString = "0xFFFFFF";
    modelInfo->color[0] = modelInfo->color[1] = modelInfo->color[2] = 1.0;  // corresponding to white: #FFFFFF
    modelInfo->renderStyle = 0;
    modelInfo->reduction = 0.0;
    modelInfo->visible = true;
    modelInfo->specular = DEFAULT_SPECULAR;
    modelInfo->selectLargestComponent = true;

    modelInfo->extractPolyFunc = NULL;

    d->getVtkImageData()->GetScalarRange(modelInfo->scalarRange);
    modelInfo->maximumThreshold = ceil(modelInfo->scalarRange[1]*0.75); // automatic upper range
    modelInfo->minimumThreshold = floor(modelInfo->scalarRange[1]*0.25);    // automatic lower range
    if(modelInfo->minimumThreshold == 0)    // zero lower threshold makes model invisible
        modelInfo->minimumThreshold += 0.01;
    modelInfo->smoothing = ceil(modelInfo->scalarRange[1]) == 1.0 ? 0.0 : DEFAULT_SMOOTHING;    // turn smoothing off for binary image

    d->getVtkImageData()->GetBounds(modelInfo->bounds);
    modelInfo->clippingValues[0] = modelInfo->bounds[0];
    modelInfo->clippingValues[1] = modelInfo->bounds[2];
    modelInfo->clippingValues[2] = modelInfo->bounds[4];

    modelInfo->clippingCube = vtkBox::New();
    modelInfo->clippingCube->SetBounds(modelInfo->bounds);
    modelInfo->invertClipping = false;

    modelInfo->planesActivated[0] = modelInfo->planesActivated[1] = modelInfo->planesActivated[2] = false;
    modelInfo->flipped[0] = modelInfo->flipped[1] = modelInfo->flipped[2] = false;
    modelInfo->dsOverlay = NULL;
    modelInfo->overlayIndex = 0;

    modelInfo->bDsPolyAdded = false;

    modelInfo->alignPlanesToPick = true;

	modelInfo->curvatureLUTIndex = 0;
	modelInfo->curvatureType = 0;

    // create lookup table for grey values
    vtkLookupTable* defaultLUT = vtkLookupTable::New();
    defaultLUT->SetValueRange(0.0, 1.0);
    defaultLUT->SetSaturationRange(0.0, 0.0);
    defaultLUT->SetAlphaRange(0.0, 1.0);
    double range[2];
    d->getVtkImageData()->GetScalarRange(range);
    defaultLUT->SetTableRange(range);
    defaultLUT->Build();

    modelInfo->defaultLUT = defaultLUT;

    // Create the three orthogonal planes for 2D image merging
    double normals[3][3] = {{1,0,0},{0,1,0},{0,0,1}};
    for (int i = 0; i < 3; i++)
    {
        vtkPlane* slicePlane = vtkPlane::New();
        slicePlane->SetNormal(normals[i]);
        slicePlane->SetOrigin(0,0,0);

        vtkImageResliceMapper* imageMapper = vtkImageResliceMapper::New();
        imageMapper->SetInput(d->getVtkImageData());
        imageMapper->SetSlicePlane(slicePlane);
        imageMapper->SetSeparateWindowLevelOperation(0);
        imageMapper->SliceFacesCameraOff();
        imageMapper->SliceAtFocalPointOff();

        vtkImageSlice *image = vtkImageSlice::New();
        image->SetMapper(imageMapper);
        image->Update();
        imageMapper->Delete();
        image->GetProperty()->SetLookupTable(defaultLUT);
        image->GetProperty()->SetUseLookupTableScalarRange(1);

        // Use the identity matrix by default
        vtkMatrix4x4 * id = vtkMatrix4x4::New();
        id->Identity();
        image->PokeMatrix(id);

        vtkObject* tfm;
        if (d->getAttributes()->getAttribute("transformation matrix", tfm ))
        {
            vtkMatrix4x4* transformationMatrix = vtkMatrix4x4::SafeDownCast(tfm);
            if (transformationMatrix == 0)
            {
                this->core()->out()->logMessage("not a valid transformation matrix");
                return;
            }
            image->PokeMatrix(transformationMatrix);
        }
        else
        {
            double* center;
            center = d->getVtkImageData()->GetCenter();
            image->SetOrigin(center[0],center[1],center[2]);
        }

        // add model to list
        modelInfo->orthogonalPlanesModels.append(image);

        // add plane to list
        modelInfo->orthogonalPlanes.append(slicePlane);

        // by default invisible
        image->VisibilityOff();

        double bounds[6];
        image->GetBounds(bounds);
        modelInfo->global_bounds[i*2] = bounds[i*2];
        modelInfo->global_bounds[i*2+1] = bounds[i*2+1];

        this->assembly->AddPart(image);
        image->Delete();
    }

    // Add to local model list
    this->modelInfoList.append(modelInfo);
}

//------------------------[ Process transfer function ]-----------------------\\

void IsosurfaceVisualization::createLookupTable(data::DataSet * d, int index)
{
	 // get color transfer function
    vtkColorTransferFunction* ctf = vtkColorTransferFunction::SafeDownCast(d->getVtkObject());
    int numberOfColorNodes = ctf->GetSize();
    //std::cout << "Number of nodes: " << numberOfColorNodes << std::endl;

    // color range
    double colorRange[2];
    ctf->GetRange(colorRange);
	int numberOfValues = ceil(colorRange[1]);

	// get piecewise function attribute
    vtkObject* cpf;
	vtkPiecewiseFunction* pf;
	int numberOfOpacityNodes;
    if(d->getAttributes()->getAttribute("piecewise function", cpf ))
	{
		// get opacity transfer function
		pf = vtkPiecewiseFunction::SafeDownCast(cpf);
		numberOfOpacityNodes = pf->GetSize();
		//pf->Print(std::cout);
		//std::cout << "Number of nodes: " << numberOfOpacityNodes << std::endl;
	}

	// no opacity transfer function available, create one (no opacity: 1.0 to 1.0)
	else
	{
		pf = vtkPiecewiseFunction::New();
		pf->AdjustRange(colorRange);

		double point1[4] = {0.0, 1.0, 0.0, 1.0};
		pf->SetNodeValue(0,point1);

		double point2[4] = {colorRange[1], 1.0, 0.0, 1.0};
		pf->SetNodeValue(0,point2);
	}

	// opacity range
    double opacityRange[2];
    pf->GetRange(opacityRange);

    // create new lookup table object
    vtkLookupTable* convLUT = vtkLookupTable::New();
    convLUT->SetScaleToLinear();
    convLUT->SetValueRange(0.0, 92.0);
    convLUT->SetNumberOfTableValues(numberOfValues);
    convLUT->SetTableRange(0,numberOfValues);
    convLUT->Build();

    // store the data from the color and opacity transfer functions in qlists
    QList<double> colorIndices;
    QList<double> colorR;
    QList<double> colorG;
    QList<double> colorB;
    QList<double> alpha;

    for(int j = 0; j<numberOfColorNodes; j++)
    {
        double colorVal[6];
        ctf->GetNodeValue(j,colorVal);

        double opacityVal[4];
        pf->GetNodeValue(j,opacityVal);

        colorIndices.append(colorVal[0]);
        colorR.append(colorVal[1]);
        colorG.append(colorVal[2]);
        colorB.append(colorVal[3]);

        alpha.append(opacityVal[1]);
    }

    // fill the lookup table with data. interpolate points between nodes of the transfer function.
    for(int j = 0; j<numberOfValues; j++)
    {
        // linear interpolation
        int ind = -1;
        for(int k = 0; k<colorIndices.length(); k++)
        {
            if(j < colorIndices.at(k))
            {
                ind = std::max(0,k-1);
                break;
            }
        }

        if(ind == -1)
            continue;

        double cl_ind = colorIndices.at(ind);
        double cl_r = colorR.at(ind);
        double cl_g = colorG.at(ind);
        double cl_b = colorB.at(ind);
        double al_v = alpha.at(ind);

        ind = ind + 1;
        double cu_ind = colorIndices.at(ind);
        double cu_r = colorR.at(ind);
        double cu_g = colorG.at(ind);
        double cu_b = colorB.at(ind);
        double au_v = alpha.at(ind);

        double ind_frac = j-cl_ind;
        ind_frac /= (double)(cu_ind - cl_ind);

        convLUT->SetTableValue(j,
                               cl_r*(1-ind_frac) + cu_r*(ind_frac),
                               cl_g*(1-ind_frac) + cu_g*(ind_frac),
                               cl_b*(1-ind_frac) + cu_b*(ind_frac),
                               al_v*(1-ind_frac) + au_v*(ind_frac));
    }

    // add to lookup tables, or replace if it already existed
    if(index == -1)
        lookUpTables.append(convLUT);
    else
        lookUpTables.replace(index,convLUT);

    // save color transfer function
    if(index == -1)
        colorTransferFunctions.append(ctf);
    else
        colorTransferFunctions.replace(index,ctf);
}

//------------------------[ Create or update the mesh ]-----------------------\\

void IsosurfaceVisualization::updateRenderingModels()
{
    // recalculate model
    if(bModelDirty)
    {
        // remove old actor (make conditional)
        if(current_modelInfo->prop != NULL)
        {
            this->assembly->RemovePart(current_modelInfo->prop);
            current_modelInfo->prop = NULL;
        }

        data::DataSet* ds = this->current_modelInfo->ds;

		// determine if we should compute curvature
		bool useCurvature = current_modelInfo->curvatureType > 0;

        // smooth dataset
        VTK_CREATE(vtkImageGaussianSmooth, smoother);
        smoother->SetInput(ds->getVtkImageData());
        smoother->SetStandardDeviation(current_modelInfo->smoothing,current_modelInfo->smoothing,current_modelInfo->smoothing);

        // marching cubes
        VTK_CREATE(vtkMarchingCubes, mcubes);
        mcubes->SetInputConnection(smoother->GetOutputPort());
        mcubes->SetNumberOfContours(1);
		if(useCurvature)
			mcubes->SetComputeScalars(1);
		else
			mcubes->SetComputeScalars(0);
		//mcubes->SetComputeNormals(1);
        mcubes->GenerateValues(1,current_modelInfo->minimumThreshold,current_modelInfo->maximumThreshold);

		// marching cubes progress bar
		this->core()->out()->createProgressBarForAlgorithm(mcubes, "Isosurface Visualization", "Creating isosurface mesh...");
		mcubes->Update();
		this->core()->out()->deleteProgressBarForAlgorithm(mcubes);

        // select largest component
        VTK_CREATE(vtkPolyDataConnectivityFilter, connectivity);
        if(current_modelInfo->selectLargestComponent)
        {
            connectivity->SetInputConnection(mcubes->GetOutputPort());
            connectivity->SetExtractionModeToLargestRegion();
        }

        /*vtkSmoothPolyDataFilter* smoother2 = vtkSmoothPolyDataFilter::New();
        smoother2->SetInputConnection(connectivity->GetOutputPort());
        smoother2->SetNumberOfIterations(50);
        smoother2->Update();*/

        // decimate
        VTK_CREATE(vtkDecimatePro, deci);
        if(current_modelInfo->reduction > 0.1)
        {
            if(current_modelInfo->selectLargestComponent)
                deci->SetInputConnection(connectivity->GetOutputPort());
            else
                deci->SetInputConnection(mcubes->GetOutputPort());
            deci->SetTargetReduction(current_modelInfo->reduction);
            deci->PreserveTopologyOn();
            deci->AccumulateErrorOn();

			// decimator progress bar
			this->core()->out()->createProgressBarForAlgorithm(deci, "Isosurface Visualization", "Decimating the mesh...");
			deci->Update();
			this->core()->out()->deleteProgressBarForAlgorithm(deci);
        }

        // Clip geometry (based on slider values)
        vtkExtractPolyDataGeometry* clipper = vtkExtractPolyDataGeometry::New();
        clipper->SetImplicitFunction(current_modelInfo->clippingCube);
        if(current_modelInfo->reduction > 0.1)
            clipper->SetInputConnection(deci->GetOutputPort());
        else
		{
			if(current_modelInfo->selectLargestComponent)
                clipper->SetInputConnection(connectivity->GetOutputPort());
            else
                clipper->SetInputConnection(mcubes->GetOutputPort());
		}
        clipper->ExtractInsideOn();
        current_modelInfo->extractPolyFunc = clipper;

		// Create polydata mapper
		vtkPolyDataMapper* polyMapper = vtkPolyDataMapper::New();
		//polyMapper->ImmediateModeRenderingOn();

		// curvature
		if(useCurvature)
		{
			// perform curvature filter
			vtkSmartPointer<vtkCurvaturesShapeIndex> curvaturesFilter =
			vtkSmartPointer<vtkCurvaturesShapeIndex>::New();
			curvaturesFilter->SetInputConnection(clipper->GetOutputPort());

			switch(current_modelInfo->curvatureType)
			{
			    case 1:
                    curvaturesFilter->SetCurvatureTypeToMaximum();
                    break;

                case 2:
                    curvaturesFilter->SetCurvatureTypeToMinimum();
                    break;

                case 3:
                    curvaturesFilter->SetCurvatureTypeToGaussian();
                    break;

                case 4:
                    curvaturesFilter->SetCurvatureTypeToMean();
                    break;

				case 5:
					curvaturesFilter->SetCurvatureTypeToShapeIndex();
                    break;

				case 6:
                    curvaturesFilter->SetCurvatureTypeToCurvedness();
                    break;
			}

			// progress of curvature
			this->core()->out()->createProgressBarForAlgorithm(curvaturesFilter, "Isosurface Visualization", "Curvature filter...");
			curvaturesFilter->Update();
			this->core()->out()->deleteProgressBarForAlgorithm(curvaturesFilter);

            // Set output of curvature
            polyMapper->SetInputConnection(curvaturesFilter->GetOutputPort());

            // Get scalar range
            double scalarRange[2];
            curvaturesFilter->GetOutput()->GetScalarRange(scalarRange);
            std::cout << "Scalar range: " << scalarRange[0] << " - " << scalarRange[1] << std::endl;

			// create default curvature lookup table
            if(this->current_modelInfo->curvatureLUTIndex == 0)
            {
                /*vtkLookupTable* convLUT;
				convLUT = vtkLookupTable::New();
				convLUT->SetNumberOfTableValues(100);
				convLUT->SetRange(0,100);
				convLUT->Build();

				for(int j = 0; j<100; j++)
				{
					double* bbb = convLUT->GetTableValue(j);

					if( j > 0 && j < 80)
						convLUT->SetTableValue(j,1.0,1.0,1.0,1.0);
					else if(j >= 80 && j < 100)
						convLUT->SetTableValue(j,1-(j-80)/20.0*0.25,1-(j-80)/20.0*0.25,1-(j-80)/20.0*0.25,1.0);
					else
						convLUT->SetTableValue(j,0.75,0.75,0.75,1.0);
				}
				convLUT->SetRampToSCurve();
				polyMapper->SetLookupTable(convLUT);*/

				double mid = 0.5 * (scalarRange[1] + scalarRange[0]);

                vtkColorTransferFunction* convLUT = vtkColorTransferFunction::New();
				convLUT->AddRGBPoint(scalarRange[0], 1.0, 0.0, 0.0);
				convLUT->AddRGBPoint(mid, 0.0, 1.0, 0.0);
				convLUT->AddRGBPoint(scalarRange[1], 0.0, 0.0, 1.0);
				convLUT->Print(std::cout);

                polyMapper->SetLookupTable(convLUT);
			}

			// Use custom transfer function
			else
			{
			    // Get LUT from saved list
			    //vtkLookupTable* convLUT;
				//convLUT = this->lookUpTables.at(current_modelInfo->curvatureLUTIndex - 1);

				vtkColorTransferFunction* convLUT;
				convLUT = this->colorTransferFunctions.at(current_modelInfo->curvatureLUTIndex - 1);

				// Use curvature output in polydatamapper and set LUT
				polyMapper->SetLookupTable(convLUT);
			}

			current_modelInfo->usingCurvature = true;

			if(scalarBar != NULL)
            {
                this->assembly->RemovePart(scalarBar);
                scalarBar = NULL;
            }
            scalarBar = vtkSmartPointer<vtkScalarBarActor>::New();
            scalarBar->SetLookupTable(polyMapper->GetLookupTable());
            scalarBar->SetTitle(
            curvaturesFilter->GetOutput()->GetPointData()->GetScalars()->GetName());
            scalarBar->SetNumberOfLabels(5);
            this->assembly->AddPart(scalarBar);
		}
		else
		{
			polyMapper->SetInputConnection(clipper->GetOutputPort());
			current_modelInfo->usingCurvature = false;
		}

		// Create actor
        VTK_CREATE(vtkLODActor,actor);
        actor->SetMapper(polyMapper);
		// No curvature, then set a color
		if(!useCurvature)
		{
			actor->GetProperty()->SetColor(current_modelInfo->color);
		}
        if(current_modelInfo->alpha < 1.0)
        {
            actor->GetProperty()->SetOpacity(current_modelInfo->alpha);
        }
        actor->GetProperty()->SetBackfaceCulling(1);
        actor->GetProperty()->SetInterpolationToPhong();
        actor->GetProperty()->SetSpecular(current_modelInfo->specular);
        actor->GetProperty()->SetSpecularPower(2.0);

        actor->SetOrientation(0,0,0);
        actor->SetPosition(0,0,0);
        actor->SetScale(1,1,1);
        actor->SetOrigin(0,0,0);

		vtkProperty *prop = actor->GetProperty();

		//vtkXMLMaterial* material = vtkXMLMaterial::CreateInstance("D:/Research/Software/vISTe/source/plugins/testing/smeesters/IsosurfaceVisualization/GLSLTwisted.xml");
		//material->Print(std::cout);

		/*glEnable(0x3000); // GL_CLIP_DISTANCE0
		glEnable(0x3000+1); // GL_CLIP_DISTANCE1
		glEnable(0x3000+2); // GL_CLIP_DISTANCE2
		//prop->LoadMaterial(material);

		//prop->LoadMaterial("D:/Research/Software/vISTe/source/plugins/testing/smeesters/IsosurfaceVisualization/ClippingPlane.xml");
		prop->LoadMaterial("/home/stephan/Research/Software/vISTe/src/smeesters/plugins/testing/smeesters/IsosurfaceVisualization/ClippingPlane.xml");
        //prop->LoadMaterial("/home/linux/Research/fMRIDTI_Visual/Software/vISTe/subversion/smeesters/plugins/testing/smeesters/IsosurfaceVisualization/ClippingPlane.xml");

		prop->ShadingOn();
		prop->AddShaderVariable("ClipX",0.0);
		prop->AddShaderVariable("ClipY",0.0);
		prop->AddShaderVariable("ClipZ",0.0);*/


        // Use the identity matrix by default
        vtkMatrix4x4 * id = vtkMatrix4x4::New();
        id->Identity();
        actor->PokeMatrix(id);

        vtkObject* tfm;
        if (ds->getAttributes()->getAttribute("transformation matrix", tfm ))
        {
            vtkMatrix4x4* transformationMatrix = vtkMatrix4x4::SafeDownCast(tfm);
            if (transformationMatrix == 0)
            {
                this->core()->out()->logMessage("not a valid transformation matrix");
                return;
            }
            actor->PokeMatrix(transformationMatrix);

            //current_modelInfo->ds_poly->getAttributes()->addAttribute("transformation matrix", vtkObject::SafeDownCast(transformationMatrix));
        }
        else
        {
            double* center;
            center = actor->GetCenter();
            actor->SetOrigin(center[0],center[1],center[2]);
        }

        id->Delete();
        tfm->Delete();

		current_modelInfo->polyDataMapper = polyMapper;
        current_modelInfo->polydata = connectivity->GetOutput();
        current_modelInfo->prop = actor;
        this->assembly->AddPart(actor);

		for(int i =0; i<3; i++)
        {
            this->updateClippingPlaneSlider(i,current_modelInfo->clippingValues[i]);
        }
        // Reset the camera of the 3D volume
        //this->fullCore()->canvas()->GetRenderer3D()->ResetCamera();

        // Reset clipping rnage
       // this->fullCore()->canvas()->GetRenderer3D()->ResetCameraClippingRange();

//        // Update the polydata of the DataSet or the manager
//        current_modelInfo->ds_poly->updateData(clipper->GetOutput());
//
//        // Add the model to the Data manager if not already done
//        if(!(current_modelInfo->bDsPolyAdded))
//        {
//            this->core()->data()->addDataSet(current_modelInfo->ds_poly);
//            std::cout << "ADD DATA!";
//            current_modelInfo->bDsPolyAdded = true;
//        }
//
//        // else Signal the core that the data set has changed
//        else
//        {
//            this->core()->data()->dataSetChanged(current_modelInfo->ds_poly);
//        }



    }

    // adjust existing model
    else
    {
        current_modelInfo->prop->GetProperty()->SetColor(current_modelInfo->color);
        current_modelInfo->prop->GetProperty()->SetOpacity(current_modelInfo->alpha);
        current_modelInfo->prop->GetProperty()->SetSpecular(current_modelInfo->specular);
    }

    bModelDirty = false;
    this->core()->render();

    //vtkDebugLeaks::PrintCurrentLeaks();
}

//------------------------[ Create text labels ]-----------------------\\

vtkActor2D* IsosurfaceVisualization::GenerateLabels(vtkSmartPointer<vtkPoints> points, vtkSmartPointer<vtkStringArray> labels)
{
	VTK_CREATE(vtkPolyData, polydata);
	polydata->SetPoints(points);

	VTK_CREATE(vtkVertexGlyphFilter, glyphFilter);
	glyphFilter->SetInputConnection(polydata->GetProducerPort());
	glyphFilter->Update();

	// Add label array.
	glyphFilter->GetOutput()->GetPointData()->AddArray(labels);

	// Create a mapper and actor for the points.
	VTK_CREATE(vtkPolyDataMapper,pointMapper);
	pointMapper->SetInputConnection(glyphFilter->GetOutputPort());

	VTK_CREATE(vtkActor, pointActor);
	pointActor->SetMapper(pointMapper);

	// Generate the label hierarchy.
	VTK_CREATE(vtkPointSetToLabelHierarchy, pointSetToLabelHierarchyFilter);
	pointSetToLabelHierarchyFilter->SetInputConnection(
	glyphFilter->GetOutputPort());
	pointSetToLabelHierarchyFilter->SetLabelArrayName("labels");
	pointSetToLabelHierarchyFilter->Update();

	// Create a mapper and actor for the labels.
	VTK_CREATE(vtkLabelPlacementMapper, labelMapper);
	labelMapper->SetInputConnection(pointSetToLabelHierarchyFilter->GetOutputPort());
	labelMapper->UseDepthBufferOff();
	labelMapper->PlaceAllLabelsOn();

	vtkActor2D* labelActor = vtkActor2D::New();
	labelActor->SetMapper(labelMapper);
	labelActor->SetLayerNumber(5);

	return labelActor;
}

///
///     CALLBACKS
///

//------------[ Callback of clipping planes position from picker ]----------------\\

void IsosurfaceVisualization::setClippingPlanesPosition(double* pos)
{
    if(this->current_modelInfo != NULL)
    {
        double local_pos[3];
        double frac_pos[3];

        for(int i = 0; i<3; i++)
        {

            double n_pos = (double)(pos[i]-current_modelInfo->global_bounds[i*2])/(abs(current_modelInfo->global_bounds[i*2] - current_modelInfo->global_bounds[i*2+1]));
            frac_pos[i] = n_pos;
            n_pos *= current_modelInfo->bounds[i*2+1];
            local_pos[i] = n_pos;
            this->clickedPoint[i] = pos[i];

            std::cout << i << " " << pos[i] << " --- " << n_pos << std::endl;
        }

        //this->clickedPoint = local_pos;

        QString sliceActorNames[3];
        sliceActorNames[0] = "X Plane";
        sliceActorNames[1] = "Y Plane";
        sliceActorNames[2] = "Z Plane";

        for(int axis = 0; axis < 3; ++axis)
        {
            std::cout << "AXIS: " << axis << std::endl;
            vtkImageSliceActor * sliceActor = static_cast<vtkImageSliceActor*>(this->core()->data()->getDataSet(sliceActorNames[axis], "sliceActor")->getVtkObject());
            //sliceActor->Print(std::cout);
            //sliceActor->GetUserMatrix()->Print(std::cout);

            vtkActor* pointer = this->pointer2DList.at(axis);
            pointer->SetUserMatrix(sliceActor->GetUserMatrix());
        }

        //vtkMedicalCanvas * canvas = this->fullCore()->canvas();

        //vtkImageSliceActor * sliceActor = this->fullCore()->data()->getDataSet()

        for(int i = 0; i<3; i++)
        {
            vtkActor* pointer = this->pointer2DList.at(i);

            switch(i)
            {
                case 0: // sagittal
                    pointer->SetPosition(frac_pos[0]*181,frac_pos[1]*217,frac_pos[2]*181);
                    break;

                case 1: // coronal
                    pointer->SetPosition(frac_pos[0]*181,300,frac_pos[2]*181);
                    break;

                case 2: // axial
                    pointer->SetPosition(frac_pos[0]*181,frac_pos[1]*217,-90);
                    break;
            }
        }

        // set medical canvas positions
        this->current_modelInfo->ds->getAttributes()->addAttribute("SlicePosX", (int)(frac_pos[0]*181));
        this->current_modelInfo->ds->getAttributes()->addAttribute("SlicePosY", (int)(frac_pos[1]*217));
        this->current_modelInfo->ds->getAttributes()->addAttribute("SlicePosZ", (int)(frac_pos[2]*181));

        if(!current_modelInfo->alignPlanesToPick)
            return;

        vtkPoints* points = vtkPoints::New();
        int radius = 5;
        int numPoints = 0;
        for(int i = -1; i<=1; i++)
        {
            for(int j = -1; j<=1; j++)
            {
                for(int k = -1; k<=1; k++)
                {
                    double neighbour_pos[3];
                    neighbour_pos[0] = local_pos[0] + i*radius;
                    neighbour_pos[1] = local_pos[1] + j*radius;
                    neighbour_pos[2] = local_pos[2] + k*radius;
                    points->InsertNextPoint(neighbour_pos);

                    //std::cout << "npoint: " << neighbour_pos[0] << " " << neighbour_pos[1] << " " << neighbour_pos[2] << std::endl;

                    numPoints++;
                }
            }
        }

        vtkPolyData* pointsPolyData = vtkPolyData::New();
        pointsPolyData->SetPoints(points);

        vtkSelectEnclosedPoints* selectEnclosedPoints = vtkSelectEnclosedPoints::New();
        selectEnclosedPoints->SetInput(pointsPolyData);
        selectEnclosedPoints->SetSurface(current_modelInfo->polydata);
        selectEnclosedPoints->Update();

        bool isInside = false;
        for(int i = 0; i<numPoints; i++)
        {
            if(selectEnclosedPoints->IsInside(i))
                isInside = true;
        }
        std::cout << "inside: " << isInside << std::endl;


        if(isInside)
        {
            for(int i = 0; i<3; i++)
            {
                if(!current_modelInfo->planesActivated[i])
                    continue;

                this->updateClippingPlaneSlider(i,(int)local_pos[i],false);

                switch(i)
                {
                    case 0:
                        this->form->horizontalSliderX->setValue(local_pos[i]);
                        break;

                    case 1:
                        this->form->horizontalSliderY->setValue(local_pos[i]);
                        break;

                    case 2:
                        this->form->horizontalSliderZ->setValue(local_pos[i]);
                        break;
                }
            }
        }

        this->core()->data()->dataSetChanged(this->current_modelInfo->ds);




        //this->core()->render();

        /*std::cout << "xrange: " << current_modelInfo->global_bounds[0] << " -- " << current_modelInfo->global_bounds[1] << std::endl;
        std::cout << "yrange: " << current_modelInfo->global_bounds[2] << " -- " << current_modelInfo->global_bounds[3] << std::endl;
        std::cout << "zrange: " << current_modelInfo->global_bounds[4] << " -- " << current_modelInfo->global_bounds[5] << std::endl;*/
    }


}

///
///     UPDATE FROM USER INPUT
///

//-----------[ Update clipping plane slider effects on mesh ]-----------------------\\

void IsosurfaceVisualization::updateClippingPlaneSlider(int direction, int value, bool render)
{
    // return if no model is chosen
    if(current_modelInfo == NULL)
        return;

    //std::cout << "Update clipping plane slider: " << direction << " value: " << value << std::endl;

    // adjust bounds of clipping cube used to cut away the mesh
    double cubeBounds[6];
    current_modelInfo->clippingCube->GetBounds(cubeBounds);

    // flipped, reverse bounds
    if(current_modelInfo->flipped[direction])
    {
        cubeBounds[direction*2] = current_modelInfo->bounds[direction*2]; // minimum bound for direction
        cubeBounds[direction*2+1] = value;
    }

    // not flipped
    else
    {
        cubeBounds[direction*2] = value;
        cubeBounds[direction*2+1] = current_modelInfo->bounds[direction*2+1]; // maximum bound for direction
    }

    current_modelInfo->clippingCube->SetBounds(cubeBounds);
    current_modelInfo->clippingValues[direction] = value;

	/*if(direction == 0)
		this->current_modelInfo->prop->GetProperty()->AddShaderVariable("ClipX",(double)(value) );
	else if(direction == 1)
		this->current_modelInfo->prop->GetProperty()->AddShaderVariable("ClipY",(double)(value) );
	else if(direction == 2)
		this->current_modelInfo->prop->GetProperty()->AddShaderVariable("ClipZ",(double)(value) );*/

    // compute origin of orthogonal clipping planes in order to move it
    double origin[3];
    origin[direction] = (value - this->current_modelInfo->bounds[direction*2+1]/2.0);

    // temporary origin correction
    if(direction == 1) { origin[direction] -= 18; }
    if(direction == 2) { origin[direction] += 18; }

    //if(direction == 0) { origin[direction] -= 6; }
    //if(direction == 1) { origin[direction] += 10; }
    //if(direction == 2) { origin[direction] -= 4; }


    vtkPlane* plane = this->current_modelInfo->orthogonalPlanes.at(direction);
    plane->SetOrigin(origin);

    // overlay plane visibility and position
    if(this->overlay_modelInfo != NULL)
    {
        vtkPlane* plane_overlay = this->overlay_modelInfo->orthogonalPlanes.at(direction);
        if(direction == 0)
            origin[direction] -= -(current_modelInfo->flipped[direction]*2 - 1)*0.1;
        else if(direction == 1)
            origin[direction] -= -(current_modelInfo->flipped[direction]*2 - 1)*0.1;
        else if(direction == 2)
            origin[direction] += -(current_modelInfo->flipped[direction]*2 - 1)*0.1;
        plane_overlay->SetOrigin(origin);
    }

    if(render)
        this->core()->render();
}

//-----------[ Update clipping plane enabled checkbox action ]-----------------------\\

void IsosurfaceVisualization::updateClippingPlaneEnabled(int direction, bool checked)
{
    current_modelInfo->planesActivated[direction] = checked;

    // orthogonal plane visibility and position
    vtkImageSlice* planeModel = this->current_modelInfo->orthogonalPlanesModels.at(direction);
    planeModel->SetVisibility(checked);

    // overlay plane visibility and position
    if(this->overlay_modelInfo != NULL)
    {
        vtkImageSlice* planeModel_overlay = this->overlay_modelInfo->orthogonalPlanesModels.at(direction);
        planeModel_overlay->SetVisibility(checked);
    }

    this->updateClippingPlaneSlider(direction,current_modelInfo->clippingValues[direction]);
}


///
///      GUI CONTROLS
///

//------------------------[ Connect Qt elements ]-----------------------\\

void IsosurfaceVisualization::connectAll()
{
    // Select data set
    connect(this->form->comboBoxDataset,SIGNAL(currentIndexChanged(int)),this,SLOT(comboBoxDataChanged()));

    // Generation
    connect(this->form->checkBoxVisible,SIGNAL(toggled(bool)),this,SLOT(checkBoxVisibleChanged(bool)));
    connect(this->form->comboBoxStyle,SIGNAL(currentIndexChanged(int)),this,SLOT(comboBoxStyleChanged()));
    connect(this->form->inputSmoothing,SIGNAL(valueChanged(double)),this,SLOT(inputSmoothingChanged(double)));
    connect(this->form->inputReduction,SIGNAL(valueChanged(double)),this,SLOT(inputReductionChanged(double)));
    connect(this->form->inputMinimumThreshold,SIGNAL(valueChanged(double)),this,SLOT(inputMinimumThresholdChanged(double)));
    connect(this->form->inputMaximumThreshold,SIGNAL(valueChanged(double)),this,SLOT(inputMaximumThresholdChanged(double)));
    connect(this->form->inputColor,SIGNAL(clicked()),this,SLOT(inputColorChanged()));
    connect(this->form->inputSpecular,SIGNAL(valueChanged(double)),this,SLOT(inputSpecularChanged(double)));
    connect(this->form->inputAlpha,SIGNAL(valueChanged(double)),this,SLOT(inputAlphaChanged(double)));
    connect(this->form->checkBoxLargestComponent,SIGNAL(toggled(bool)),this,SLOT(checkBoxLargestComponentChanged(bool)));
    connect(this->form->comboBoxCurvatureType,SIGNAL(currentIndexChanged(int)),this,SLOT(comboBoxCurvatureTypeChanged()));
	connect(this->form->comboBoxCurvatureLUT,SIGNAL(currentIndexChanged(int)),this,SLOT(comboBoxCurvatureLUTChanged()));
    connect(this->form->buttonUpdate,SIGNAL(clicked()),this,SLOT(buttonUpdateClicked()));

    // Clipping planes
    connect(this->form->checkBoxX,SIGNAL(toggled(bool)),this,SLOT(checkBoxXChanged(bool)));
    connect(this->form->checkBoxFlipX,SIGNAL(toggled(bool)),this,SLOT(checkBoxFlipXChanged(bool)));
    connect(this->form->spinX,SIGNAL(valueChanged(int)),this,SLOT(spinXChanged(int)));
    connect(this->form->horizontalSliderX,SIGNAL(valueChanged(int)),this,SLOT(horizontalSliderXChanged(int)));

    connect(this->form->checkBoxY,SIGNAL(toggled(bool)),this,SLOT(checkBoxYChanged(bool)));
    connect(this->form->checkBoxFlipY,SIGNAL(toggled(bool)),this,SLOT(checkBoxFlipYChanged(bool)));
    connect(this->form->spinY,SIGNAL(valueChanged(int)),this,SLOT(spinYChanged(int)));
    connect(this->form->horizontalSliderY,SIGNAL(valueChanged(int)),this,SLOT(horizontalSliderYChanged(int)));

    connect(this->form->checkBoxZ,SIGNAL(toggled(bool)),this,SLOT(checkBoxZChanged(bool)));
    connect(this->form->checkBoxFlipZ,SIGNAL(toggled(bool)),this,SLOT(checkBoxFlipZChanged(bool)));
    connect(this->form->spinZ,SIGNAL(valueChanged(int)),this,SLOT(spinZChanged(int)));
    connect(this->form->horizontalSliderZ,SIGNAL(valueChanged(int)),this,SLOT(horizontalSliderZChanged(int)));

    connect(this->form->checkBoxAlignPlanesToPick,SIGNAL(toggled(bool)),this,SLOT(checkBoxAlignPlanesToPickChanged(bool)));
    connect(this->form->checkBoxInvertClipping,SIGNAL(toggled(bool)),this,SLOT(checkBoxInvertClippingChanged(bool)));

    connect(this->form->comboBoxBaseLayerLUT,SIGNAL(currentIndexChanged(int)),this,SLOT(comboBoxBaseLayerLUTChanged()));
    connect(this->form->comboBoxOverlay,SIGNAL(currentIndexChanged(int)),this,SLOT(comboBoxOverlayChanged()));
    connect(this->form->comboBoxOverlayLUT,SIGNAL(currentIndexChanged(int)),this,SLOT(comboBoxOverlayLUTChanged()));

    // Export
    connect(this->form->buttonSaveMesh,SIGNAL(clicked()),this,SLOT(buttonSaveMeshClicked()));

    // Measurement
    connect(this->form->buttonSetPointA,SIGNAL(clicked()),this,SLOT(buttonSetPointAClicked()));
    connect(this->form->lineEditNamePointA,SIGNAL(textChanged(QString)),this,SLOT(lineEditNamePointAChanged(QString)));

    connect(this->form->buttonSetPointB,SIGNAL(clicked()),this,SLOT(buttonSetPointBClicked()));
    connect(this->form->lineEditNamePointB,SIGNAL(textChanged(QString)),this,SLOT(lineEditNamePointBChanged(QString)));
    connect(this->form->comboBoxFiberData,SIGNAL(currentIndexChanged(int)),this,SLOT(comboBoxFiberDataChanged()));
    connect(this->form->spinFiberChoice,SIGNAL(valueChanged(int)),this,SLOT(fiberSelectUpdate(int)));
	connect(this->form->sliderFiberChoice,SIGNAL(valueChanged(int)),this,SLOT(fiberSelectUpdate(int)));
	connect(this->form->spinFiberRefinement,SIGNAL(valueChanged(double)),this,SLOT(fiberRefinementUpdate(double)));
	connect(this->form->sliderFiberRefinement,SIGNAL(valueChanged(int)),this,SLOT(fiberRefinementUpdate(int)));

    connect(this->form->buttonSetLineColor,SIGNAL(clicked()),this,SLOT(buttonSetLineColorClicked()));
    //connect(this->form->buttonSaveMeasurement,SIGNAL(clicked()),this,SLOT(buttonSaveMeasurementClicked()));
}

//------------------------[ Disconnect Qt elements ]-----------------------\\

void IsosurfaceVisualization::disconnectAll()
{
    disconnect(this->form->comboBoxDataset,SIGNAL(currentIndexChanged(int)),this,SLOT(comboBoxDataChanged()));
    disconnect(this->form->checkBoxVisible,SIGNAL(toggled(bool)),this,SLOT(checkBoxVisibleChanged(bool)));
    disconnect(this->form->buttonUpdate,SIGNAL(clicked()),this,SLOT(buttonUpdateClicked()));
    disconnect(this->form->inputMaximumThreshold,SIGNAL(valueChanged(double)),this,SLOT(inputMaximumThresholdChanged(double)));
    disconnect(this->form->inputMinimumThreshold,SIGNAL(valueChanged(double)),this,SLOT(inputMinimumThresholdChanged(double)));
    disconnect(this->form->inputSmoothing,SIGNAL(valueChanged(double)),this,SLOT(inputSmoothingChanged(double)));
    disconnect(this->form->inputColor,SIGNAL(textChanged(QString)),this,SLOT(inputColorChanged(QString)));
    disconnect(this->form->inputAlpha,SIGNAL(valueChanged(double)),this,SLOT(inputAlphaChanged(double)));
    disconnect(this->form->inputReduction,SIGNAL(valueChanged(double)),this,SLOT(inputReductionChanged(double)));
    disconnect(this->form->inputSpecular,SIGNAL(valueChanged(double)),this,SLOT(inputSpecularChanged(double)));
    disconnect(this->form->horizontalSliderX,SIGNAL(valueChanged(int)),this,SLOT(horizontalSliderXChanged(int)));
    disconnect(this->form->horizontalSliderY,SIGNAL(valueChanged(int)),this,SLOT(horizontalSliderYChanged(int)));
    disconnect(this->form->horizontalSliderZ,SIGNAL(valueChanged(int)),this,SLOT(horizontalSliderZChanged(int)));
    disconnect(this->form->checkBoxX,SIGNAL(toggled(bool)),this,SLOT(checkBoxXChanged(bool)));
    disconnect(this->form->checkBoxY,SIGNAL(toggled(bool)),this,SLOT(checkBoxYChanged(bool)));
    disconnect(this->form->checkBoxZ,SIGNAL(toggled(bool)),this,SLOT(checkBoxZChanged(bool)));
    disconnect(this->form->comboBoxStyle,SIGNAL(currentIndexChanged(int)),this,SLOT(comboBoxStyleChanged()));
    disconnect(this->form->comboBoxOverlay,SIGNAL(currentIndexChanged(int)),this,SLOT(comboBoxOverlayChanged()));
    disconnect(this->form->comboBoxOverlayLUT,SIGNAL(currentIndexChanged(int)),this,SLOT(comboBoxOverlayLUTChanged()));
    disconnect(this->form->comboBoxBaseLayerLUT,SIGNAL(currentIndexChanged(int)),this,SLOT(comboBoxBaseLayerLUTChanged()));
    disconnect(this->form->checkBoxFlipX,SIGNAL(toggled(bool)),this,SLOT(checkBoxFlipXChanged(bool)));
    disconnect(this->form->checkBoxFlipY,SIGNAL(toggled(bool)),this,SLOT(checkBoxFlipYChanged(bool)));
    disconnect(this->form->checkBoxFlipZ,SIGNAL(toggled(bool)),this,SLOT(checkBoxFlipZChanged(bool)));
    disconnect(this->form->buttonSaveMesh,SIGNAL(clicked()),this,SLOT(buttonSaveMeshClicked()));
    disconnect(this->form->spinX,SIGNAL(valueChanged(int)),this,SLOT(spinXChanged(int)));
    disconnect(this->form->spinY,SIGNAL(valueChanged(int)),this,SLOT(spinYChanged(int)));
    disconnect(this->form->spinZ,SIGNAL(valueChanged(int)),this,SLOT(spinZChanged(int)));
}

//------------------------[ Select another dataset from combobox ]-----------------------\\

void IsosurfaceVisualization::comboBoxDataChanged()
{
    // select model info matching the dataset
    int index = this->form->comboBoxDataset->currentIndex() - 1;

    // "None" value
    if(index < 0)
    {
        // disable gui
		this->form->toolBox->setEnabled(false);
    }
    else
    {
        // enable gui
		this->form->toolBox->setEnabled(true);

        // select model info struct
        this->current_modelInfo = this->modelInfoList.at(index);

        // update gui values
        this->form->horizontalSliderX->setMinimum(current_modelInfo->bounds[0]);
        this->form->horizontalSliderX->setMaximum(current_modelInfo->bounds[1]);
        this->form->horizontalSliderY->setMinimum(current_modelInfo->bounds[2]);
        this->form->horizontalSliderY->setMaximum(current_modelInfo->bounds[3]);
        this->form->horizontalSliderZ->setMinimum(current_modelInfo->bounds[4]);
        this->form->horizontalSliderZ->setMaximum(current_modelInfo->bounds[5]);
        this->form->horizontalSliderX->setValue(current_modelInfo->clippingValues[0]);
        this->form->horizontalSliderY->setValue(current_modelInfo->clippingValues[1]);
        this->form->horizontalSliderZ->setValue(current_modelInfo->clippingValues[2]);

        this->form->spinX->setMinimum(current_modelInfo->bounds[0]);
        this->form->spinX->setMaximum(current_modelInfo->bounds[1]);
        this->form->spinY->setMinimum(current_modelInfo->bounds[2]);
        this->form->spinY->setMaximum(current_modelInfo->bounds[3]);
        this->form->spinZ->setMinimum(current_modelInfo->bounds[4]);
        this->form->spinZ->setMaximum(current_modelInfo->bounds[5]);
        this->form->spinX->setValue(current_modelInfo->clippingValues[0]);
        this->form->spinY->setValue(current_modelInfo->clippingValues[1]);
        this->form->spinZ->setValue(current_modelInfo->clippingValues[2]);

        this->form->inputAlpha->setValue(current_modelInfo->alpha);
        this->form->inputSmoothing->setValue(current_modelInfo->smoothing);
        this->form->inputSpecular->setValue(current_modelInfo->specular);
        this->form->checkBoxVisible->setChecked(current_modelInfo->visible);
        this->form->inputReduction->setValue(current_modelInfo->reduction);
        this->form->comboBoxStyle->setCurrentIndex(current_modelInfo->renderStyle);
        this->form->checkBoxLargestComponent->setChecked(current_modelInfo->selectLargestComponent);
        this->form->checkBoxInvertClipping->setChecked(current_modelInfo->invertClipping);

        this->form->inputMinimumThreshold->setMinimum(current_modelInfo->scalarRange[0]);
        this->form->inputMaximumThreshold->setMaximum(current_modelInfo->scalarRange[1]);
        this->form->inputMaximumThreshold->setValue(current_modelInfo->maximumThreshold);
        this->form->inputMinimumThreshold->setValue(current_modelInfo->minimumThreshold);

        this->form->checkBoxX->setChecked(current_modelInfo->planesActivated[0]);
        this->form->checkBoxY->setChecked(current_modelInfo->planesActivated[1]);
        this->form->checkBoxZ->setChecked(current_modelInfo->planesActivated[2]);

        this->form->checkBoxFlipX->setChecked(current_modelInfo->flipped[0]);
        this->form->checkBoxFlipY->setChecked(current_modelInfo->flipped[1]);
        this->form->checkBoxFlipZ->setChecked(current_modelInfo->flipped[2]);

        // enable all items in overlay combobox
        int numberOfItems = this->form->comboBoxOverlay->count();
        for(int i = 1; i<numberOfItems; i++)
        {
            QModelIndex ind = this->form->comboBoxOverlay->model()->index(i, 0);
            QVariant v(1 | 32);
            this->form->comboBoxOverlay->model()->setData(ind, v, Qt::UserRole - 1);
        }

        // disable current dataset item in overlay combobox
        QModelIndex ind = this->form->comboBoxOverlay->model()->index(index + 1, 0); // index + 1 to skip "None"
        QVariant v(0);
        this->form->comboBoxOverlay->model()->setData(ind, v, Qt::UserRole - 1);

        // set overlay combobox index
        this->form->comboBoxOverlay->setCurrentIndex(this->current_modelInfo->overlayIndex);
        this->form->comboBoxOverlayLUT->setCurrentIndex(this->current_modelInfo->overlayLUTIndex);

        // set base layer combobox index
        this->form->comboBoxBaseLayer->setCurrentIndex(index);
        this->form->comboBoxBaseLayerLUT->setCurrentIndex(this->current_modelInfo->baseLayerLUTIndex);

        // mark model update required
		bModelDirty = true;
    }
}

//---------------[ Select overlay dataset for clipping planes ]-----------------\\

void IsosurfaceVisualization::comboBoxOverlayChanged()
{
    // select model info matching the dataset
    int index = this->form->comboBoxOverlay->currentIndex();

    this->current_modelInfo->overlayIndex = index;

    if((index - 1) >= 0) // Skip "None" value
    {
        this->overlay_modelInfo = this->modelInfoList.at(index - 1);

        this->updateClippingPlaneEnabled(0,current_modelInfo->planesActivated[0]);
        this->updateClippingPlaneEnabled(1,current_modelInfo->planesActivated[1]);
        this->updateClippingPlaneEnabled(2,current_modelInfo->planesActivated[2]);
    }
    else
    {
        this->overlay_modelInfo = NULL;
    }

    this->core()->render();
}

//---------------[ Select overlay LUT for clipping planes ]-----------------\\

void IsosurfaceVisualization::comboBoxOverlayLUTChanged()
{
    // select model info matching the dataset
    int index = this->form->comboBoxOverlayLUT->currentIndex();

    if(this->overlay_modelInfo != NULL)
    {
        this->current_modelInfo->overlayLUTIndex = index;

        if(index > 0)
        {
            vtkLookupTable* ctf = this->lookUpTables.at(index - 1);
            for(int i = 0; i<3; i++)
            {
                this->overlay_modelInfo->orthogonalPlanesModels[i]->GetProperty()->SetLookupTable(ctf);
            }
        }

        // use default LUT
        else if(index == 0)
        {
            for(int i = 0; i<3; i++)
            {
                this->overlay_modelInfo->orthogonalPlanesModels[i]->GetProperty()->SetLookupTable(this->overlay_modelInfo->defaultLUT);
            }
        }

		this->core()->render();
    }
}

//---------------[ Select base layer LUT for clipping planes ]-----------------\\

void IsosurfaceVisualization::comboBoxBaseLayerLUTChanged()
{
    // select model info matching the dataset
    int index = this->form->comboBoxBaseLayerLUT->currentIndex();

    if(this->current_modelInfo != NULL)
    {
        this->current_modelInfo->baseLayerLUTIndex = index;

		if(index > 0)
        {
            vtkLookupTable* ctf = this->lookUpTables.at(index - 1);
            for(int i = 0; i<3; i++)
            {
                this->current_modelInfo->orthogonalPlanesModels[i]->GetProperty()->SetLookupTable(ctf);
            }
        }

        // use default LUT
        else if(index == 0)
        {
            for(int i = 0; i<3; i++)
            {
                this->current_modelInfo->orthogonalPlanesModels[i]->GetProperty()->SetLookupTable(this->current_modelInfo->defaultLUT);
            }
        }

		this->core()->render();
    }
}

//---------------[ Select curvature LUT for mesh ]-----------------\\

void IsosurfaceVisualization::comboBoxCurvatureLUTChanged()
{
    // select model info matching the dataset
    int index = this->form->comboBoxCurvatureLUT->currentIndex();

    if(this->current_modelInfo->prop != NULL)
    {
		this->current_modelInfo->curvatureLUTIndex = index;

        if(index > 0)
        {
			// if curvature was already computed. (otherwise update is necessary)
			if(current_modelInfo->usingCurvature)
			{
				vtkLookupTable* ctf = this->lookUpTables.at(index - 1);
				current_modelInfo->polyDataMapper->SetLookupTable(ctf);
				current_modelInfo->polyDataMapper->Update();
			}
			else
				bModelDirty = true;
        }

        // disable curvature
        else if(index == 0)
        {
			// set blank LUT
			vtkLookupTable* ctf = vtkLookupTable::New();
            current_modelInfo->polyDataMapper->SetLookupTable(ctf);
        }

		this->core()->render();
    }
}

//---------------[ Select curvature type for mesh ]-----------------\\

void IsosurfaceVisualization::comboBoxCurvatureTypeChanged()
{
    this->current_modelInfo->curvatureType = this->form->comboBoxCurvatureType->currentIndex();
    bModelDirty = true;
}

//------------------------[ Change visibility of model ]-----------------------\\

void IsosurfaceVisualization::checkBoxVisibleChanged(bool checked)
{
    current_modelInfo->visible = checked;
    if(current_modelInfo->prop)
    {
        current_modelInfo->prop->SetVisibility(checked);
        this->core()->render();
    }
}

//------------------------[ Change rendering style of model ]-----------------------\\

void IsosurfaceVisualization::comboBoxStyleChanged()
{
    // select model info matching the dataset
    int index = this->form->comboBoxStyle->currentIndex();

    current_modelInfo->renderStyle = index;
    if(current_modelInfo->prop != NULL)
    {
        switch(index)
        {
        case 0:
            current_modelInfo->prop->GetProperty()->SetRepresentationToSurface();
            break;

        case 1:
            current_modelInfo->prop->GetProperty()->SetRepresentationToWireframe();
            break;

        case 2:
            current_modelInfo->prop->GetProperty()->SetRepresentationToPoints();
            break;
        }
        this->core()->render();
    }
}

//------------------[ Click on the update button to start processing mesh ]-----------------------\\

void IsosurfaceVisualization::buttonUpdateClicked()
{
    if(this->form->inputMinimumThreshold->value() > this->form->inputMaximumThreshold->value())
    {
        this->form->inputMinimumThreshold->setValue(this->form->inputMaximumThreshold->value());
        this->core()->render();
    }
    updateRenderingModels();
}

//------------------[ Change maximum threshold ]-----------------------\\

void IsosurfaceVisualization::inputMaximumThresholdChanged(double value)
{
    current_modelInfo->maximumThreshold = value;
    bModelDirty = true;
    this->form->inputMinimumThreshold->setMaximum(value);
}

//------------------[ Change minimum threshold ]-----------------------\\

void IsosurfaceVisualization::inputMinimumThresholdChanged(double value)
{
    current_modelInfo->minimumThreshold = value;
    bModelDirty = true;
    this->form->inputMaximumThreshold->setMinimum(value);
}

//------------------[ Change smoothing value ]-----------------------\\

void IsosurfaceVisualization::inputSmoothingChanged(double value)
{
    current_modelInfo->smoothing = value;
    bModelDirty = true;
}

//------------------[ Change color ]-----------------------\

void IsosurfaceVisualization::inputColorChanged()
{
    /*int color_hex = value.toInt(0,16);
    current_modelInfo->color[0] = ((color_hex >> 16) & 0xff) / 255.0; // red
    current_modelInfo->color[1] = ((color_hex >> 8) & 0xff) / 255.0; // green
    current_modelInfo->color[2] = (color_hex & 0xff) / 255.0; // blue
    current_modelInfo->colorString = value;*/

    QColor oldColor;
    oldColor.setRgbF(current_modelInfo->color[0], current_modelInfo->color[1], current_modelInfo->color[2]);
    QColor newColor = QColorDialog::getColor(oldColor, 0);

    current_modelInfo->color[0] = newColor.redF();
    current_modelInfo->color[1] = newColor.greenF();
    current_modelInfo->color[2] = newColor.blueF();
}

//------------------[ Change alpha ]-----------------------\

void IsosurfaceVisualization::inputAlphaChanged(double value)
{
    current_modelInfo->alpha = value;
}

//------------------[ Change specularity ]-----------------------\

void IsosurfaceVisualization::inputSpecularChanged(double value)
{
    current_modelInfo->specular = value;
}

//------------------[ Change mesh reduction ]-----------------------\

void IsosurfaceVisualization::inputReductionChanged(double value)
{
    current_modelInfo->reduction = value;
    bModelDirty = true;
}

//---------------[ Change clipping plane slider in x direction ]-------------------\

void IsosurfaceVisualization::horizontalSliderXChanged(int value)
{
    this->form->spinX->setValue(value);
    this->updateClippingPlaneSlider(0,value);

	//this->current_modelInfo->prop->GetProperty()->AddShaderVariable("ClipX",(double)(value) );
	//this->core()->render();
}

//---------------[ Change clipping plane slider in y direction ]-------------------\

void IsosurfaceVisualization::horizontalSliderYChanged(int value)
{
    this->form->spinY->setValue(value);
    this->updateClippingPlaneSlider(1,value);

		//this->current_modelInfo->prop->GetProperty()->AddShaderVariable("ClipY",(double)(value) );
	//this->core()->render();
}

//---------------[ Change clipping plane slider in z direction ]-------------------\

void IsosurfaceVisualization::horizontalSliderZChanged(int value)
{
    this->form->spinZ->setValue(value);
    this->updateClippingPlaneSlider(2,value);

		//this->current_modelInfo->prop->GetProperty()->AddShaderVariable("ClipZ",(double)(value) );
	//this->core()->render();
}

//---------------[ Change clipping plane checkbox in x direction ]-------------------\

void IsosurfaceVisualization::checkBoxXChanged(bool checked)
{
    this->updateClippingPlaneEnabled(0,checked);
}

//---------------[ Change clipping plane checkbox in y direction ]-------------------\

void IsosurfaceVisualization::checkBoxYChanged(bool checked)
{
    this->updateClippingPlaneEnabled(1,checked);
}

//---------------[ Change clipping plane checkbox in z direction ]-------------------\

void IsosurfaceVisualization::checkBoxZChanged(bool checked)
{
    this->updateClippingPlaneEnabled(2,checked);
}

//---------------[ Change clipping plane flip checkbox in x direction ]-------------------\

void IsosurfaceVisualization::checkBoxFlipXChanged(bool checked)
{
    current_modelInfo->flipped[0] = checked;
    this->updateClippingPlaneSlider(0,current_modelInfo->clippingValues[0]);
}

//---------------[ Change clipping plane flip checkbox in y direction ]-------------------\

void IsosurfaceVisualization::checkBoxFlipYChanged(bool checked)
{
    current_modelInfo->flipped[1] = checked;
    this->updateClippingPlaneSlider(1,current_modelInfo->clippingValues[1]);
}

//---------------[ Change clipping plane flip checkbox in z direction ]-------------------\

void IsosurfaceVisualization::checkBoxFlipZChanged(bool checked)
{
    current_modelInfo->flipped[2] = checked;
    this->updateClippingPlaneSlider(2,current_modelInfo->clippingValues[2]);
}

//---------------[ Change clipping plane spin value in X direction ]-------------------\

void IsosurfaceVisualization::spinXChanged(int value)
{
   // this->form->horizontalSliderX->setValue(value);
    //this->updateClippingPlaneSlider(0,value);
}

//---------------[ Change clipping plane spin value in Y direction ]-------------------\

void IsosurfaceVisualization::spinYChanged(int value)
{
    //this->form->horizontalSliderY->setValue(value);
   // this->updateClippingPlaneSlider(1,value);
}

//---------------[ Change clipping plane spin value in Z direction ]-------------------\

void IsosurfaceVisualization::spinZChanged(int value)
{
   // this->form->horizontalSliderZ->setValue(value);
    //this->updateClippingPlaneSlider(2,value);
}

//---------------[ Save 3D mesh ]-------------------\

void IsosurfaceVisualization::buttonSaveMeshClicked()
{
    if(current_modelInfo->polydata != NULL)
    {
        QString fileName = QFileDialog::getSaveFileName(0, "Save 3D mesh","","VTK (*.vtk)");

		if(fileName.isEmpty())
			return;

        QByteArray ba = fileName.toLocal8Bit();
        char* fileName_c = ba.data();

        // export model
        vtkPolyDataWriter* polyDataWriter = vtkPolyDataWriter::New();
        polyDataWriter->SetInput(current_modelInfo->polydata);
        polyDataWriter->SetFileName(fileName_c);
        polyDataWriter->Write();
    }
}

void IsosurfaceVisualization::checkBoxInvertClippingChanged(bool checked)
{
    if(current_modelInfo->extractPolyFunc != NULL)
    {
        current_modelInfo->invertClipping = checked;

        if(checked == false)
            current_modelInfo->extractPolyFunc->ExtractInsideOn();
        else
            current_modelInfo->extractPolyFunc->ExtractInsideOff();

        this->core()->render();
    }
}

void IsosurfaceVisualization::checkBoxLargestComponentChanged(bool checked)
{
    current_modelInfo->selectLargestComponent = checked;
    bModelDirty = true;
}

void IsosurfaceVisualization::checkBoxAlignPlanesToPickChanged(bool checked)
{
    current_modelInfo->alignPlanesToPick = checked;
}

void IsosurfaceVisualization::buttonSetPointAClicked()
{
    setMeasuredPoint(0);
}

void IsosurfaceVisualization::buttonSetPointBClicked()
{
    setMeasuredPoint(1);
}

void IsosurfaceVisualization::setMeasuredPoint(int id)
{
    MeasuredPoint* point = this->measuredPointList.at(id);

    if(!point->set)
    {
        this->assembly->AddPart(point->sphere);
    }

    point->set = true;

    point->x = this->clickedPoint[0];
    point->y = this->clickedPoint[1];
    point->z = this->clickedPoint[2];

    point->sphere->SetPosition(point->x,point->y,point->z);

    if(id == 0)
    {
        this->form->inputXPointA->setValue(point->x);
        this->form->inputYPointA->setValue(point->y);
        this->form->inputZPointA->setValue(point->z);
    }
    else
    {
        this->form->inputXPointB->setValue(point->x);
        this->form->inputYPointB->setValue(point->y);
        this->form->inputZPointB->setValue(point->z);
    }

    calculateDistance();

    this->core()->render();
}

void IsosurfaceVisualization::calculateDistance()
{
    MeasuredPoint* pointA = this->measuredPointList.at(0);
    MeasuredPoint* pointB = this->measuredPointList.at(1);

    if(!pointA->set || !pointB->set)
        return;

    double distance =   sqrt( (pointA->x - pointB->x)*(pointA->x - pointB->x) +
                        (pointA->y - pointB->y)*(pointA->y - pointB->y) +
                        (pointA->z - pointB->z)*(pointA->z - pointB->z) );
    QString labeltext = QString("Measured distance: %1 mm").arg(distance,0,'f',2);
	QString labeltext_short = QString("%1 mm").arg(distance,0,'f',2);
    this->form->measuredDistanceLabel->setText(labeltext);

    // renew line
    vtkSmartPointer<vtkLineSource> lineSource =
        vtkSmartPointer<vtkLineSource>::New();
    lineSource->SetPoint1(pointA->x,pointA->y,pointA->z);
    lineSource->SetPoint2(pointB->x,pointB->y,pointB->z);
    vtkSmartPointer<vtkPolyDataMapper> lineMapper =
        vtkSmartPointer<vtkPolyDataMapper>::New();
    lineMapper->SetInputConnection(lineSource->GetOutputPort());

    vtkActor* lineActor =
        vtkActor::New();
    lineActor->SetMapper(lineMapper);
    lineActor->GetProperty()->SetLineStipplePattern(0xFF00);

    this->assembly->AddPart(lineActor);
    if(measuredLine != NULL)
        this->assembly->RemovePart(measuredLine);
    measuredLine = lineActor;



	measuredLabelPoints = vtkPoints::New();
	measuredLabelStrings = vtkStringArray::New();

	measuredLabelStrings->SetName("labels");
	measuredLabelPoints->InsertNextPoint(pointA->x,pointA->y,pointA->z+5);
	measuredLabelStrings->InsertNextValue(this->form->lineEditNamePointA->text().toLocal8Bit().constData());

	measuredLabelPoints->InsertNextPoint(pointB->x,pointB->y,pointB->z+5);
	measuredLabelStrings->InsertNextValue(this->form->lineEditNamePointB->text().toLocal8Bit().constData());

	measuredLabelPoints->InsertNextPoint(pointA->x + (pointB->x-pointA->x)/2.0,
							pointA->y + (pointB->y-pointA->y)/2.0,
							pointA->z + (pointB->z-pointA->z)/2.0+5);
	measuredLabelStrings->InsertNextValue(labeltext_short.toLocal8Bit().constData());

    if(measuredLabels != NULL)
        this->assembly->RemovePart(measuredLabels);

	measuredLabels = this->GenerateLabels(measuredLabelPoints,measuredLabelStrings);
	this->assembly->AddPart(measuredLabels);


    if (this->currentElectrodesColor.isValid())
    {
        measuredLine->GetProperty()->SetColor(this->currentElectrodesColor.redF(), this->currentElectrodesColor.greenF(), this->currentElectrodesColor.blueF());
        pointA->sphere->GetProperty()->SetColor(this->currentElectrodesColor.redF(), this->currentElectrodesColor.greenF(), this->currentElectrodesColor.blueF());
        pointB->sphere->GetProperty()->SetColor(this->currentElectrodesColor.redF(), this->currentElectrodesColor.greenF(), this->currentElectrodesColor.blueF());
    }


    // remove old blobs
    /*QList<vtkActor*>::iterator i;
    for(i = this->depthElectrodeBlobs.begin(); i!=this->depthElectrodeBlobs.end(); i++)
    {
        this->assembly->RemovePart((*i));
    }
    this->depthElectrodeBlobs.clear();

    // create new blobs
    float stepSize = 4; // 1.5 mm depth electrode seperation
    int steps = ceil(distance/stepSize);

    vec3* line = new vec3;
    line->x = (pointA->x - pointB->x)/distance;
    line->y = (pointA->y - pointB->y)/distance;
    line->z = (pointA->z - pointB->z)/distance;

    for(int j = 0; j<steps; j++)
    {
        vtkSmartPointer<vtkSphereSource> diskSource =
            vtkSmartPointer<vtkSphereSource>::New();
        diskSource->SetRadius(2);

        vtkSmartPointer<vtkPolyDataMapper> diskMapper =
            vtkSmartPointer<vtkPolyDataMapper>::New();
        diskMapper->SetInputConnection(diskSource->GetOutputPort());

        vtkActor* diskActor =
            vtkActor::New();
        diskActor->SetMapper(diskMapper);

        diskActor->SetPosition(-line->x * j * stepSize + pointA->x, -line->y * j * stepSize + pointA->y, -line->z * j * stepSize + pointA->z );

        if (this->currentElectrodesColor.isValid())
        {
            diskActor->GetProperty()->SetColor(this->currentElectrodesColor.redF(), this->currentElectrodesColor.greenF(), this->currentElectrodesColor.blueF());
        }

        this->assembly->AddPart(diskActor);
        this->depthElectrodeBlobs.append(diskActor);
    }*/
}

void IsosurfaceVisualization::buttonSetLineColorClicked()
{
    double oldColorRGB[3];
    QColor oldColor;
    //vtkActor* electrodeActor = depthElectrodeBlobs.at(0);
    measuredLine->GetProperty()->GetColor(oldColorRGB);
    oldColor.setRgbF(oldColorRGB[0], oldColorRGB[1], oldColorRGB[2]);

    // Use a color dialog to get the new color
    this->currentElectrodesColor = QColorDialog::getColor(oldColor, 0);

    measuredLine->GetProperty()->SetColor(this->currentElectrodesColor.redF(), this->currentElectrodesColor.greenF(), this->currentElectrodesColor.blueF());

    MeasuredPoint* pointA = this->measuredPointList.at(0);
    MeasuredPoint* pointB = this->measuredPointList.at(1);
    pointA->sphere->GetProperty()->SetColor(this->currentElectrodesColor.redF(), this->currentElectrodesColor.greenF(), this->currentElectrodesColor.blueF());
    pointB->sphere->GetProperty()->SetColor(this->currentElectrodesColor.redF(), this->currentElectrodesColor.greenF(), this->currentElectrodesColor.blueF());

    // set current electrodes to that color
    /*QList<vtkActor*>::iterator i;
    for(i = this->depthElectrodeBlobs.begin(); i!=this->depthElectrodeBlobs.end(); i++)
    {
        if (this->currentElectrodesColor.isValid())
        {
            (*i)->GetProperty()->SetColor(this->currentElectrodesColor.redF(), this->currentElectrodesColor.greenF(), this->currentElectrodesColor.blueF());
        }
    }*/
}

void IsosurfaceVisualization::buttonSaveMeasurementClicked()
{
    this->depthElectrodeBlobs.clear();
}

void IsosurfaceVisualization::lineEditNamePointAChanged(QString value)
{
	this->measuredLabelStrings->SetValue(0,value.toLocal8Bit().constData());
	this->core()->render();
}

void IsosurfaceVisualization::lineEditNamePointBChanged(QString value)
{
	this->measuredLabelStrings->SetValue(1,value.toLocal8Bit().constData());
	this->core()->render();
}

void IsosurfaceVisualization::comboBoxFiberDataChanged()
{
    // select model info matching the dataset
    int index = this->form->comboBoxFiberData->currentIndex() - 1;

    // Update selected fiber data in settings
    current_modelInfo->selectedFiberData = index;

    // "None" value
    if(index < 0)
    {
        // disable gui
        this->form->sliderFiberChoice->setEnabled(false);
        this->form->spinFiberChoice->setEnabled(false);
        this->form->sliderFiberRefinement->setEnabled(false);
        this->form->spinFiberRefinement->setEnabled(false);
    }
    else
    {
        // enable gui
        this->form->sliderFiberChoice->setEnabled(true);
        this->form->spinFiberChoice->setEnabled(true);
        this->form->sliderFiberRefinement->setEnabled(true);
        this->form->spinFiberRefinement->setEnabled(true);

        // Get data struct
        SortedFibers* sortedFibers = this->sortedFibersList.at(index);

        // process data set
        processFiberAnteriorSorting(sortedFibers);

        // Set GUI value limits
        int amountOfFibers = sortedFibers->selectedLines.length() - 1; // 100 or less
        this->form->sliderFiberChoice->setMaximum(amountOfFibers);
        this->form->spinFiberChoice->setMaximum(amountOfFibers);

        // update selected point b
        fiberSelectUpdate(sortedFibers->userSelectedLine);







        secondWindow.setGeometry(0, 0, 1150, 600);


        // QVTK set up and initialization
        QVTKWidget *qvtkWidget = new QVTKWidget(&secondWindow);
        QVTKWidget *qvtkWidget2 = new QVTKWidget(&secondWindow);

        // Set up my 2D world...
        VTK_CREATE(vtkContextView, view); // This contains a chart object
        view->SetInteractor(qvtkWidget->GetInteractor());
        qvtkWidget->SetRenderWindow(view->GetRenderWindow());

        VTK_CREATE(vtkContextView, view2); // This contains a chart object
        view2->SetInteractor(qvtkWidget2->GetInteractor());
        qvtkWidget2->SetRenderWindow(view2->GetRenderWindow());

        // Create a table with some points in it...
        VTK_CREATE(vtkTable, table);
        VTK_CREATE(vtkFloatArray, arrX);
        arrX->SetName("X Axis");
        table->AddColumn(arrX);
        VTK_CREATE(vtkFloatArray, arrC);
        arrC->SetName("Average score");
        table->AddColumn(arrC);
        VTK_CREATE(vtkFloatArray, arrS);
        arrS->SetName("Local average score");
        table->AddColumn(arrS);


        MeasuredPoint* pointA = this->measuredPointList.at(0);

        int length = sortedFibers->selectedLines.length();
        table->SetNumberOfRows(length);
        for(int j = 0; j<length; j++)
        {
            int fiberLength = sortedFibers->selectedLines.at(j)->scalarData.length();

            int anteriorIndex = sortedFibers->selectedLines.at(j)->anteriorPointIndex;
            Vec3* pointB = sortedFibers->selectedLines.at(j)->data.at(anteriorIndex);

            double distance =   sqrt( (pointA->x - pointB->x)*(pointA->x - pointB->x) +
                                (pointA->y - pointB->y)*(pointA->y - pointB->y) +
                                (pointA->z - pointB->z)*(pointA->z - pointB->z) );

            // distance axis
            table->SetValue(j, 0, distance);

            // average score
            double avgScore = 0.0;
            for(int i = 0; i<fiberLength; i++)
            {
                avgScore += sortedFibers->selectedLines.at(j)->scalarData.at(i);
            }
            avgScore /= fiberLength;
            table->SetValue(j, 1, avgScore);

            // local average score
            double localAvgScore = 0.0;
            int pointRange = 25;
            int minBound = std::max(0,anteriorIndex - pointRange);
            int maxBound = std::min(fiberLength,anteriorIndex + pointRange);
            for(int i = minBound; i<maxBound; i++)
            {
                localAvgScore += sortedFibers->selectedLines.at(j)->scalarData.at(i);
            }
            localAvgScore /= maxBound - minBound;
            table->SetValue(j, 2, localAvgScore);


            //std::cout << " j: " << j << " distance: " << distance
            //<< " avg score: " << avgScore
            //<< " local avg score: " << localAvgScore << std::endl;

        }
        table->Update();

        // Add multiple line plots, setting the colors etc
        vtkSmartPointer<vtkChartXY> chart = vtkSmartPointer<vtkChartXY>::New();
        view->GetScene()->AddItem(chart);
        vtkPlot *line = chart->AddPlot(vtkChart::POINTS);
        line->SetInput(table, 0, 1);
        line->SetColor(255, 0, 0, 255);
        line->SetWidth(2.0);
        chart->GetAxis(vtkAxis::LEFT)->SetTitle("Connectivity measure (-)");
        chart->GetAxis(vtkAxis::BOTTOM)->SetTitle("Distance (mm)");
        chart->SetTitle("Average score");
        //chart->SetActionToButton(vtkChart::SELECT,vtkContextMouseEvent::LEFT_BUTTON);
        //chart->SetActionToButton(vtkChart::ZOOM,vtkContextMouseEvent::MIDDLE_BUTTON);
        //chart->SetActionToButton(vtkChart::PAN,vtkContextMouseEvent::RIGHT_BUTTON);
        chart->SetClickActionToButton(vtkChart::SELECT, vtkContextMouseEvent::LEFT_BUTTON);
        chart->SetClickActionToButton(vtkChart::NOTIFY, vtkContextMouseEvent::RIGHT_BUTTON);

        vtkSmartPointer<vtkChartXY> chart2 = vtkSmartPointer<vtkChartXY>::New();
        view2->GetScene()->AddItem(chart2);
        line = chart2->AddPlot(vtkChart::POINTS);
        line->SetInput(table, 0, 2);
        line->SetColor(0, 255, 0, 255);
        line->SetWidth(2.0);
        chart2->GetAxis(vtkAxis::LEFT)->SetTitle("Connectivity measure (-)");
        chart2->GetAxis(vtkAxis::BOTTOM)->SetTitle("Distance (mm)");
        chart2->SetTitle("Average local score");
        //chart2->SetClickActionToButton(vtkChart::SELECT, vtkContextMouseEvent::LEFT_BUTTON);
        //chart2->SetClickActionToButton(vtkChart::NOTIFY, vtkContextMouseEvent::RIGHT_BUTTON);
        chart2->SetSelectionMode(vtkContextScene::SELECTION_TOGGLE);

        // Now lets try to add a table view
        QVBoxLayout *layout = new QVBoxLayout(&secondWindow);
        layout->addWidget(qvtkWidget);
        layout->addWidget(qvtkWidget2);

        secondWindow.raise();
        secondWindow.show();





    }
}

void IsosurfaceVisualization::processFiberAnteriorSorting(SortedFibers* sortedFibers)
{
    // Look if dataset is already processed ...
    if(sortedFibers->selectedLines.length() != 0)
        return;

    // Get the polydata from the data set
    vtkPolyData * polydata = sortedFibers->ds->getVtkPolyData();

	// Get the transformation matrix
    vtkObject* tfm;
    vtkMatrix4x4* transformationMatrix;
    if (sortedFibers->ds->getAttributes()->getAttribute("transformation matrix", tfm ))
    {
        transformationMatrix = vtkMatrix4x4::SafeDownCast(tfm);
        if (transformationMatrix == 0)
        {
            this->core()->out()->logMessage("not a valid transformation matrix");
            return;
        }
    }

    // Get fiber tracts
    vtkCellArray * fibers = polydata->GetLines();
    vtkIdType numberOfFibers = fibers->GetNumberOfCells();

    vtkPointData* pointdata = polydata->GetPointData();
    pointdata->Print(std::cout);

    // Get scalar values (if present)
    // Check if the fibers contain point data
    vtkDataArray * scalars;
    bool hasScalars = false;
    if (polydata->GetPointData())
    {
        if (scalars = polydata->GetPointData()->GetScalars())
        {
            // Scalar array should have as many points as the input fiber set, and at least one component
            if (scalars->GetNumberOfTuples() == polydata->GetNumberOfPoints() && scalars->GetNumberOfComponents() > 0)
            {
                hasScalars = true;
            }
        }
    }

    // Map used to store fiber indices (value) and their anterior point (key)
	QMap<double, FiberData*> fiberMap;

    vtkIdType numberOfPoints;
	vtkIdType * pointList;
	double currentPoint[3];
    fibers->InitTraversal();

    // Loop through all fibers
	for (vtkIdType fiberId = 0; fiberId < numberOfFibers; ++fiberId)
	{
		// Get the number of points and the list of point IDs of the current fiber
		fibers->GetNextCell(numberOfPoints, pointList);

		// Do nothing if the fiber is empty
		if (numberOfPoints == 0)
			continue;

        double mostAnteriorPoint = -1e32; // low value
        int anteriorPointIndex = -1;

        // Create fiber data struct
        FiberData* fiberData = new FiberData;

		// Loop through all points of the current fiber
		for (vtkIdType pointId = 0; pointId < numberOfPoints; ++pointId)
		{
            // Get the coordinates of the current point
            polydata->GetPoint(pointList[pointId], currentPoint);

			// Update most anterior point
            if(currentPoint[1] > mostAnteriorPoint)
            {
                mostAnteriorPoint = currentPoint[1];
                anteriorPointIndex = pointId;
            }

            // Do the point transformation
            double vec[4] = {currentPoint[0],currentPoint[1],currentPoint[2],1};
			transformationMatrix->MultiplyPoint(vec,vec);

            // Save fiber point as Vec3 for easier management
            Vec3* vec3 = new Vec3;
            vec3->x = vec[0];
            vec3->y = vec[1];
            vec3->z = vec[2];

            // Save to fiber data struct
            fiberData->data.append(vec3);
		}

		// If it contains scalar data, copy it
		if(hasScalars)
		{
		    for (vtkIdType pointId = 0; pointId < numberOfPoints; ++pointId)
            {
                fiberData->scalarData.append(scalars->GetTuple1(pointList[pointId]));
            }
		}

        // Set anterior point index in struct
        fiberData->anteriorPointIndex = anteriorPointIndex;

        // Set refinement point at zero
        fiberData->userPointRefinement = 0;

        // Add fiber data in QMap for sorting
		fiberMap.insert(mostAnteriorPoint, fiberData);
	}

    // Select most anterior fibers
    int rFiberIndex = 0;
    int numberOfOutputFibers = 100;
    QMap<double, FiberData*>::iterator rIter = fiberMap.end();
    while (rIter != fiberMap.begin())
	{
		// Decrement the iterator
		rIter--;

        // Add the fiber data struct to the sorted fibers list
        sortedFibers->selectedLines.append(rIter.value());

		// Break if we've reached the desired number of fibers
		if (++rFiberIndex == numberOfOutputFibers)
			break;
	}
}

void IsosurfaceVisualization::fiberSelectUpdate(int value)
{
    // Get data struct
    SortedFibers* sortedFibers = this->sortedFibersList.at(current_modelInfo->selectedFiberData);

    // get selected fiber index
    sortedFibers->userSelectedLine = std::max(value,0);

    // fiber data
    FiberData* fiberData = sortedFibers->selectedLines.at(sortedFibers->userSelectedLine);

    // set gui values
    this->form->sliderFiberChoice->setValue(value);
    this->form->spinFiberChoice->setValue(value);
    this->form->sliderFiberRefinement->setValue(fiberData->userPointRefinement);
    this->form->spinFiberRefinement->setValue(fiberData->userPointRefinement);

    // Set GUI limits on fiber refinement
    int refinementValue_min = -fiberData->anteriorPointIndex;
    int refinementValue_max = fiberData->data.length() - fiberData->anteriorPointIndex - 1;
    this->form->sliderFiberRefinement->setMinimum(refinementValue_min);
    this->form->spinFiberRefinement->setMinimum(refinementValue_min);
    this->form->sliderFiberRefinement->setMaximum(refinementValue_max);
    this->form->spinFiberRefinement->setMaximum(refinementValue_max);

    // update fiber point
    fiberPointSelect();
}

void IsosurfaceVisualization::fiberRefinementUpdate(double value)
{
    // Get data struct
    SortedFibers* sortedFibers = this->sortedFibersList.at(current_modelInfo->selectedFiberData);

    // fiber data
    FiberData* fiberData = sortedFibers->selectedLines.at(sortedFibers->userSelectedLine);

    // Update refinement value
    fiberData->userPointRefinement = value;

    // set GUI values
    this->form->sliderFiberRefinement->setValue(fiberData->userPointRefinement);
    this->form->spinFiberRefinement->setValue(fiberData->userPointRefinement);

    // update fiber point
    fiberPointSelect();
}

void IsosurfaceVisualization::fiberRefinementUpdate(int value)
{
    fiberRefinementUpdate((double)value);
}

void IsosurfaceVisualization::fiberPointSelect()
{
    // Get data struct
    SortedFibers* sortedFibers = this->sortedFibersList.at(current_modelInfo->selectedFiberData);

    // fiber data
    FiberData* fiberData = sortedFibers->selectedLines.at(sortedFibers->userSelectedLine);

    Vec3* pointA = fiberData->data.at(floor(fiberData->anteriorPointIndex + fiberData->userPointRefinement));
    Vec3* pointB = fiberData->data.at(ceil(fiberData->anteriorPointIndex + fiberData->userPointRefinement));

    // set point B
    // do a linear interpolation for the refinement setting
    double loc_pos = fmod(fiberData->anteriorPointIndex + fiberData->userPointRefinement,1);
    this->clickedPoint[0] = pointB->x * loc_pos + pointA->x * (1 - loc_pos);
    this->clickedPoint[1] = pointB->y * loc_pos + pointA->y * (1 - loc_pos);
    this->clickedPoint[2] = pointB->z * loc_pos + pointA->z * (1 - loc_pos);

    // update point b
    setMeasuredPoint(1);
}

///
///     vISTe communication
///

//-----------[ Returns visualization component as VTK object ]---------------\\
//
vtkProp * IsosurfaceVisualization::getVtkProp()
{
    return this->assembly;
}

//-----------------[ Returns GUI component as Qt widget ]---------------\\
//
QWidget * IsosurfaceVisualization::getGUI()
{
    return this->widget;
}

}

Q_EXPORT_PLUGIN2(libIsosurfaceVisualization, bmia::IsosurfaceVisualization)
