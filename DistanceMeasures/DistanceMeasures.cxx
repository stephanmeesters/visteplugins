#include "DistanceMeasures.h"

#define VTK_CREATE(type, name) \
  vtkSmartPointer<type> name = vtkSmartPointer<type>::New()

namespace bmia
{

///
///      INITIALIZATION
///

//------------------------[ Plugin constructor ]-----------------------\\

DistanceMeasures::DistanceMeasures() : plugin::AdvancedPlugin("DistanceMeasures")
{
    this->widget = NULL;
    this->form   = NULL;
}

//------------------------[ Plugin destructor ]-----------------------\\

DistanceMeasures::~DistanceMeasures()
{
    delete this->widget;
    delete this->form;
}

//------------------------[ Initialization ]-----------------------\\

void DistanceMeasures::init()
{
    this->widget = new QWidget();
    this->form = new Ui::DistanceMeasuresForm();
    this->form->setupUi(this->widget);

    // Link events in the GUI to function calls
    this->connectAll();
    this->assembly = vtkPropAssembly::New();

    this->SetupPointers();
}

void DistanceMeasures::SetupPointers()
{
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
}

///
///      DATA I/O
///

//------------------------[ Dataset added ]-----------------------\\

void DistanceMeasures::dataSetAdded(data::DataSet * d)
{
    // Assert the data set pointer (should never be NULL)
    Q_ASSERT(d);

	// Get the kind of the data set
    QString kind = d->getKind();

    // Load fiber dataset
    if (kind == "fibers")
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

	// Load settings (should be only one)
	else if(kind == "settings")
	{
		printf("SETTIGNS ADDED!");
		this->settings = d;
	}
}

//------------------------[ Dataset changed ]-----------------------\\

void DistanceMeasures::dataSetChanged(data::DataSet * d)
{
    // Assert the data set pointer (should never be NULL)
    Q_ASSERT(d);

	// Get the kind of the data set
    QString kind = d->getKind();

	// Update fibers
	if (kind == "fibers")
	{
		// to-do
	}

	// Update settings
	else if(kind == "settings")
	{
		double pos[3];
		d->getAttributes()->getAttribute("SlicePosX",pos[0]);
		d->getAttributes()->getAttribute("SlicePosY",pos[1]);
		d->getAttributes()->getAttribute("SlicePosZ",pos[2]);

		this->settings = d;

		printf("new settings: x:%f y:%f z:%f \n",pos[0],pos[1],pos[2]);
	}
}

//------------------------[ Dataset removed ]-----------------------\\

void DistanceMeasures::dataSetRemoved(data::DataSet * d)
{
    // Assert the data set pointer (should never be NULL)
    Q_ASSERT(d);

	// Get the kind of the data set
    QString kind = d->getKind();

    // Remove fiber dataset
    if (kind == "fibers")
	{
	    // Check if the data set exists
		int dsIndex = this->FindInputDataSet(d);

        // Does not exist, return
		if (dsIndex == -1)
			return;

        // Remove from UI combobox for selection of overlay
        this->form->comboBoxFiberData->removeItem(dsIndex);

        // Remove from collection
        this->sortedFibersList.removeAt(dsIndex);
	}
}

int DistanceMeasures::FindInputDataSet(data::DataSet * ds)
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

//------------------------[ Create text labels ]-----------------------\\

vtkActor2D* DistanceMeasures::GenerateLabels(vtkSmartPointer<vtkPoints> points, vtkSmartPointer<vtkStringArray> labels)
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
///     GUI CALLBACKS
///

void DistanceMeasures::setMeasuredPoint(int id)
{
    // Get measured point struct
    MeasuredPoint* point = this->measuredPointList.at(id);
    if(!point->set)
    {
        this->assembly->AddPart(point->sphere);
        point->set = true;
    }

    // Get point coordinates from dataset
    double pos[3];
    this->settings->getAttributes()->getAttribute("SlicePosX",pos[0]);
    this->settings->getAttributes()->getAttribute("SlicePosY",pos[1]);
    this->settings->getAttributes()->getAttribute("SlicePosZ",pos[2]);

    // Save new point
    point->x = pos[0];
    point->y = pos[1];
    point->z = pos[2];

    // Update sphere location
    point->sphere->SetPosition(point->x,point->y,point->z);

    // Update forms
    if(id == 0)
    {
        this->form->inputXPointA->setValue(pos[0]);
        this->form->inputYPointA->setValue(pos[1]);
        this->form->inputZPointA->setValue(pos[2]);
    }
    else
    {
        this->form->inputXPointB->setValue(pos[0]);
        this->form->inputYPointB->setValue(pos[1]);
        this->form->inputZPointB->setValue(pos[2]);
    }

    // Calculate distance and update graphics
    calculateDistance();

    // Render scene
    this->core()->render();
}

void DistanceMeasures::calculateDistance()
{
    // Get points
    MeasuredPoint* pointA = this->measuredPointList.at(0);
    MeasuredPoint* pointB = this->measuredPointList.at(1);

    // If one of the points isnt set, abort
    if(!pointA->set || !pointB->set)
        return;

    // Compute Euclidean distance
    double distance =   sqrt( (pointA->x - pointB->x)*(pointA->x - pointB->x) +
                        (pointA->y - pointB->y)*(pointA->y - pointB->y) +
                        (pointA->z - pointB->z)*(pointA->z - pointB->z) );

    // Set distance text string
    QString labeltext = QString("Measured distance: %1 mm").arg(distance,0,'f',2);
	QString labeltext_short = QString("%1 mm").arg(distance,0,'f',2);
    this->form->measuredDistanceLabel->setText(labeltext);

    // Create line
    vtkSmartPointer<vtkLineSource> lineSource =
        vtkSmartPointer<vtkLineSource>::New();
    vtkSmartPointer<vtkPolyDataMapper> lineMapper =
        vtkSmartPointer<vtkPolyDataMapper>::New();
    vtkActor* lineActor =
        vtkActor::New();
    lineSource->SetPoint1(pointA->x,pointA->y,pointA->z);
    lineSource->SetPoint2(pointB->x,pointB->y,pointB->z);
    lineMapper->SetInputConnection(lineSource->GetOutputPort());
    lineActor->SetMapper(lineMapper);
    lineActor->GetProperty()->SetLineStipplePattern(0xFF00);

    // Add line to scene
    this->assembly->AddPart(lineActor);
    if(measuredLine != NULL)
        this->assembly->RemovePart(measuredLine);
    measuredLine = lineActor;

    // Create labels
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

    // Add labels to scene
    if(measuredLabels != NULL)
        this->assembly->RemovePart(measuredLabels);
	measuredLabels = this->GenerateLabels(measuredLabelPoints,measuredLabelStrings);
	this->assembly->AddPart(measuredLabels);
}

void DistanceMeasures::comboBoxFiberDataChanged()
{
    // select model info matching the dataset
    int index = this->form->comboBoxFiberData->currentIndex() - 1;

    // Selected fiber data
    this->selectedFiberData = index;

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
        int amountOfFibers = sortedFibers->selectedLines.length() - 1;
        this->form->sliderFiberChoice->setMaximum(amountOfFibers);
        this->form->spinFiberChoice->setMaximum(amountOfFibers);

        // update selected point b
        fiberSelectUpdate(sortedFibers->userSelectedLine);

        // Clear scalar type list
        this->form->comboBoxScalar->blockSignals(true);
        for(int i = this->form->comboBoxScalar->count()-1; i>=0; i--)
        {
            this->form->comboBoxScalar->removeItem(i);
        }

        // Add new scalar types to list
        vtkPolyData * polydata = sortedFibers->ds->getVtkPolyData();

        // Get number of scalar types
        int numberOfScalarTypes = polydata->GetPointData()->GetNumberOfArrays();
        bool hasScalars = numberOfScalarTypes > 0;

        // Fill scalar list with names
        for(int i = 0; i < numberOfScalarTypes; i++)
        {
            this->form->comboBoxScalar->addItem(polydata->GetPointData()->GetArray(i)->GetName());
        }
        this->form->comboBoxScalar->blockSignals(false);
    }
}

void DistanceMeasures::processFiberAnteriorSorting(SortedFibers* sortedFibers)
{
    // Look if dataset is already processed ...
    //if(sortedFibers->selectedLines.length() != 0)
    //    return;
    // clear previous lines
    sortedFibers->selectedLines.empty();

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

//--------------------------------------------------------------------------

// print results to text file (temporary)
using namespace std;
ofstream dataFile("/home/linux/Stephan/fiberScoring_output.csv");

//--------------------------------------------------------------------------

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

//--------------------------------------------------------------------------

            float avgVal = 0.0;
            int samplesize = 10;
            for (vtkIdType pointId = std::max(0,anteriorPointIndex-samplesize); pointId < std::min((int)numberOfPoints,anteriorPointIndex+samplesize); ++pointId)
            {
                avgVal += scalars->GetTuple1(pointList[pointId]);
            }
            avgVal /= std::min((int)numberOfPoints,anteriorPointIndex+samplesize) - std::max(0,anteriorPointIndex-samplesize);

            dataFile << avgVal << "\n";

//--------------------------------------------------------------------------

		}

        // Set anterior point index in struct
        fiberData->anteriorPointIndex = anteriorPointIndex;

        // Set refinement point at zero
        fiberData->userPointRefinement = 0;

        // Add fiber data in QMap for sorting
		fiberMap.insert(mostAnteriorPoint, fiberData);
	}

//--------------------------------------------------------------------------

dataFile.close();

//--------------------------------------------------------------------------

    // Select most anterior fibers
    int rFiberIndex = 0;
    int numberOfOutputFibers = numberOfFibers;
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

void DistanceMeasures::fiberSelectUpdate(int value)
{
    // Get data struct
    SortedFibers* sortedFibers = this->sortedFibersList.at(this->selectedFiberData);

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

void DistanceMeasures::fiberRefinementUpdate(double value)
{
    // Get data struct
    SortedFibers* sortedFibers = this->sortedFibersList.at(this->selectedFiberData);

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

void DistanceMeasures::fiberRefinementUpdate(int value)
{
    fiberRefinementUpdate((double)value);
}

void DistanceMeasures::fiberPointSelect()
{
    // Get data struct
    SortedFibers* sortedFibers = this->sortedFibersList.at(this->selectedFiberData);

    // fiber data
    FiberData* fiberData = sortedFibers->selectedLines.at(sortedFibers->userSelectedLine);

    Vec3* pointA = fiberData->data.at(floor(fiberData->anteriorPointIndex + fiberData->userPointRefinement));
    Vec3* pointB = fiberData->data.at(ceil(fiberData->anteriorPointIndex + fiberData->userPointRefinement));

    // set point B
    // do a linear interpolation for the refinement setting
    double loc_pos = fmod(fiberData->anteriorPointIndex + fiberData->userPointRefinement,1);
    this->settings->getAttributes()->addAttribute("SlicePosX", pointB->x * loc_pos + pointA->x * (1 - loc_pos));
    this->settings->getAttributes()->addAttribute("SlicePosY", pointB->y * loc_pos + pointA->y * (1 - loc_pos));
    this->settings->getAttributes()->addAttribute("SlicePosZ", pointB->z * loc_pos + pointA->z * (1 - loc_pos));
    this->core()->data()->dataSetChanged(this->settings);

    // update point b
    setMeasuredPoint(1);
}

void DistanceMeasures::buttonPlotConnectivityClicked()
{
    // Get data struct
    SortedFibers* sortedFibers = this->sortedFibersList.at(this->selectedFiberData);

    // Get points
    MeasuredPoint* pointA = this->measuredPointList.at(0);
    MeasuredPoint* pointB = this->measuredPointList.at(1);

    // If one of the points isnt set, abort
    if(!pointA->set || !pointB->set)
        return;

    QWidget *secondWindow = new QWidget();
    secondWindow->setGeometry(0, 0, 1150, 600);

    // QVTK set up and initialization
    QVTKWidget *qvtkWidget = new QVTKWidget(secondWindow);
    QVTKWidget *qvtkWidget2 = new QVTKWidget(secondWindow);

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


//--------------------------------------------------------------------------

    // print results to text file (temporary)
    using namespace std;
    ofstream dataFile("/home/linux/Stephan/connectivityData.csv");
    //dataFile.open, ios::out, ios::app);
    dataFile << "Distance, AvgScore, LocalScore\n";
    for(int j = 0; j<length; j++)
    {
        dataFile << table->GetValue(j,0) << "," << table->GetValue(j,1) << "," << table->GetValue(j,2) << "\n";
    }
    dataFile.close();

//--------------------------------------------------------------------------

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
    //chart->GetAxis(vtkAxis::LEFT)->SetRange(-0.1, 0.7);
    //chart->GetAxis(vtkAxis::BOTTOM)->SetRange(27,56);
    chart->GetAxis(vtkAxis::LEFT)->SetRange(-0.4, 0.9);
    chart->GetAxis(vtkAxis::BOTTOM)->SetRange(30,72);
    chart->GetAxis(vtkAxis::BOTTOM)->SetBehavior(vtkAxis::FIXED);
    chart->GetAxis(vtkAxis::LEFT)->SetBehavior(vtkAxis::FIXED);

    vtkSmartPointer<vtkChartXY> chart2 = vtkSmartPointer<vtkChartXY>::New();
    view2->GetScene()->AddItem(chart2);
    line = chart2->AddPlot(vtkChart::POINTS);
    line->SetInput(table, 0, 2);
    line->SetColor(0, 255, 0, 255);
    line->SetWidth(2.0);
    chart2->GetAxis(vtkAxis::LEFT)->SetTitle("Connectivity measure (-)");
    chart2->GetAxis(vtkAxis::BOTTOM)->SetTitle("Distance (mm)");
    chart2->SetTitle("Local average score");
    chart2->SetClickActionToButton(vtkChart::SELECT, vtkContextMouseEvent::LEFT_BUTTON);
    chart2->SetClickActionToButton(vtkChart::NOTIFY, vtkContextMouseEvent::RIGHT_BUTTON);
    //chart2->SetSelectionMode(vtkContextScene::SELECTION_TOGGLE);
    //chart2->GetAxis(vtkAxis::LEFT)->SetRange(-2.8, 0.8);
    //chart2->GetAxis(vtkAxis::BOTTOM)->SetRange(27,56);
    chart2->GetAxis(vtkAxis::LEFT)->SetRange(-2.8, 1);
    chart2->GetAxis(vtkAxis::BOTTOM)->SetRange(30,72);
    chart2->GetAxis(vtkAxis::BOTTOM)->SetBehavior(vtkAxis::FIXED);
    chart2->GetAxis(vtkAxis::LEFT)->SetBehavior(vtkAxis::FIXED);

    // Now lets try to add a table view
    QVBoxLayout *layout = new QVBoxLayout(secondWindow);
    layout->addWidget(qvtkWidget);
    layout->addWidget(qvtkWidget2);

    secondWindow->raise();
    secondWindow->show();
}

///
///      GUI CONTROLS
///

//------------------------[ Connect Qt elements ]-----------------------\\

void DistanceMeasures::connectAll()
{
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

    connect(this->form->buttonPlotConnectivity,SIGNAL(clicked()),this,SLOT(buttonPlotConnectivityClicked()));

    connect(this->form->inputXPointA,SIGNAL(clicked()),this,SLOT(updatePointFromSpinBox()));
    connect(this->form->inputYPointA,SIGNAL(clicked()),this,SLOT(updatePointFromSpinBox()));
    connect(this->form->inputZPointA,SIGNAL(clicked()),this,SLOT(updatePointFromSpinBox()));
}

//------------------------[ Disconnect Qt elements ]-----------------------\\

void DistanceMeasures::disconnectAll()
{

}

void DistanceMeasures::updatePointFromSpinBox()
{
    std::cout << 123123123;
    printf("sdsdsdsd");
    this->settings->getAttributes()->addAttribute("SlicePosX", this->form->inputXPointA->value());
    this->settings->getAttributes()->addAttribute("SlicePosY", this->form->inputYPointA->value());
    this->settings->getAttributes()->addAttribute("SlicePosZ", this->form->inputZPointA->value());
    setMeasuredPoint(0);
}

//------------------------[ Set point A ]-----------------------\\

void DistanceMeasures::buttonSetPointAClicked()
{
    setMeasuredPoint(0);
}

//------------------------[ Set point B ]-----------------------\\

void DistanceMeasures::buttonSetPointBClicked()
{
    setMeasuredPoint(1);
}

void DistanceMeasures::lineEditNamePointAChanged(QString value)
{
	this->measuredLabelStrings->SetValue(0,value.toLocal8Bit().constData());
	this->core()->render();
}

void DistanceMeasures::lineEditNamePointBChanged(QString value)
{
	this->measuredLabelStrings->SetValue(1,value.toLocal8Bit().constData());
	this->core()->render();
}

void DistanceMeasures::buttonSetLineColorClicked()
{
    double oldColorRGB[3];
    QColor oldColor;
    measuredLine->GetProperty()->GetColor(oldColorRGB);
    oldColor.setRgbF(oldColorRGB[0], oldColorRGB[1], oldColorRGB[2]);

    // Use a color dialog to get the new color
    QColor currentColor = QColorDialog::getColor(oldColor, 0);

    measuredLine->GetProperty()->SetColor(currentColor.redF(), currentColor.greenF(), currentColor.blueF());

    MeasuredPoint* pointA = this->measuredPointList.at(0);
    MeasuredPoint* pointB = this->measuredPointList.at(1);
    if(!pointA->set || !pointB->set)
        return;
    pointA->sphere->GetProperty()->SetColor(currentColor.redF(), currentColor.greenF(), currentColor.blueF());
    pointB->sphere->GetProperty()->SetColor(currentColor.redF(), currentColor.greenF(), currentColor.blueF());
}

///
///     vISTe communication
///

//-----------[ Returns visualization component as VTK object ]---------------\\
//
vtkProp * DistanceMeasures::getVtkProp()
{
    return this->assembly;
}

//-----------------[ Returns GUI component as Qt widget ]---------------\\
//
QWidget * DistanceMeasures::getGUI()
{
    return this->widget;
}

}

Q_EXPORT_PLUGIN2(libDistanceMeasures, bmia::DistanceMeasures)
