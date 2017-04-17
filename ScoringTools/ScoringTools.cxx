#include "ScoringTools.h"

#define VTK_CREATE(type, name) \
  vtkSmartPointer<type> name = vtkSmartPointer<type>::New()

#define SLIDER_SUBSTEPS 100

namespace bmia
{

///
///      INITIALIZATION
///

//------------------------[ Plugin constructor ]-----------------------\\

ScoringTools::ScoringTools() : plugin::AdvancedPlugin("ScoringTools")
{
    this->widget = NULL;
    this->form   = NULL;
}

//------------------------[ Plugin destructor ]-----------------------\\

ScoringTools::~ScoringTools()
{
    delete this->widget;
    delete this->form;
}

//------------------------[ Initialization ]-----------------------\\

void ScoringTools::init()
{
    this->widget = new QWidget();
    this->form = new Ui::ScoringToolsForm();
    this->form->setupUi(this->widget);

    // Link events in the GUI to function calls
    this->connectAll();
    this->assembly = vtkPropAssembly::New();

    // default selected fiber (none)
    selectedFiberDataset = -1;

    // disable GUI by default
    DisableGUI();
}

///
///      DATA I/O
///

//------------------------[ Dataset added ]-----------------------\\

void ScoringTools::dataSetAdded(data::DataSet * d)
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

        if (d->getAttributes()->hasIntAttribute("isScoreThresholded"))
			return;

        // Create new fiber struct
        SortedFibers* sortedFibers = new SortedFibers;

        // Initialize struct
        sortedFibers->ds = d;
        sortedFibers->ds_processed = NULL;
		sortedFibers->userSelectedLine = 0;
		sortedFibers->selectedScalarType = 0;
		sortedFibers->outputFiberDataName = d->getName().append("_thresholded");
		sortedFibers->processed = false;
		sortedFibers->hasScalars = false;
		sortedFibers->prunePercentage = 100;

        // Add the new data set to the list of currently available fiber sets
        this->sortedFibersList.append(sortedFibers);

        // Add to UI combobox for distance measurements to fibers
        this->form->fibersCombo->addItem(d->getName());
	}

	// Load scalar volume
    else if(kind == "scalar volume")
    {
        // Keep track of the datasets used by this plugin
        this->roiDataSets.append(d);

        // Add to UI combobox for selection of data
        this->form->selectROICombo->addItem(d->getName());
    }
}

//------------------------[ Dataset changed ]-----------------------\\

void ScoringTools::dataSetChanged(data::DataSet * d)
{
    // Assert the data set pointer (should never be NULL)
    Q_ASSERT(d);

	// Get the kind of the data set
    QString kind = d->getKind();

    // Update fiber dataset
//    if (kind == "fibers")
//	{
//	    // Check if the data set exists
//		int dsIndex = this->FindInputDataSet(d);
//
//        // Does not exist, return
//		if (dsIndex == -1)
//			return;
//
//        if (d->getAttributes()->hasIntAttribute("isScoreThresholded"))
//			return;
//
//        this->core()->out()->logMessage(d->getName());
//
//        // Get sorted fibers
//        SortedFibers* sortedFibers = this->sortedFibersList.at(dsIndex );
//
//        // Flag for reprocessing
//        sortedFibers->processed = false;
//
//        // Reload if this is the active fiber dataset
//        if(dsIndex == selectedFiberDataset)
//            SelectFiberDataSet(selectedFiberDataset+1);
//	}
}

//------------------------[ Dataset removed ]-----------------------\\

void ScoringTools::dataSetRemoved(data::DataSet * d)
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

        // Select 'none' in combobox
        this->form->fibersCombo->setCurrentIndex(0);
        this->SelectFiberDataSet(0);

        // Remove from UI combobox for selection of overlay
        this->form->fibersCombo->removeItem(dsIndex+1);

        // Clean up struct
        SortedFibers* sortedFibers = this->sortedFibersList.at(dsIndex);
        sortedFibers->ds = NULL;
        sortedFibers->ds_processed = NULL;
        sortedFibers->selectedLines.clear();
        sortedFibers->scalarThresholdSettings.clear();

        // Remove from collection
        this->sortedFibersList.removeAt(dsIndex);
	}

	// Load scalar volume
    else if(kind == "scalar volume")
    {
        // Check if the data set has been added to this plugin
        if (!(this->roiDataSets.contains(d)))
            return;

        // Get index
        int dsIndex = this->roiDataSets.indexOf(d);

        // Remove from UI combobox for selection of data
        this->form->selectROICombo->removeItem(dsIndex+1);

        // Remove from datasets list
        this->roiDataSets.removeAt(dsIndex);
    }
}

int ScoringTools::FindInputDataSet(data::DataSet * ds)
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

void ScoringTools::SelectFiberDataSet(int index)
{
    // set selected fiber index
    this->selectedFiberDataset = index - 1;

    // no fiber selected
    if(this->selectedFiberDataset == -1)
    {
        DisableGUI();
        return;
    }

    // Clear scalar type list
    this->form->scalarTypeCombo->blockSignals(true);
    for(int i = this->form->scalarTypeCombo->count()-1; i>=0; i--)
    {
        this->form->scalarTypeCombo->removeItem(i);
    }

    // Add new scalar types to list
    SortedFibers* sortedFibers = this->sortedFibersList.at(this->selectedFiberDataset);
    vtkPolyData * polydata = sortedFibers->ds->getVtkPolyData();

    // Get number of scalar types
    sortedFibers->numberOfScalarTypes = polydata->GetPointData()->GetNumberOfArrays();
    sortedFibers->hasScalars = sortedFibers->numberOfScalarTypes > 0;

    // Fill scalar list with names
    for(int i = 0; i < sortedFibers->numberOfScalarTypes; i++)
    {
        this->form->scalarTypeCombo->addItem(polydata->GetPointData()->GetArray(i)->GetName());
    }
    this->form->scalarTypeCombo->blockSignals(false);

    // Create threshold settings structs for scalar types
    if(!sortedFibers->processed)
    {
        // remove old threshold settings in case there were any
        sortedFibers->scalarThresholdSettings.clear();
        for(int i = 0; i<sortedFibers->numberOfScalarTypes; i++)
        {
            // Create struct
            ThresholdSettings* ts = new ThresholdSettings;

            // Set scalar range
            vtkDoubleArray* scalarData = static_cast<vtkDoubleArray*>(polydata->GetPointData()->GetArray(i));
            double scalarRange[2];
            scalarData->GetValueRange(scalarRange);
            memcpy(ts->scalarRange,scalarRange,sizeof(scalarRange));
            memcpy(ts->globalSetting,scalarRange,sizeof(scalarRange));
            memcpy(ts->averageScore,scalarRange,sizeof(scalarRange));
            memcpy(ts->minkowskiAverageScore,scalarRange,sizeof(scalarRange));
            ts->minkowskiOrder = 5;

            // Add to settings list
            sortedFibers->scalarThresholdSettings.append(ts);
        }
    }

    // Set output data name
    this->form->outputLineEdit->setText(sortedFibers->outputFiberDataName);

    // Compute length
    if(!sortedFibers->processed)
    {
        ComputeFiberLengthRange();

        // Set default fiber length values
        double* fiberLengthRange = sortedFibers->lengthOfFiberRange;
        sortedFibers->lengthOfFiberSetting[0] = sortedFibers->lengthOfFiberRange[0];
        sortedFibers->lengthOfFiberSetting[1] = sortedFibers->lengthOfFiberRange[1];

        //printf("length of fiber min: %d max: %d", sortedFibers->lengthOfFiberRange[0], sortedFibers->lengthOfFiberRange[1]);

        // Update fiber length slider
        BlockSignals();
        this->form->fiberLengthSlider->setRange(sortedFibers->lengthOfFiberSetting[0]*100,sortedFibers->lengthOfFiberSetting[1]*100);
        this->form->fiberLengthSpinBox->setRange(sortedFibers->lengthOfFiberSetting[0],sortedFibers->lengthOfFiberSetting[1]);
    }

    // Update GUI
    UpdateGUI();

    // Set processed flag
    sortedFibers->processed = true;

    // Select the standard scalar
    SelectScalarType(sortedFibers->selectedScalarType);

    // Enable GUI
    EnableGUI();
}

void ScoringTools::SelectScalarType(int index)
{
    // return if fiber is none
    if(this->selectedFiberDataset == -1)
        return;

    // Get selected scalar type data
    SortedFibers* sortedFibers = this->sortedFibersList.at(this->selectedFiberDataset);

     // Update selected scalar type
    sortedFibers->selectedScalarType = index;

    // Check if amount of scalar types is larger than 0
    if(sortedFibers->numberOfScalarTypes == 0)
        return;

    // Get scalar data
    vtkPolyData * polydata = sortedFibers->ds->getVtkPolyData();
    vtkDoubleArray* scalarData = static_cast<vtkDoubleArray*>(polydata->GetPointData()->GetArray(index));

    // get threshold settings struct
    ThresholdSettings* thresholdSettings = GetThresholdSettings();

    // Set average value slider ranges
    BlockSignals();
    double* scalarRange = thresholdSettings->scalarRange;
    this->form->averageValueMinSlider->setRange(scalarRange[0]*100,scalarRange[1]*100);
    this->form->averageValueMinSpinBox->setRange(scalarRange[0],scalarRange[1]);
    this->form->averageValueMaxSlider->setRange(scalarRange[0]*100,scalarRange[1]*100);
    this->form->averageValueMaxSpinBox->setRange(scalarRange[0],scalarRange[1]);

    // Set minkowski average value slider ranges
    this->form->minkowskiAverageValueMinSlider->setRange(scalarRange[0]*100,scalarRange[1]*100);
    this->form->minkowskiAverageValueMinSpinBox->setRange(scalarRange[0],scalarRange[1]);
    this->form->minkowskiAverageValueMaxSlider->setRange(scalarRange[0]*100,scalarRange[1]*100);
    this->form->minkowskiAverageValueMaxSpinBox->setRange(scalarRange[0],scalarRange[1]);

    // Set global value slider ranges
    this->form->globalMinimumSlider->setRange(scalarRange[0]*100,scalarRange[1]*100);
    this->form->globalMinimumSpinBox->setRange(scalarRange[0],scalarRange[1]);
    this->form->globalMaximumSlider->setRange(scalarRange[0]*100,scalarRange[1]*100);
    this->form->globalMaximumSpinBox->setRange(scalarRange[0],scalarRange[1]);

    // Update GUI
    UpdateGUI();

    //printf("average score: %f %f\n",thresholdSettings->averageScore[0],thresholdSettings->averageScore[1]);


}

void ScoringTools::SetActiveScalars()
{
    // return if fiber is none
    if(this->selectedFiberDataset == -1)
        return;

    // Get selected scalar type data
    SortedFibers* sortedFibers = this->sortedFibersList.at(this->selectedFiberDataset);

    // Get scalar data
    vtkPolyData * polydata = sortedFibers->ds->getVtkPolyData();

    // Set active scalars of fiber dataset
    polydata->GetPointData()->SetActiveScalars(polydata->GetPointData()->GetArray(sortedFibers->selectedScalarType)->GetName());
    this->core()->data()->dataSetChanged(sortedFibers->ds);

    if(sortedFibers->ds_processed != NULL)
    {
        vtkPolyData * polydata_processed = sortedFibers->ds_processed->getVtkPolyData();
        polydata_processed->GetPointData()->SetActiveScalars(polydata->GetPointData()->GetArray(sortedFibers->selectedScalarType)->GetName());
        this->core()->data()->dataSetChanged(sortedFibers->ds_processed);
    }
}

void ScoringTools::ComputeFibers()
{
    // return if fiber is none
    if(this->selectedFiberDataset == -1)
        return;

    // Get polydata of original fibers
    SortedFibers* sortedFibers = this->sortedFibersList.at(this->selectedFiberDataset);
    //ThresholdSettings* thresholdSettings = sortedFibers->scalarThresholdSettings.at(sortedFibers->selectedScalarType);
    vtkPolyData * polydata = sortedFibers->ds->getVtkPolyData();

    // critera
    //printf("average score: %f %f\n",thresholdSettings->averageScore[0],thresholdSettings->averageScore[1]);

    // Perform fiber selection filter
    vtkFiberSelectionFilter* selectionFilter = vtkFiberSelectionFilter::New();
	selectionFilter->SetInput(polydata);
	selectionFilter->SetMaximumFiberLength(sortedFibers->lengthOfFiberSetting[1]);
	selectionFilter->SetPrunePercentage(sortedFibers->prunePercentage);
	//selectionFilter->SetThresholdSettings((QList<ThresholdSettings*>)sortedFibers->scalarThresholdSettings);
	for(int i =0; i<sortedFibers->scalarThresholdSettings.length(); i++)
	{
	    ThresholdSettings* thresholdSettings = sortedFibers->scalarThresholdSettings.at(i);
	    //selectionFilter->AddThresholdSetting(thresholdSettings->set, thresholdSettings->averageScore, thresholdSettings->globalSetting);
	    selectionFilter->AddThresholdSetting(thresholdSettings);
	}

	//selectionFilter->SetAverageScoreRange(thresholdSettings->averageScore);
	//selectionFilter->SetScalarType(sortedFibers->selectedScalarType);

	// Run the filter
	this->core()->out()->createProgressBarForAlgorithm(selectionFilter, "Fiber selection");
	selectionFilter->Update();
	this->core()->out()->deleteProgressBarForAlgorithm(selectionFilter);

    // Prepare output polydata
    vtkPolyData* outputPoly = vtkPolyData::New();

	// In case ROI cutting is enabled, run the filter
    if(this->form->selectROICombo->currentIndex() > 0)
    {
        vtkFiberROICutting* roiCutting = vtkFiberROICutting::New();
        roiCutting->SetInput(selectionFilter->GetOutput());
        roiCutting->SetROIData(this->roiDataSets.at(this->form->selectROICombo->currentIndex()-1));
        roiCutting->SetCutFibersAtROIEnds(this->form->checkCutFibersAtROI->isChecked());

        // Get the transformation matrix
        vtkObject* tfm;
        vtkMatrix4x4* transformationMatrix;
        if (sortedFibers->ds->getAttributes()->getAttribute("transformation matrix", tfm ))
        {
            transformationMatrix = vtkMatrix4x4::SafeDownCast(tfm);
            if (transformationMatrix == 0)
            {
                return;
            }
        }
        roiCutting->SetFiberTransformationMatrix(transformationMatrix);

        // Run the filter
        this->core()->out()->createProgressBarForAlgorithm(roiCutting, "ROI cutting");
        roiCutting->Update();
        this->core()->out()->deleteProgressBarForAlgorithm(roiCutting);

        outputPoly->ShallowCopy(roiCutting->GetOutput());  // disconnect from filter to prevent unwanted future updates
    }
    else
    	outputPoly->ShallowCopy(selectionFilter->GetOutput());  // disconnect from filter to prevent unwanted future updates

    // Set active scalars of output data equal to input data
    if(sortedFibers->hasScalars)
        outputPoly->GetPointData()->SetActiveScalars(polydata->GetPointData()->GetArray(sortedFibers->selectedScalarType)->GetName());

    // Construst vIST/e dataset
    data::DataSet* ds = sortedFibers->ds_processed;
	if (ds != NULL)
	{
		ds->updateData(outputPoly);
		ds->setName(sortedFibers->outputFiberDataName);

		// Fibers should be visible, and the visualization pipeline should be updated
		ds->getAttributes()->addAttribute("isVisible", 1.0);
		ds->getAttributes()->addAttribute("updatePipeline", 1.0);

		// We add this attribute to make sure that output data sets are not added to the input data sets
		ds->getAttributes()->addAttribute("isScoreThresholded", 1);

		// Copy the transformation matrix to the output
		ds->getAttributes()->copyTransformationMatrix(sortedFibers->ds);

		this->core()->data()->dataSetChanged(ds);
	}

	// Otherwise, create a new data set
	else
	{
		ds = new data::DataSet(sortedFibers->outputFiberDataName, "fibers", outputPoly);

		// Fibers should be visible, and the visualization pipeline should be updated
		ds->getAttributes()->addAttribute("isVisible", 1.0);
		ds->getAttributes()->addAttribute("updatePipeline", 1.0);

		// We add this attribute to make sure that output data sets are not added to the input data sets
		ds->getAttributes()->addAttribute("isScoreThresholded", 1);

		// Copy the transformation matrix to the output
		ds->getAttributes()->copyTransformationMatrix(sortedFibers->ds);

		this->core()->data()->addDataSet(ds);
	}

    // Update ds locally
	sortedFibers->ds_processed = ds;

	// Hide the input data set
	sortedFibers->ds->getAttributes()->addAttribute("isVisible", -1.0);
	this->core()->data()->dataSetChanged(sortedFibers->ds);

    // Update gui enablers
	EnableGUI();
}

void ScoringTools::ShowHistogram(HistogramType histType)
{
    // return if fiber is none
    if(this->selectedFiberDataset == -1)
        return;

    // Get polydata of fibers
    SortedFibers* sortedFibers = this->sortedFibersList.at(this->selectedFiberDataset);

    vtkPolyData * polydata;
    if(histType == HISTOGRAM_INPUT)
        polydata = sortedFibers->ds->getVtkPolyData();
    else if(histType == HISTOGRAM_OUTPUT)
        polydata = sortedFibers->ds_processed->getVtkPolyData();

    vtkPointData * inputPD = polydata->GetPointData();
    vtkCellArray * inputLines = polydata->GetLines();
    ThresholdSettings* thresholdSettings = GetThresholdSettings();

    // Compute scalar histogram
    const int numberOfBins = 200;
    int averageHist[numberOfBins];
    for(int i = 0; i<numberOfBins; i++)
    {
        averageHist[i] = 0;
    }

    int minkowskiAverageHist[numberOfBins];
    for(int i = 0; i<numberOfBins; i++)
    {
        minkowskiAverageHist[i] = 0;
    }

    int perPointHist[numberOfBins];
    for(int i = 0; i<numberOfBins; i++)
    {
        perPointHist[i] = 0;
    }

    double scalarNormDenominator = thresholdSettings->scalarRange[1] - thresholdSettings->scalarRange[0];
    double scalarStep = scalarNormDenominator/(double)numberOfBins;
    double xAxis[numberOfBins];
    for(int i = 0; i<numberOfBins; i++)
    {
        xAxis[i] = thresholdSettings->scalarRange[0] + scalarStep * i;
    }

    int numberOfCells = inputLines->GetNumberOfCells();
    int selectedScalarType = sortedFibers->selectedScalarType;
    // Loop through all input fibers
	for (vtkIdType lineId = 0; lineId < numberOfCells; ++lineId)
	{
        // Get the data of the current fiber
        vtkCell * currentCell = polydata->GetCell(lineId);
        int numberOfFiberPoints = currentCell->GetNumberOfPoints();

        double averageValue = 0.0;
        double minkowskiAverageValue = 0.0;

        // Loop through all points in the fiber
        for (int pointId = 0; pointId < numberOfFiberPoints; ++pointId)
        {
            // Get the point ID of the current fiber point
            vtkIdType currentPointId = currentCell->GetPointId(pointId);

            // Get the active scalar of the fiber point
            double scalar = inputPD->GetArray(selectedScalarType)->GetTuple1(currentPointId);

            // Average value of fiber
            averageValue += scalar;

            // Minkowski average value of fiber
            if(thresholdSettings->minkowskiOrder > 1)
                minkowskiAverageValue += pow(scalar, thresholdSettings->minkowskiOrder);
            else
                minkowskiAverageValue += scalar;

            // per point histogram
            int binval = floor((scalar - thresholdSettings->scalarRange[0]) / scalarNormDenominator * (double) numberOfBins + 0.5);
            perPointHist[binval] = perPointHist[binval] + 1;
        }

        // finish average value
        averageValue /= numberOfFiberPoints;

        // finish minkowski average value
        minkowskiAverageValue /= numberOfFiberPoints;
        minkowskiAverageValue = pow(minkowskiAverageValue, 1.0/(double)thresholdSettings->minkowskiOrder);

        // normalize average value and place in histogram bin
        int binval = floor((averageValue - thresholdSettings->scalarRange[0]) / scalarNormDenominator * (double) numberOfBins + 0.5);
        std::cout << binval << std::endl;
        if(binval > 0)
            averageHist[binval] = averageHist[binval] + 1;

        // normalize minkowski average value and place in histogram bin
        binval = floor((minkowskiAverageValue - thresholdSettings->scalarRange[0]) / scalarNormDenominator * (double) numberOfBins + 0.5);
        std::cout << binval << std::endl;
        if(binval > 0)
            minkowskiAverageHist[binval] = minkowskiAverageHist[binval] + 1;
	}

	// compute maximum values
	int averageHistMaximum = 0;
	int minkowskiAverageHistMaximum = 0;
	int perPointHistMaximum = 0;
    for(int i = 0; i<numberOfBins; i++)
    {
        if(averageHist[i] > averageHistMaximum)
            averageHistMaximum = averageHist[i];
        if(minkowskiAverageHist[i] > minkowskiAverageHistMaximum)
            minkowskiAverageHistMaximum = minkowskiAverageHist[i];
        if(perPointHist[i] > perPointHistMaximum)
            perPointHistMaximum = perPointHist[i];
    }

    // Create a table with histogram values in it
    VTK_CREATE(vtkTable, table);

    VTK_CREATE(vtkFloatArray, arrXAxis);
    arrXAxis->SetName("X axis");
    table->AddColumn(arrXAxis);

    VTK_CREATE(vtkFloatArray, arrAvgHist);
    arrAvgHist->SetName("AvgHist");
    table->AddColumn(arrAvgHist);

    VTK_CREATE(vtkFloatArray, arrMinkowskiAvgHist);
    arrMinkowskiAvgHist->SetName("MinkowskiAvgHist");
    table->AddColumn(arrMinkowskiAvgHist);

    VTK_CREATE(vtkFloatArray,arrPerPointHist);
    arrPerPointHist->SetName("PerPointHist");
    table->AddColumn(arrPerPointHist);

    table->SetNumberOfRows(numberOfBins);
    for (int i = 0; i < numberOfBins; i++)
    {
        table->SetValue(i,0,xAxis[i]);
        table->SetValue(i,1,(float)averageHist[i]/(float)averageHistMaximum);
        table->SetValue(i,2,(float)minkowskiAverageHist[i]/(float)minkowskiAverageHistMaximum);
        table->SetValue(i,3,(float)perPointHist[i]/(float)perPointHistMaximum);
    }

    // Set up QT histogram window
    QWidget *histogramWindow = new QWidget();
    histogramWindow->setGeometry(0, 0, 600, 600);

    QVTKWidget *qvtkWidget = new QVTKWidget(histogramWindow);
    VTK_CREATE(vtkContextView, view); // This contains a chart object
    view->SetInteractor(qvtkWidget->GetInteractor());
    qvtkWidget->SetRenderWindow(view->GetRenderWindow());

    VTK_CREATE(vtkChartXY, chart);
    chart->GetAxis(vtkAxis::LEFT)->SetTitle("");
    chart->GetAxis(vtkAxis::BOTTOM)->SetTitle(polydata->GetPointData()->GetArray(sortedFibers->selectedScalarType)->GetName());
    QString title = QString("Scalar histogram: %1").arg(sortedFibers->ds->getName());
    chart->SetTitle(title.toUtf8().constData());
    view->GetScene()->AddItem(chart);
    chart->SetShowLegend(true);

    // Add multiple line plots, setting the colors etc
    vtkPlot *line = 0;
    line = chart->AddPlot(vtkChart::LINE);
    line->SetInput(table, 0, 1);
    line->SetColor(0, 255, 0, 255);

    //line = chart->AddPlot(vtkChart::LINE);
    //line->SetInput(table, 0, 2);
    //line->SetColor(0, 0, 255, 255);

    line = chart->AddPlot(vtkChart::LINE);
    line->SetInput(table, 0, 3);
    line->SetColor(255, 0, 0, 255);

    QVBoxLayout *layout = new QVBoxLayout(histogramWindow);
    layout->addWidget(qvtkWidget);

    histogramWindow->raise();
    histogramWindow->show();
}

///
///      FIBER CALCULATIONS
///

void ScoringTools::ComputeFiberLengthRange()
{
    // return if fiber is none
    if(this->selectedFiberDataset == -1)
        return;

    // Get polydata of original fibers
    SortedFibers* sortedFibers = this->sortedFibersList.at(this->selectedFiberDataset);
    vtkPolyData * polydata = sortedFibers->ds->getVtkPolyData();
    vtkCellArray * inputLines = polydata->GetLines();

    int numberOfCells = inputLines->GetNumberOfCells();
    int maxLength = 0;
    // Loop through all input fibers
	for (vtkIdType lineId = 0; lineId < numberOfCells; ++lineId)
	{
        // Get the data of the current fiber
        vtkCell * currentCell = polydata->GetCell(lineId);
        int numberOfFiberPoints = currentCell->GetNumberOfPoints();
        if(numberOfFiberPoints > maxLength)
            maxLength = numberOfFiberPoints;
	}

	// Update max length value
	sortedFibers->lengthOfFiberRange[0] = 0;
	sortedFibers->lengthOfFiberRange[1] = maxLength;
}

///
///     GUI CALLBACKS
///

void ScoringTools::EnableGUI()
{
    if(GetSortedFibers()->hasScalars)
        this->form->scalarGroupBox->setEnabled(true);
    else
        this->form->scalarGroupBox->setEnabled(false);
    this->form->shapeGroupBox->setEnabled(true);
    this->form->updateButton->setEnabled(true);
    if(GetSortedFibers()->ds_processed != NULL)
        this->form->displayOutputHistogramButton->setEnabled(true);
    else
        this->form->displayOutputHistogramButton->setEnabled(false);
}

void ScoringTools::DisableGUI()
{
    this->form->scalarGroupBox->setEnabled(false);
    this->form->shapeGroupBox->setEnabled(false);
    this->form->updateButton->setEnabled(false);
}

void ScoringTools::BlockSignals()
{
    this->form->averageValueMinSlider->blockSignals(true);
    this->form->averageValueMinSpinBox->blockSignals(true);
    this->form->averageValueMaxSlider->blockSignals(true);
    this->form->averageValueMaxSpinBox->blockSignals(true);
    this->form->globalMinimumSlider->blockSignals(true);
    this->form->globalMinimumSpinBox->blockSignals(true);
    this->form->globalMaximumSlider->blockSignals(true);
    this->form->globalMaximumSpinBox->blockSignals(true);
    this->form->minkowskiAverageValueMinSlider->blockSignals(true);
    this->form->minkowskiAverageValueMinSpinBox->blockSignals(true);
    this->form->minkowskiAverageValueMaxSlider->blockSignals(true);
    this->form->minkowskiAverageValueMaxSpinBox->blockSignals(true);
}

void ScoringTools::AllowSignals()
{
    this->form->averageValueMinSlider->blockSignals(false);
    this->form->averageValueMinSpinBox->blockSignals(false);
    this->form->averageValueMaxSlider->blockSignals(false);
    this->form->averageValueMaxSpinBox->blockSignals(false);
    this->form->globalMinimumSlider->blockSignals(false);
    this->form->globalMinimumSpinBox->blockSignals(false);
    this->form->globalMaximumSlider->blockSignals(false);
    this->form->globalMaximumSpinBox->blockSignals(false);
    this->form->minkowskiAverageValueMinSlider->blockSignals(false);
    this->form->minkowskiAverageValueMinSpinBox->blockSignals(false);
    this->form->minkowskiAverageValueMaxSlider->blockSignals(false);
    this->form->minkowskiAverageValueMaxSpinBox->blockSignals(false);
}

SortedFibers* ScoringTools::GetSortedFibers()
{
    SortedFibers* sortedFibers = this->sortedFibersList.at(this->selectedFiberDataset);
    return sortedFibers;
}

ThresholdSettings* ScoringTools::GetThresholdSettings()
{
    SortedFibers* sortedFibers = this->sortedFibersList.at(this->selectedFiberDataset);
    ThresholdSettings* thresholdSettings = sortedFibers->scalarThresholdSettings.at(sortedFibers->selectedScalarType);
    return thresholdSettings;
}

void ScoringTools::UpdateGUI()
{
    // block signal propagation
    BlockSignals();

    // get sorted fibers
    SortedFibers* sortedFibers = this->sortedFibersList.at(this->selectedFiberDataset);

    // set fiber length GUI settings
    this->form->fiberLengthSlider->setValue(sortedFibers->lengthOfFiberSetting[1]*SLIDER_SUBSTEPS);
    this->form->fiberLengthSpinBox->setValue(sortedFibers->lengthOfFiberSetting[1]);

    if(sortedFibers->hasScalars)
    {
        // get threshold settings struct
        ThresholdSettings* thresholdSettings = GetThresholdSettings();

        // set average value GUI settings
        this->form->averageValueMinSlider->setValue(thresholdSettings->averageScore[0]*SLIDER_SUBSTEPS);
        this->form->averageValueMinSpinBox->setValue(thresholdSettings->averageScore[0]);
        this->form->averageValueMaxSlider->setValue(thresholdSettings->averageScore[1]*SLIDER_SUBSTEPS);
        this->form->averageValueMaxSpinBox->setValue(thresholdSettings->averageScore[1]);
        this->form->averageValueMinSlider->setMaximum(this->form->averageValueMaxSlider->value()-1);
        this->form->averageValueMaxSlider->setMinimum(this->form->averageValueMinSlider->value()+1);

        // minkowski average value GUI settings
        this->form->minkowskiAverageValueMinSlider->setValue(thresholdSettings->minkowskiAverageScore[0]*SLIDER_SUBSTEPS);
        this->form->minkowskiAverageValueMinSpinBox->setValue(thresholdSettings->minkowskiAverageScore[0]);
        this->form->minkowskiAverageValueMaxSlider->setValue(thresholdSettings->minkowskiAverageScore[1]*SLIDER_SUBSTEPS);
        this->form->minkowskiAverageValueMaxSpinBox->setValue(thresholdSettings->minkowskiAverageScore[1]);
        this->form->minkowskiAverageValueMinSlider->setMaximum(this->form->minkowskiAverageValueMaxSlider->value()-1);
        this->form->minkowskiAverageValueMaxSlider->setMinimum(this->form->minkowskiAverageValueMinSlider->value()+1);
        this->form->minkowskiOrderSpinBox->setValue(thresholdSettings->minkowskiOrder);

        // global minimum
        this->form->globalMinimumSlider->setValue(thresholdSettings->globalSetting[0]*SLIDER_SUBSTEPS);
        this->form->globalMaximumSlider->setValue(thresholdSettings->globalSetting[1]*SLIDER_SUBSTEPS);
        this->form->globalMinimumSpinBox->setValue(thresholdSettings->globalSetting[0]);
        this->form->globalMaximumSpinBox->setValue(thresholdSettings->globalSetting[1]);

        // set scalar combobox
        this->form->scalarTypeCombo->setCurrentIndex(sortedFibers->selectedScalarType);

        // prune percentage
        this->form->prunePercentageSpinBox->setValue(sortedFibers->prunePercentage);
    }

    // re-enable signals
    AllowSignals();
}

//
//  slots
//

void ScoringTools::fibersComboChanged(int index)
{
    SelectFiberDataSet(index);
}

void ScoringTools::scalarTypeComboChanged(int index)
{
    SelectScalarType(index);
}

void ScoringTools::averageValueMinSliderChanged(int value)
{
    GetThresholdSettings()->averageScore[0] = (double)value / (double)SLIDER_SUBSTEPS;
    UpdateGUI();
}

void ScoringTools::averageValueMaxSliderChanged(int value)
{
    GetThresholdSettings()->averageScore[1] = (double)value / (double)SLIDER_SUBSTEPS;
    UpdateGUI();
}

void ScoringTools::averageValueMinSpinBoxChanged(double value)
{
    GetThresholdSettings()->averageScore[0] = value;
    UpdateGUI();
}

void ScoringTools::averageValueMaxSpinBoxChanged(double value)
{
    GetThresholdSettings()->averageScore[1] = value;
    UpdateGUI();
}

void ScoringTools::fiberLengthSliderChanged(int value)
{
    GetSortedFibers()->lengthOfFiberSetting[1] = (double)value / (double)SLIDER_SUBSTEPS;
    UpdateGUI();
}

void ScoringTools::fiberLengthSpinBoxChanged(double value)
{
    GetSortedFibers()->lengthOfFiberSetting[1] = value;
    UpdateGUI();
}

void ScoringTools::globalMinimumSliderChanged(int value)
{
    GetThresholdSettings()->globalSetting[0] = (double)value / (double)SLIDER_SUBSTEPS;
    UpdateGUI();
}

void ScoringTools::globalMinimumSpinBoxChanged(double value)
{
    GetThresholdSettings()->globalSetting[0] = value;
    UpdateGUI();
}

void ScoringTools::globalMaximumSliderChanged(int value)
{
    GetThresholdSettings()->globalSetting[1] = (double)value / (double)SLIDER_SUBSTEPS;
    UpdateGUI();
}

void ScoringTools::globalMaximumSpinBoxChanged(double value)
{
    GetThresholdSettings()->globalSetting[1] = value;
    UpdateGUI();
}

void ScoringTools::minkowskiAverageValueMinSliderChanged(int value)
{
    GetThresholdSettings()->minkowskiAverageScore[0] = (double)value / (double)SLIDER_SUBSTEPS;
    UpdateGUI();
}

void ScoringTools::minkowskiAverageValueMaxSliderChanged(int value)
{
    GetThresholdSettings()->minkowskiAverageScore[1] = (double)value / (double)SLIDER_SUBSTEPS;
    UpdateGUI();
}

void ScoringTools::minkowskiAverageValueMinSpinBoxChanged(double value)
{
    GetThresholdSettings()->minkowskiAverageScore[0] = value;
    UpdateGUI();
}

void ScoringTools::minkowskiAverageValueMaxSpinBoxChanged(double value)
{
    GetThresholdSettings()->minkowskiAverageScore[1] = value;
    UpdateGUI();
}

void ScoringTools::minkowskiOrderSpinBoxChanged(int value)
{
    GetThresholdSettings()->minkowskiOrder = value;
    UpdateGUI();
}

void ScoringTools::prunePercentageSpinBoxChanged(int value)
{
    GetSortedFibers()->prunePercentage = value;
    UpdateGUI();
}

void ScoringTools::updateButtonClicked()
{
    ComputeFibers();
}

void ScoringTools::displayHistogramButtonClicked()
{
    ShowHistogram(HISTOGRAM_INPUT);
}

void ScoringTools::displayOutputHistogramButtonClicked()
{
    ShowHistogram(HISTOGRAM_OUTPUT);
}

void ScoringTools::setActiveScalarsButtonClicked()
{
    SetActiveScalars();
}

void ScoringTools::outputLineEditChanged(QString text)
{
    // return if fiber is none
    if(this->selectedFiberDataset == -1)
        return;

    SortedFibers* sortedFibers = this->sortedFibersList.at(this->selectedFiberDataset);
    sortedFibers->outputFiberDataName = text;
}

///
///      GUI CONTROLS
///

//------------------------[ Connect Qt elements ]-----------------------\\

void ScoringTools::connectAll()
{
    connect(this->form->fibersCombo,SIGNAL(currentIndexChanged(int)),this,SLOT(fibersComboChanged(int)));
    connect(this->form->scalarTypeCombo,SIGNAL(currentIndexChanged(int)),this,SLOT(scalarTypeComboChanged(int)));
    connect(this->form->averageValueMinSlider,SIGNAL(valueChanged(int)),this,SLOT(averageValueMinSliderChanged(int)));
    connect(this->form->averageValueMinSpinBox,SIGNAL(valueChanged(double)),this,SLOT(averageValueMinSpinBoxChanged(double)));
    connect(this->form->averageValueMaxSlider,SIGNAL(valueChanged(int)),this,SLOT(averageValueMaxSliderChanged(int)));
    connect(this->form->averageValueMaxSpinBox,SIGNAL(valueChanged(double)),this,SLOT(averageValueMaxSpinBoxChanged(double)));
    connect(this->form->updateButton,SIGNAL(clicked()),this,SLOT(updateButtonClicked()));
    connect(this->form->displayHistogramButton,SIGNAL(clicked()),this,SLOT(displayHistogramButtonClicked()));
    connect(this->form->displayOutputHistogramButton,SIGNAL(clicked()),this,SLOT(displayOutputHistogramButtonClicked()));
    connect(this->form->setActiveScalarsButton,SIGNAL(clicked()),this,SLOT(setActiveScalarsButtonClicked()));
    connect(this->form->outputLineEdit,SIGNAL(textChanged(QString)),this,SLOT(outputLineEditChanged(QString)));
    connect(this->form->fiberLengthSlider,SIGNAL(valueChanged(int)),this,SLOT(fiberLengthSliderChanged(int)));
    connect(this->form->fiberLengthSpinBox,SIGNAL(valueChanged(double)),this,SLOT(fiberLengthSpinBoxChanged(double)));
    connect(this->form->globalMinimumSlider,SIGNAL(valueChanged(int)),this,SLOT(globalMinimumSliderChanged(int)));
    connect(this->form->globalMinimumSpinBox,SIGNAL(valueChanged(double)),this,SLOT(globalMinimumSpinBoxChanged(double)));
    connect(this->form->globalMaximumSlider,SIGNAL(valueChanged(int)),this,SLOT(globalMaximumSliderChanged(int)));
    connect(this->form->globalMaximumSpinBox,SIGNAL(valueChanged(double)),this,SLOT(globalMaximumSpinBoxChanged(double)));
    connect(this->form->minkowskiAverageValueMinSlider,SIGNAL(valueChanged(int)),this,SLOT(minkowskiAverageValueMinSliderChanged(int)));
    connect(this->form->minkowskiAverageValueMinSpinBox,SIGNAL(valueChanged(double)),this,SLOT(minkowskiAverageValueMinSpinBoxChanged(double)));
    connect(this->form->minkowskiAverageValueMaxSlider,SIGNAL(valueChanged(int)),this,SLOT(minkowskiAverageValueMaxSliderChanged(int)));
    connect(this->form->minkowskiAverageValueMaxSpinBox,SIGNAL(valueChanged(double)),this,SLOT(minkowskiAverageValueMaxSpinBoxChanged(double)));
    connect(this->form->minkowskiOrderSpinBox,SIGNAL(valueChanged(int)),this,SLOT(minkowskiOrderSpinBoxChanged(int)));
    connect(this->form->prunePercentageSpinBox,SIGNAL(valueChanged(int)),this,SLOT(prunePercentageSpinBoxChanged(int)));
}

//------------------------[ Disconnect Qt elements ]-----------------------\\

void ScoringTools::disconnectAll()
{

}



///
///     vISTe communication
///

//-----------[ Returns visualization component as VTK object ]---------------\\
//
vtkProp * ScoringTools::getVtkProp()
{
    return this->assembly;
}

//-----------------[ Returns GUI component as Qt widget ]---------------\\
//
QWidget * ScoringTools::getGUI()
{
    return this->widget;
}

}

Q_EXPORT_PLUGIN2(libScoringTools, bmia::ScoringTools)
