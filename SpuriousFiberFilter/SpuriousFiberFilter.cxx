#include "SpuriousFiberFilter.h"

#define VTK_CREATE(type, name) \
  vtkSmartPointer<type> name = vtkSmartPointer<type>::New()

#define SLIDER_SUBSTEPS 100

#include <iostream>

namespace bmia
{

///
///      INITIALIZATION
///

//------------------------[ Plugin constructor ]-----------------------\\

SpuriousFiberFilter::SpuriousFiberFilter() : plugin::AdvancedPlugin("SpuriousFiberFilter")
{
    this->widget = NULL;
    this->form   = NULL;
}

//------------------------[ Plugin destructor ]-----------------------\\

SpuriousFiberFilter::~SpuriousFiberFilter()
{
    delete this->widget;
    delete this->form;
}

//------------------------[ Initialization ]-----------------------\\

void SpuriousFiberFilter::init()
{
    this->widget = new QWidget();
    this->form = new Ui::SpuriousFiberFilterForm();
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

void SpuriousFiberFilter::dataSetAdded(data::DataSet * d)
{
    // Assert the data set pointer (should never be NULL)
    Q_ASSERT(d);

	// Get the kind of the data set
    QString kind = d->getKind();

    // Load fiber dataset
    if (kind == "fibers")
	{
        // Check if dataset is already added
        if(FindInputDataSet(d) != -1)
            return;

	    // Check if fiber has polydata
	    if (d->getVtkPolyData() == NULL)
			return;

        // Check if the fiber is created by this plugin
        if (d->getAttributes()->hasIntAttribute("isSpuriousFiltered"))
			return;

        // Create new fiber struct
        SortedFibers* sortedFibers = new SortedFibers;

        // Initialize struct
        sortedFibers->ds = d;
        sortedFibers->ds_processed = NULL;
		sortedFibers->outputFiberDataName = d->getName().append("_[SF]");
		sortedFibers->processed = false;

		// Create parameter settings struct
		ParameterSettings* ps = new ParameterSettings;
		ps->applyKernelInBothDirs = false;
		ps->D33 = 1.0;
		ps->D44 = 0.04;
		ps->t = 1.4;
		ps->windowSize = 8;
		ps->epsilon = 0;
		ps->minDist = 15;
		sortedFibers->ps = ps;
		ps->requireRecompute = true;
		ps->cutoff = 3.5;

		ps->fiberMinScores = NULL;
        ps->fiberScores = NULL;
        ps->fiberStartIds = NULL;
        ps->avgScoreTotal = 0.0;

        ps->applySubsampling = false;
        ps->samplingStep = 1.0;
        ps->goodness = -1;

        ps->applySelectAnterior = false;
        ps->numberOfAnteriorFibers = 3000;

        // Add the new data set to the list of currently available fiber sets
        this->sortedFibersList.append(sortedFibers);

        // Add to UI combobox for distance measurements to fibers
        this->form->fibersCombo->addItem(d->getName());

        // If first fiber set, select by default
        if(this->sortedFibersList.count() == 1)
            SelectFiberDataSet(0);
	}
}

//------------------------[ Dataset changed ]-----------------------\\

void SpuriousFiberFilter::dataSetChanged(data::DataSet * d)
{
    // Assert the data set pointer (should never be NULL)
    Q_ASSERT(d);

	// Get the kind of the data set
    QString kind = d->getKind();

    // to-do

}

//------------------------[ Dataset removed ]-----------------------\\

void SpuriousFiberFilter::dataSetRemoved(data::DataSet * d)
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
        sortedFibers->preprocessedPolyData = NULL;

        // Remove from collection
        this->sortedFibersList.removeAt(dsIndex);
	}
}

int SpuriousFiberFilter::FindInputDataSet(data::DataSet * ds)
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

void SpuriousFiberFilter::SelectFiberDataSet(int index)
{
    // set selected fiber index
    this->selectedFiberDataset = index - 1;

    // no fiber selected
    if(this->selectedFiberDataset == -1)
    {
        DisableGUI();
        return;
    }

    // Set output data name
    SortedFibers* sortedFibers = this->sortedFibersList.at(this->selectedFiberDataset);
    this->form->outputLineEdit->setText(sortedFibers->outputFiberDataName);

    // Update GUI
    UpdateGUI();

    // Enable GUI
    EnableGUI();
}

void SpuriousFiberFilter::ComputeScore()
{
    // return if fiber is none
    if(this->selectedFiberDataset == -1)
        return;

    // print
    printf("Starting processing of fibers with id:%d\n", this->selectedFiberDataset);

    // Get fiber info struct
    SortedFibers* sortedFibers = this->sortedFibersList.at(this->selectedFiberDataset);
    ParameterSettings* ps = sortedFibers->ps;

    // set name
    ps->outputFiberDataName = sortedFibers->outputFiberDataName;

    // Get polydata of original fibers
    vtkPolyData * polydata = sortedFibers->ds->getVtkPolyData();

    // Select anterior fibers
    vtkFiberSelectAnterior* anteriorFilter;
    if(ps->applySelectAnterior)
    {
        ps->numberOfAnteriorFibers = std::min(ps->numberOfAnteriorFibers,(int)polydata->GetLines()->GetNumberOfCells());
        printf("NUMBER OF ANTEROR FIBERS: %d\n",ps->numberOfAnteriorFibers);

        anteriorFilter = vtkFiberSelectAnterior::New();
        anteriorFilter->SetInput(polydata);
        anteriorFilter->SetNumberOfAnteriorFibers(ps->numberOfAnteriorFibers);

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
        anteriorFilter->SetFiberTransformationMatrix(transformationMatrix);


        this->core()->out()->createProgressBarForAlgorithm(anteriorFilter, "Spurious fiber filter", "Selecting anterior fibers...");
        anteriorFilter->Update();
        this->core()->out()->deleteProgressBarForAlgorithm(anteriorFilter);
       // preprocessedPolyData = anteriorFilter->GetOutput();
       //
    }

    // Spline sample
    vtkSplineFilter* splineFilter;
    if(ps->applySubsampling)
    {
        splineFilter = vtkSplineFilter::New();
        if(ps->applySelectAnterior)
            splineFilter->SetInput(anteriorFilter->GetOutput());
        else
            splineFilter->SetInput(polydata);
        splineFilter->SetSubdivideToLength();
        splineFilter->SetLength(1);
        this->core()->out()->createProgressBarForAlgorithm(splineFilter, "Spurious fiber filter", "Spline sampling the fibers...");
        splineFilter->Update();
        this->core()->out()->deleteProgressBarForAlgorithm(splineFilter);

        //sortedFibers->preprocessedPolyData = preprocessedPolyData = splineFilter->GetOutput();
       // preprocessedPolyData = splineFilter->GetOutput();
        //splineFilter->Delete();
    }

    // Create filter
    vtkFiberSpuriousFilter* scoringFilter = vtkFiberSpuriousFilter::New();
    if(ps->applySubsampling)
        scoringFilter->SetInput(splineFilter->GetOutput());
    else if(ps->applySelectAnterior)
        scoringFilter->SetInput(anteriorFilter->GetOutput());
    else
        scoringFilter->SetInput(polydata);
    scoringFilter->SetParameters(sortedFibers->ps);

    // Run the filter
	this->core()->out()->createProgressBarForAlgorithm(scoringFilter, "Spurious fiber filter");
	scoringFilter->Update();
	this->core()->out()->deleteProgressBarForAlgorithm(scoringFilter);

    // Copy the created polydata
	vtkPolyData* outputPoly = vtkPolyData::New();
	outputPoly->ShallowCopy(scoringFilter->GetOutput());  // disconnect from filter to prevent unwanted future updates



		// Create a polydata writer
		//vtkPolyDataWriter * writer = vtkPolyDataWriter::New();

  //      QFileInfo asdd(ps->outputFiberDataName);
  //      QString asd = QString("/home/linux/Stephan/FIBER_TO_BUNDLE_COHERENCE/output/%1.vtk").arg(asdd.fileName());
  //      qDebug() << asd;

		//// Configure the writer
		//writer->SetFileName(asd.toLocal8Bit().constData());
		//writer->SetInput(outputPoly);
		//writer->SetFileTypeToASCII();

		//// Enable progress bar for the writer
		//this->core()->out()->createProgressBarForAlgorithm(writer, "Fiber Visualization", "Writing fibers to file...");

		//// Write output file
		//writer->Write();

		//// Disable progress bar for the writer
		//this->core()->out()->deleteProgressBarForAlgorithm(writer);

		//// Delete the writer
		//writer->Delete();



    // Construst vIST/e dataset
    data::DataSet* ds = sortedFibers->ds_processed;
	if (ds != NULL)
	{
		ds->updateData(outputPoly);
		ds->setName(sortedFibers->outputFiberDataName);

		// Fibers should be visible, and the visualization pipeline should be updated
		ds->getAttributes()->addAttribute("isVisible", 1.0);
		ds->getAttributes()->addAttribute("updatePipeline", 1.0);

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
		ds->getAttributes()->addAttribute("isSpuriousFiltered", 1);

		// Copy the transformation matrix to the output
		ds->getAttributes()->copyTransformationMatrix(sortedFibers->ds);

		this->core()->data()->addDataSet(ds);
	}

    // Update ds locally
	sortedFibers->ds_processed = ds;

	// Hide the input data set
	sortedFibers->ds->getAttributes()->addAttribute("isVisible", -1.0);
	this->core()->data()->dataSetChanged(sortedFibers->ds);

	// Clean up
	if(ps->applySelectAnterior)
        anteriorFilter->Delete();
	if(ps->applySubsampling)
        splineFilter->Delete();
    scoringFilter->Delete();

    // update gui
    UpdateGUI();
}

///
///     GUI CALLBACKS
///

void SpuriousFiberFilter::EnableGUI()
{
    this->form->fiberPreprocGroupBox->setEnabled(true);
    this->form->kernelSettingsGroupBox->setEnabled(true);
    this->form->thresholdSettingsGroupBox->setEnabled(true);
    this->form->updateButton->setEnabled(true);
    this->form->outputLineEdit->setEnabled(true);
}

void SpuriousFiberFilter::DisableGUI()
{
    this->form->fiberPreprocGroupBox->setEnabled(false);
    this->form->kernelSettingsGroupBox->setEnabled(false);
    this->form->thresholdSettingsGroupBox->setEnabled(false);
    this->form->updateButton->setEnabled(false);
    this->form->outputLineEdit->setEnabled(false);
}

void SpuriousFiberFilter::BlockSignals()
{
    this->form->fibersCombo->blockSignals(true);
    this->form->anteriorSelectSpinBox->blockSignals(true);
}

void SpuriousFiberFilter::AllowSignals()
{
    this->form->fibersCombo->blockSignals(false);
    this->form->anteriorSelectSpinBox->blockSignals(false);
}

void SpuriousFiberFilter::UpdateGUI()
{
    // block signal propagation
    BlockSignals();

    // get sorted fibers
    SortedFibers* sortedFibers = this->sortedFibersList.at(this->selectedFiberDataset);
    ParameterSettings* ps = sortedFibers->ps;

    // set GUI values
    this->form->fibersCombo->setCurrentIndex(this->selectedFiberDataset == -1?0:this->selectedFiberDataset+1);
    this->form->D33Slider->setValue(ps->D33*SLIDER_SUBSTEPS);
    this->form->D33SpinBox->setValue(ps->D33);
    this->form->D44Slider->setValue(ps->D44*SLIDER_SUBSTEPS);
    this->form->D44SpinBox->setValue(ps->D44);
    this->form->timeSlider->setValue(ps->t*SLIDER_SUBSTEPS);
    this->form->timeSpinBox->setValue(ps->t);
    this->form->epsilonSlider->setValue(ps->epsilon*SLIDER_SUBSTEPS);
    this->form->epsilonSpinBox->setValue(ps->epsilon);
    this->form->windowSizeSlider->setValue(ps->windowSize*SLIDER_SUBSTEPS);
    this->form->windowSizeSpinBox->setValue(ps->windowSize);
    this->form->minDistSlider->setValue(ps->minDist*SLIDER_SUBSTEPS);
    this->form->minDistSpinBox->setValue(ps->minDist);
    this->form->cutoffSlider->setValue(ps->cutoff*SLIDER_SUBSTEPS);
    this->form->cutoffSpinBox->setValue(ps->cutoff);

    this->form->applyKernelInBothDirsCheckBox->setChecked(ps->applyKernelInBothDirs);

    this->form->applySubsamplingCheckBox->setChecked(ps->applySubsampling);
    this->form->subsampleSpinBox->setValue(ps->samplingStep);
    this->form->applySelectAnteriorCheckBox->setChecked(ps->applySelectAnterior);

    vtkPolyData* pd = sortedFibers->ds->getVtkPolyData();
    if(pd)
    {
        this->form->anteriorSelectSpinBox->setMaximum(pd->GetLines()->GetNumberOfCells());
    }
    this->form->anteriorSelectSpinBox->setValue(ps->numberOfAnteriorFibers);

    // show goodness result
    QString goodness;
    if(sortedFibers->ps->goodness==-1)
        goodness = QString("Measure of fiber goodness: --");
    else
        goodness = QString("Measure of fiber goodness: %1").arg(sortedFibers->ps->goodness);
    this->form->labelMeasureFiberGoodness->setText(goodness);

    // re-enable signals
    AllowSignals();
}

SortedFibers* SpuriousFiberFilter::GetSortedFibers()
{
    SortedFibers* sortedFibers = this->sortedFibersList.at(this->selectedFiberDataset);
    return sortedFibers;
}

ParameterSettings* SpuriousFiberFilter::GetParameterSettings()
{
    SortedFibers* sortedFibers = this->sortedFibersList.at(this->selectedFiberDataset);
    return sortedFibers->ps;
}

///
///     SLOTS
///

/** Fiber data **/
void SpuriousFiberFilter::fibersComboChanged(int index)
{
    SelectFiberDataSet(index);
}

/** Update button **/
void SpuriousFiberFilter::updateButtonClicked()
{
    ComputeScore();
}

/** Process batch button **/
void SpuriousFiberFilter::processBatchButtonClicked()
{
    int originalSelectedFiberData = this->selectedFiberDataset;
    for(int i = 0; i<this->sortedFibersList.length(); i++)
    {
        SelectFiberDataSet(i+1);
        GetParameterSettings()->requireRecompute = true;

        // TEST FOR MAX SCORE THRESHOLD RANGE
        SortedFibers* sortedFibers = this->sortedFibersList.at(this->selectedFiberDataset);
        double originalThres = GetParameterSettings()->cutoff;
        QString originalOutputName = sortedFibers->outputFiberDataName;
        //for(double j = 0.001; j<=0.101; j+=0.01)
        //for(double j = 0.0001; j<=0.101; j+=0.005)
        for(double j = 0.0001; j<=0.071; j+=0.01)
        {
            //j=0.101;
            GetParameterSettings()->cutoff = j*0.053474; // compensating for default kernel value
            sortedFibers->outputFiberDataName = QString("%1_%2").arg(originalOutputName).arg(j);
            UpdateGUI();
            ComputeScore();
            //break;
        }
        sortedFibers->outputFiberDataName = originalOutputName;
        GetParameterSettings()->cutoff = originalThres;

        //ComputeScore();
    }

    // copy output to dropbox
   // system("cp /home/linux/Stephan/FIBER_TO_BUNDLE_COHERENCE/output/* /home/linux/Stephan/Dropbox/newresults");
}

/** D33 **/
void SpuriousFiberFilter::D33SliderChanged(int value)
{
    GetParameterSettings()->D33 = (double)value / (double)SLIDER_SUBSTEPS;
    UpdateGUI();
    GetParameterSettings()->requireRecompute = true;
}

void SpuriousFiberFilter::D33SpinBoxChanged(double value)
{
    GetParameterSettings()->D33 = value;
    UpdateGUI();
    GetParameterSettings()->requireRecompute = true;
}

/** D44 **/
void SpuriousFiberFilter::D44SliderChanged(int value)
{
    GetParameterSettings()->D44 = (double)value / (double)SLIDER_SUBSTEPS;
    UpdateGUI();
    GetParameterSettings()->requireRecompute = true;
}

void SpuriousFiberFilter::D44SpinBoxChanged(double value)
{
    GetParameterSettings()->D44 = value;
    UpdateGUI();
    GetParameterSettings()->requireRecompute = true;
}

/** Time **/
void SpuriousFiberFilter::TimeSliderChanged(int value)
{
    GetParameterSettings()->t = (double)value / (double)SLIDER_SUBSTEPS;
    UpdateGUI();
    GetParameterSettings()->requireRecompute = true;
}

void SpuriousFiberFilter::TimeSpinBoxChanged(double value)
{
    GetParameterSettings()->t = value;
    UpdateGUI();
    GetParameterSettings()->requireRecompute = true;
}

/** Min dist **/
void SpuriousFiberFilter::minDistSliderChanged(int value)
{
    GetParameterSettings()->minDist = (double)value / (double)SLIDER_SUBSTEPS;
    UpdateGUI();
    GetParameterSettings()->requireRecompute = true;
}

void SpuriousFiberFilter::minDistSpinBoxChanged(double value)
{
    GetParameterSettings()->minDist = value;
    UpdateGUI();
    GetParameterSettings()->requireRecompute = true;
}

/** Window size **/
void SpuriousFiberFilter::WindowSizeSliderChanged(int value)
{
    GetParameterSettings()->windowSize = (double)value / (double)SLIDER_SUBSTEPS;
    UpdateGUI();
    GetParameterSettings()->requireRecompute = true;
}

void SpuriousFiberFilter::WindowSizeSpinBoxChanged(double value)
{
    GetParameterSettings()->windowSize = value;
    UpdateGUI();
    GetParameterSettings()->requireRecompute = true;
}

/** Epsilon **/
void SpuriousFiberFilter::EpsilonSliderChanged(int value)
{
    GetParameterSettings()->epsilon = (double)value / (double)SLIDER_SUBSTEPS;
    UpdateGUI();
}

void SpuriousFiberFilter::EpsilonSpinBoxChanged(double value)
{
    GetParameterSettings()->epsilon = value;
    UpdateGUI();
}

/** Apply kernel in both dirs? **/
void SpuriousFiberFilter::applyKernelInBothDirsCheckBoxChanged(bool checked)
{
    GetSortedFibers()->ps->applyKernelInBothDirs = checked;
    GetParameterSettings()->requireRecompute = true;
}

/** Apply fiber subsampling? **/
void SpuriousFiberFilter::applySubsamplingCheckBoxChanged(bool checked)
{
    GetSortedFibers()->ps->applySubsampling = checked;
    GetParameterSettings()->requireRecompute = true;
}

/** Fiber sampling step **/
void SpuriousFiberFilter::subsampleSpinBoxChanged(double value)
{
    GetSortedFibers()->ps->samplingStep = value;
    GetParameterSettings()->requireRecompute = true;
}

/** Apply select anterior fibers? **/
void SpuriousFiberFilter::applySelectAnteriorCheckBoxChanged(bool checked)
{
    GetSortedFibers()->ps->applySelectAnterior = checked;
    GetParameterSettings()->requireRecompute = true;
}

/** How many fibers are selected **/
void SpuriousFiberFilter::anteriorSelectSpinBoxChanged(int value)
{
    GetSortedFibers()->ps->numberOfAnteriorFibers = value;
    GetParameterSettings()->requireRecompute = true;
    //printf("%d\n",value);
}

/** Cut-off value (exclude closest fibers within -- mm) **/
void SpuriousFiberFilter::cutoffSliderChanged(int value)
{
    GetParameterSettings()->cutoff = (double)value / (double)SLIDER_SUBSTEPS;
    UpdateGUI();
    GetParameterSettings()->requireRecompute = true;
}

/** Cut-off value (exclude closest fibers within -- mm) **/
void SpuriousFiberFilter::cutoffSpinBoxChanged(double value)
{
    GetParameterSettings()->cutoff = value;
    UpdateGUI();
    GetParameterSettings()->requireRecompute = true;
}

/** Change fiber output name **/
void SpuriousFiberFilter::outputLineEditChanged(QString value)
{
    this->sortedFibersList.at(this->selectedFiberDataset)->outputFiberDataName = value;
}


///
///      GUI CONTROLS
///

//------------------------[ Connect Qt elements ]-----------------------\\

void SpuriousFiberFilter::connectAll()
{
    connect(this->form->fibersCombo,SIGNAL(currentIndexChanged(int)),this,SLOT(fibersComboChanged(int)));
    connect(this->form->updateButton,SIGNAL(clicked()),this,SLOT(updateButtonClicked()));
    connect(this->form->processBatchButton,SIGNAL(clicked()),this,SLOT(processBatchButtonClicked()));

    connect(this->form->applyKernelInBothDirsCheckBox,SIGNAL(toggled(bool)),this,SLOT(applyKernelInBothDirsCheckBoxChanged(bool)));

    connect(this->form->D33Slider,SIGNAL(valueChanged(int)),this,SLOT(D33SliderChanged(int)));
    connect(this->form->D33SpinBox,SIGNAL(valueChanged(double)),this,SLOT(D33SpinBoxChanged(double)));
    connect(this->form->D44Slider,SIGNAL(valueChanged(int)),this,SLOT(D44SliderChanged(int)));
    connect(this->form->D44SpinBox,SIGNAL(valueChanged(double)),this,SLOT(D44SpinBoxChanged(double)));
    connect(this->form->windowSizeSlider,SIGNAL(valueChanged(int)),this,SLOT(WindowSizeSliderChanged(int)));
    connect(this->form->windowSizeSpinBox,SIGNAL(valueChanged(double)),this,SLOT(WindowSizeSpinBoxChanged(double)));
    connect(this->form->epsilonSlider,SIGNAL(valueChanged(int)),this,SLOT(EpsilonSliderChanged(int)));
    connect(this->form->epsilonSpinBox,SIGNAL(valueChanged(double)),this,SLOT(EpsilonSpinBoxChanged(double)));
    connect(this->form->timeSlider,SIGNAL(valueChanged(int)),this,SLOT(TimeSliderChanged(int)));
    connect(this->form->timeSpinBox,SIGNAL(valueChanged(double)),this,SLOT(TimeSpinBoxChanged(double)));
    connect(this->form->minDistSlider,SIGNAL(valueChanged(int)),this,SLOT(minDistSliderChanged(int)));
    connect(this->form->minDistSpinBox,SIGNAL(valueChanged(double)),this,SLOT(minDistSpinBoxChanged(double)));
    connect(this->form->cutoffSlider,SIGNAL(valueChanged(int)),this,SLOT(cutoffSliderChanged(int)));
    connect(this->form->cutoffSpinBox,SIGNAL(valueChanged(double)),this,SLOT(cutoffSpinBoxChanged(double)));

    connect(this->form->applySubsamplingCheckBox,SIGNAL(toggled(bool)),this,SLOT(applySubsamplingCheckBoxChanged(bool)));
    connect(this->form->subsampleSpinBox,SIGNAL(valueChanged(double)),this,SLOT(subsampleSpinBoxChanged(double)));
    connect(this->form->applySelectAnteriorCheckBox,SIGNAL(toggled(bool)),this,SLOT(applySelectAnteriorCheckBoxChanged(bool)));
    connect(this->form->anteriorSelectSpinBox,SIGNAL(valueChanged(int)),this,SLOT(anteriorSelectSpinBoxChanged(int)));

    connect(this->form->outputLineEdit,SIGNAL(textEdited(QString)),this,SLOT(outputLineEditChanged(QString)));

}

//------------------------[ Disconnect Qt elements ]-----------------------\\

void SpuriousFiberFilter::disconnectAll()
{

}



///
///     vISTe communication
///

//-----------[ Returns visualization component as VTK object ]---------------\\
//
vtkProp * SpuriousFiberFilter::getVtkProp()
{
    return this->assembly;
}

//-----------------[ Returns GUI component as Qt widget ]---------------\\
//
QWidget * SpuriousFiberFilter::getGUI()
{
    return this->widget;
}

}

Q_EXPORT_PLUGIN2(libSpuriousFiberFilter, bmia::SpuriousFiberFilter)
