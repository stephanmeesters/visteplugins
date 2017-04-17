#include "ScoringMeasures.h"

#define VTK_CREATE(type, name) \
  vtkSmartPointer<type> name = vtkSmartPointer<type>::New()

#define SLIDER_SUBSTEPS 100

namespace bmia
{

///
///      INITIALIZATION
///

//------------------------[ Plugin constructor ]-----------------------\\

ScoringMeasures::ScoringMeasures() : plugin::AdvancedPlugin("ScoringMeasures")
{
    this->widget = NULL;
    this->form   = NULL;
}

//------------------------[ Plugin destructor ]-----------------------\\

ScoringMeasures::~ScoringMeasures()
{
    delete this->widget;
    delete this->form;
}

//------------------------[ Initialization ]-----------------------\\

void ScoringMeasures::init()
{
    this->widget = new QWidget();
    this->form = new Ui::ScoringMeasuresForm();
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

void ScoringMeasures::dataSetAdded(data::DataSet * d)
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
        if (d->getAttributes()->hasIntAttribute("isSM"))
			return;

        // Create new fiber struct
        SortedFibers* sortedFibers = new SortedFibers;

        // Initialize struct
        sortedFibers->ds = d;
        sortedFibers->ds_processed = NULL;
		sortedFibers->outputFiberDataName = d->getName().append("_[SM]");
		sortedFibers->processed = false;
		sortedFibers->selectedGlyphData = -1;
		sortedFibers->preprocessedPolyData = NULL;

		// Create parameter settings struct
		ParameterSettings* ps = new ParameterSettings;
		ps->useGlyphData = false;
		ps->lambda = 0.0;
		ps->beta = 0.0;
		ps->muu = 0.0;
		ps->typeOfCurve = vtkFiberScoringMeasuresFilter::CURVE_TYPE_GEODESIC;
		ps->standardizeScalars = true;
		ps->normalizeGlyphData = false;
		ps->applyLog = true;
		sortedFibers->ps = ps;

        // Add the new data set to the list of currently available fiber sets
        this->sortedFibersList.append(sortedFibers);

        // Add to UI combobox for distance measurements to fibers
        this->form->fibersCombo->addItem(d->getName());

        // If first fiber set, select by default
        if(this->sortedFibersList.count() == 1)
            SelectFiberDataSet(0);
	}

	// Discrete sphere functions
	else if (kind == "discrete sphere")
	{
	    // Check if dataset is already added
	    if(this->glyphDataSets.contains(d))
            return;

		// Check if the data set contains an image
		vtkImageData * image = d->getVtkImageData();

		if (!image)
			return;

		// Check if the image contains point data
		vtkPointData * imagePD = image->GetPointData();

		if (!imagePD)
			return;

		// Check if the point data contains a spherical directions array
		if (!(imagePD->GetArray("Spherical Directions")))
			return;

		// We can use this data set, so add it to the list and the GUI
		this->glyphDataSets.append(d);
		this->form->glyphDataCombo->addItem(d->getName());
	}
}

//------------------------[ Dataset changed ]-----------------------\\

void ScoringMeasures::dataSetChanged(data::DataSet * d)
{
    // Assert the data set pointer (should never be NULL)
    Q_ASSERT(d);

	// Get the kind of the data set
    QString kind = d->getKind();

    // to-do

}

//------------------------[ Dataset removed ]-----------------------\\

void ScoringMeasures::dataSetRemoved(data::DataSet * d)
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

int ScoringMeasures::FindInputDataSet(data::DataSet * ds)
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

void ScoringMeasures::SelectFiberDataSet(int index)
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

void ScoringMeasures::SelectGlyphDataSet(int index)
{
    // return if fiber is none
    if(this->selectedFiberDataset == -1)
        return;

    // Get selected scalar type data
    SortedFibers* sortedFibers = this->sortedFibersList.at(this->selectedFiberDataset);

    // Update selected scalar type
    sortedFibers->selectedGlyphData = index - 1;


}

void ScoringMeasures::ComputeScore()
{
    // return if fiber is none
    if(this->selectedFiberDataset == -1)
        return;

    // Get fiber info struct
    SortedFibers* sortedFibers = this->sortedFibersList.at(this->selectedFiberDataset);

    // check if glyph dataset is none
    bool useGlyphData = sortedFibers->selectedGlyphData != -1;
    bool useInternalEnergy = sortedFibers->ps->lambda != 0 || !useGlyphData;

    // update in parameter settings if no glyph data is used
    sortedFibers->ps->useGlyphData = useGlyphData;

    // Get polydata of original fibers
    vtkPolyData * polydata = sortedFibers->ds->getVtkPolyData();
    vtkImageData * image;
    if(useGlyphData)
    {
        data::DataSet * DSF = glyphDataSets.at(sortedFibers->selectedGlyphData);
        image = DSF->getVtkImageData();
    }
    else
        image = NULL;

    // Perform fiber preprocessing if not already done
    vtkPolyData* preprocessedPolyData = sortedFibers->preprocessedPolyData;
    if(useInternalEnergy && preprocessedPolyData == NULL)
    {
        // Smooth the lines
        vtkSmoothPolyDataFilter* smoothPoly = vtkSmoothPolyDataFilter::New();
        smoothPoly->SetInput(polydata);
        smoothPoly->SetNumberOfIterations(50);
        this->core()->out()->createProgressBarForAlgorithm(smoothPoly, "Scoring measure", "Smoothing the fibers...");
        smoothPoly->Update();
        this->core()->out()->deleteProgressBarForAlgorithm(smoothPoly);

        // Spline sample
        vtkSplineFilter* splineFilter = vtkSplineFilter::New();
        splineFilter->SetInput(smoothPoly->GetOutput());
        splineFilter->SetSubdivideToLength();
        splineFilter->SetLength(1);
        this->core()->out()->createProgressBarForAlgorithm(splineFilter, "Scoring measure", "Spline sampling the fibers...");
        splineFilter->Update();
        this->core()->out()->deleteProgressBarForAlgorithm(splineFilter);

        sortedFibers->preprocessedPolyData = preprocessedPolyData = splineFilter->GetOutput();
    }

    // Create filter
    vtkFiberScoringMeasuresFilter* scoringFilter = vtkFiberScoringMeasuresFilter::New();
    if(useInternalEnergy)
        scoringFilter->SetInput(preprocessedPolyData);
    else
        scoringFilter->SetInput(polydata);
    scoringFilter->SetInputVolume(image);
    scoringFilter->SetParameters(sortedFibers->ps);

    // Run the filter
	this->core()->out()->createProgressBarForAlgorithm(scoringFilter, "Fiber scoring");
	scoringFilter->Update();
	this->core()->out()->deleteProgressBarForAlgorithm(scoringFilter);

    // Copy the created polydata
	vtkPolyData* outputPoly = vtkPolyData::New();
	outputPoly->ShallowCopy(scoringFilter->GetOutput());  // disconnect from filter to prevent unwanted future updates

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
		ds->getAttributes()->addAttribute("isSM", 1);

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
	//deci->Delete();
    scoringFilter->Delete();
}

///
///     GUI CALLBACKS
///

void ScoringMeasures::EnableGUI()
{
    if(GetSortedFibers()->selectedGlyphData != -1)
    {
        this->form->lambdaSlider->setEnabled(true);
        this->form->lambdaSpinBox->setEnabled(true);
        this->form->lambdaLabel->setEnabled(true);
        this->form->lambdaTopLabel->setEnabled(true);
    }
    else
    {
        this->form->lambdaSlider->setEnabled(false);
        this->form->lambdaSpinBox->setEnabled(false);
        this->form->lambdaLabel->setEnabled(false);
        this->form->lambdaTopLabel->setEnabled(false);
    }
    this->form->dataDependentGroupBox->setEnabled(true);
    this->form->dataIndependentGroupBox->setEnabled(true);
    this->form->updateButton->setEnabled(true);
    this->form->outputLineEdit->setEnabled(true);

    this->form->settingsGroupBox->setEnabled(true);
}

void ScoringMeasures::DisableGUI()
{
    this->form->dataDependentGroupBox->setEnabled(false);
    this->form->dataIndependentGroupBox->setEnabled(false);
    this->form->updateButton->setEnabled(false);
    this->form->outputLineEdit->setEnabled(false);

    this->form->settingsGroupBox->setEnabled(false);

}

void ScoringMeasures::BlockSignals()
{
    this->form->glyphDataCombo->blockSignals(true);
}

void ScoringMeasures::AllowSignals()
{
    this->form->glyphDataCombo->blockSignals(false);
}

void ScoringMeasures::UpdateGUI()
{
    // block signal propagation
    BlockSignals();

    // get sorted fibers
    SortedFibers* sortedFibers = this->sortedFibersList.at(this->selectedFiberDataset);
    ParameterSettings* ps = sortedFibers->ps;

    // set GUI values
    this->form->lambdaSlider->setValue(ps->lambda*SLIDER_SUBSTEPS);
    this->form->lambdaSpinBox->setValue(ps->lambda);
    this->form->betaSlider->setValue(ps->beta*SLIDER_SUBSTEPS);
    this->form->betaSpinBox->setValue(ps->beta);
    this->form->muuSlider->setValue(ps->muu*SLIDER_SUBSTEPS);
    this->form->muuSpinBox->setValue(ps->muu);

    this->form->glyphDataCombo->setCurrentIndex(sortedFibers->selectedGlyphData + 1);

    this->form->standardizeScalarsCheckBox->setChecked(ps->standardizeScalars);
    this->form->normalizeGlyphDataCheckBox->setChecked(ps->normalizeGlyphData);
    this->form->applyLogarithmCheckBox->setChecked(ps->applyLog);

    // re-enable signals
    AllowSignals();
}

SortedFibers* ScoringMeasures::GetSortedFibers()
{
    SortedFibers* sortedFibers = this->sortedFibersList.at(this->selectedFiberDataset);
    return sortedFibers;
}

ParameterSettings* ScoringMeasures::GetParameterSettings()
{
    SortedFibers* sortedFibers = this->sortedFibersList.at(this->selectedFiberDataset);
    return sortedFibers->ps;
}

///
///     SLOTS
///

void ScoringMeasures::fibersComboChanged(int index)
{
    SelectFiberDataSet(index);
}

void ScoringMeasures::glyphDataComboChanged(int index)
{
    SelectGlyphDataSet(index);
    EnableGUI();
}

void ScoringMeasures::updateButtonClicked()
{
    ComputeScore();
}

void ScoringMeasures::standardizeScalarsCheckBoxChanged(bool checked)
{
    GetSortedFibers()->ps->standardizeScalars = checked;
}

void ScoringMeasures::normalizeGlyphDataCheckBoxChanged(bool checked)
{
    GetSortedFibers()->ps->normalizeGlyphData = checked;
}

void ScoringMeasures::applyLogarithmCheckBoxChanged(bool checked)
{
    GetSortedFibers()->ps->applyLog = checked;
}

void ScoringMeasures::lambdaSliderChanged(int value)
{
    GetParameterSettings()->lambda = (double)value / (double)SLIDER_SUBSTEPS;
    UpdateGUI();
}

void ScoringMeasures::lambdaSpinBoxChanged(double value)
{
    GetParameterSettings()->lambda = value;
    UpdateGUI();
}

void ScoringMeasures::betaSliderChanged(int value)
{
    GetParameterSettings()->beta = (double)value / (double)SLIDER_SUBSTEPS;
    UpdateGUI();
}

void ScoringMeasures::betaSpinBoxChanged(double value)
{
    GetParameterSettings()->beta = value;
    UpdateGUI();
}

void ScoringMeasures::muuSliderChanged(int value)
{
    GetParameterSettings()->muu = (double)value / (double)SLIDER_SUBSTEPS;
    UpdateGUI();
}

void ScoringMeasures::muuSpinBoxChanged(double value)
{
    GetParameterSettings()->muu = value;
    UpdateGUI();
}

///
///      GUI CONTROLS
///

//------------------------[ Connect Qt elements ]-----------------------\\

void ScoringMeasures::connectAll()
{
    connect(this->form->fibersCombo,SIGNAL(currentIndexChanged(int)),this,SLOT(fibersComboChanged(int)));
    connect(this->form->glyphDataCombo,SIGNAL(currentIndexChanged(int)),this,SLOT(glyphDataComboChanged(int)));
    connect(this->form->updateButton,SIGNAL(clicked()),this,SLOT(updateButtonClicked()));

    connect(this->form->standardizeScalarsCheckBox,SIGNAL(toggled(bool)),this,SLOT(standardizeScalarsCheckBoxChanged(bool)));
    connect(this->form->normalizeGlyphDataCheckBox,SIGNAL(toggled(bool)),this,SLOT(normalizeGlyphDataCheckBoxChanged(bool)));
    connect(this->form->applyLogarithmCheckBox,SIGNAL(toggled(bool)),this,SLOT(applyLogarithmCheckBoxChanged(bool)));

    connect(this->form->lambdaSlider,SIGNAL(valueChanged(int)),this,SLOT(lambdaSliderChanged(int)));
    connect(this->form->lambdaSpinBox,SIGNAL(valueChanged(double)),this,SLOT(lambdaSpinBoxChanged(double)));
    connect(this->form->betaSlider,SIGNAL(valueChanged(int)),this,SLOT(betaSliderChanged(int)));
    connect(this->form->betaSpinBox,SIGNAL(valueChanged(double)),this,SLOT(betaSpinBoxChanged(double)));
    connect(this->form->muuSlider,SIGNAL(valueChanged(int)),this,SLOT(muuSliderChanged(int)));
    connect(this->form->muuSpinBox,SIGNAL(valueChanged(double)),this,SLOT(muuSpinBoxChanged(double)));
}

//------------------------[ Disconnect Qt elements ]-----------------------\\

void ScoringMeasures::disconnectAll()
{

}



///
///     vISTe communication
///

//-----------[ Returns visualization component as VTK object ]---------------\\
//
vtkProp * ScoringMeasures::getVtkProp()
{
    return this->assembly;
}

//-----------------[ Returns GUI component as Qt widget ]---------------\\
//
QWidget * ScoringMeasures::getGUI()
{
    return this->widget;
}

}

Q_EXPORT_PLUGIN2(libScoringMeasures, bmia::ScoringMeasures)
