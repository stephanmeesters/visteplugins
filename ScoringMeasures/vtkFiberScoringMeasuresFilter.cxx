
/** Includes */

#include "vtkFiberScoringMeasuresFilter.h"


namespace bmia {


vtkStandardNewMacro(vtkFiberScoringMeasuresFilter);


//-----------------------------[ Constructor ]-----------------------------\\

vtkFiberScoringMeasuresFilter::vtkFiberScoringMeasuresFilter()
{
	// Set default options
    this->inputVolume = NULL;
}


//------------------------------[ Destructor ]-----------------------------\\

vtkFiberScoringMeasuresFilter::~vtkFiberScoringMeasuresFilter()
{

}

//-------------------------------[ SetInputVolume ]-------------------------------\\

void vtkFiberScoringMeasuresFilter::SetInputVolume(vtkImageData * image)
{
	// Store the pointer
	this->inputVolume = image;
}

//-------------------------------[ SetParameters ]-------------------------------\\

void vtkFiberScoringMeasuresFilter::SetParameters(ParameterSettings* ps)
{
    // Store the pointer
    this->ps = ps;
}

//-------------------------------[ Execute ]-------------------------------\\

double* Difference(double* vec, double* vec2)
{
    double* dp = (double*) malloc(3*sizeof(double));
    dp[0] = vec[0] - vec2[0];
    dp[1] = vec[1] - vec2[1];
    dp[2] = vec[2] - vec2[2];
    return dp;
}

double* HalvedDifference(double* vec, double* vec2)
{
    double* dp = (double*) malloc(3*sizeof(double));
    dp[0] = (vec[0] - vec2[0])/2.0;
    dp[1] = (vec[1] - vec2[1])/2.0;
    dp[2] = (vec[2] - vec2[2])/2.0;
    return dp;
}

double* Normalize(double* vec)
{
    double length = sqrt(vec[0]*vec[0] + vec[1]*vec[1] + vec[2]*vec[2]);
    double* nvec = (double*) malloc(3*sizeof(double));
    nvec[0] = vec[0]/length;
    nvec[1] = vec[1]/length;
    nvec[2] = vec[2]/length;
    return nvec;
}

double* Cross(double* vec, double* vec2)
{
    double* c = (double*) malloc(3*sizeof(double));
    c[0] = vec[1]*vec2[2] - vec[2]*vec2[1];
    c[1] = vec[2]*vec2[0] - vec[0]*vec2[2];
    c[2] = vec[0]*vec2[1] - vec[1]*vec2[0];
    return c;
}

double Norm(double* vec)
{
    return vec[0]*vec[0] + vec[1]*vec[1] + vec[2]*vec[2];
}

void PrintVector(double* vec)
{
    if(vec == NULL)
        return;
    printf("%f, %f, %f \n", vec[0], vec[1], vec[2]);
}

void vtkFiberScoringMeasuresFilter::Execute()
{
    //
    //      CHECK POLYDATA INPUT
    //

	// Get the input
	vtkPolyData * input = this->GetInput();
	if(!input)
	{
		vtkErrorMacro(<< "Input has not been set.");
		return;
	}

	// Check if the input contains point data
	vtkPointData * inputPD = input->GetPointData();
	if (!inputPD)
	{
		vtkErrorMacro(<< "Input does not have point data.");
		return;
	}

	// Get the points of the input
	vtkPoints * inputPoints = input->GetPoints();
	if (!inputPoints)
	{
		vtkErrorMacro(<< "Input does not have points.");
		return;
	}

	// Get the lines array of the input
	vtkCellArray * inputLines = input->GetLines();
	if (!inputLines)
	{
		vtkErrorMacro(<< "Input does not have lines.");
		return;
	}

	// Get the output
	vtkPolyData * output = this->GetOutput();
	if (!output)
	{
		vtkErrorMacro(<< "Output has not been set.");
		return;
	}

	// Check if the output contains point data
	vtkPointData * outputPD = output->GetPointData();
	if (!outputPD)
	{
		vtkErrorMacro(<< "Output does not have point data.");
		return;
	}

    //
    //      CHECK GLYPH DATA INPUT
    //
    vtkPointData * imagePD;
    vtkDoubleArray * radiiArray;
    vtkDoubleArray * anglesArray;
    int numberOfAngles;
	if(ps->useGlyphData)
	{
        imagePD = this->inputVolume->GetPointData();
        if (!imagePD)
        {
            vtkErrorMacro(<< "Point data in glyph data have not been set.");
            return;
        }

        // Get the array containing the glyphs radii
        radiiArray = vtkDoubleArray::SafeDownCast(imagePD->GetArray("Vectors"));
        if (!radiiArray)
        {
            vtkErrorMacro(<< "Vectors in glyph data have not been set.");
            return;
        }

        // Get the array containing the angles for each glyph vertex
        anglesArray = vtkDoubleArray::SafeDownCast(imagePD->GetArray("Spherical Directions"));
        if (!anglesArray)
        {
            vtkErrorMacro(<< "Spherical directions in glyph data have not been set.");
            return;
        }

        // Angles array should have two components. Furthermore, the number of sets of
        // angles should match the number of radii.
        if (anglesArray->GetNumberOfComponents() != 2 || anglesArray->GetNumberOfTuples() != radiiArray->GetNumberOfComponents())
        {
            vtkErrorMacro(<< "Angles and radii arrays in glyph data do not match.");
            return;
        }

        numberOfAngles = anglesArray->GetNumberOfTuples();
        /*for(int i = 0; i<numberOfAngles; i++)
        {
            // Get the two angles (azimuth and zenith)
            double * angles = anglesArray->GetTuple2(i);

            double r = radiiArray->GetComponent(0, i);

            printf("i:%d, anglesarray:[%f,%f], radiiarray:%f \n",i,angles[0],angles[1],r);
        }*/
	}

    //
    //      PREPARE OUTPUT POLYDATA
    //

    // Prepare scalars list
    QList<vtkDoubleArray*> outputScalarsList;
    int numberOfScalarTypes = inputPD->GetNumberOfArrays();
    for(int i = 0; i < numberOfScalarTypes; i++)
    {
        vtkDoubleArray* outputScalars = vtkDoubleArray::New();
        //char* name = inputPD->GetArray(i)->GetName();
        outputScalars->SetName(inputPD->GetArray(i)->GetName());
        outputScalarsList.append(outputScalars);
    }

    // Add new scalar list for SM
    /*vtkDoubleArray* SMScalars_external = vtkDoubleArray::New();
    SMScalars->SetName("SM_external");
    SMScalars->SetNumberOfComponents(1);

    vtkDoubleArray* SMScalars_internal = vtkDoubleArray::New();
    SMScalars->SetName("SM_internal");
    SMScalars->SetNumberOfComponents(1);*/

    vtkDoubleArray* SMScalars = vtkDoubleArray::New();
    SMScalars->SetName("SM");
    SMScalars->SetNumberOfComponents(1);

	// Create a point set for the output
	vtkPoints * outputPoints = vtkPoints::New();
	output->SetPoints(outputPoints);
	outputPoints->Delete();

	// Create a line array for the output
	vtkCellArray * outputLines = vtkCellArray::New();
	output->SetLines(outputLines);
	outputLines->Delete();

	//
    //      PROCESS
    //

	// Number of points in the current fiber, and a list of its point IDs
	vtkIdType numberOfPoints;
	vtkIdType * pointList;

    // Setup progress bar
    int numberOfCells = inputLines->GetNumberOfCells();
	int progressStep = numberOfCells / 25;
	progressStep += (progressStep == 0) ? 1 : 0;
	this->SetProgressText("Scoring fibers...");
	this->UpdateProgress(0.0);

	// Parameters
	double lambda;
	if(!ps->useGlyphData)
        lambda = 1.0;
    else
        lambda = ps->lambda;
    double beta = ps->beta;
    double muu = ps->muu;

    //printf("LAMBDA: %f\n",lambda);

	// Loop through all input fibers
	for (vtkIdType lineId = 0; lineId < numberOfCells; ++lineId)
	{
	    // Update the progress bar
		if ((lineId % progressStep) == 0)
		{
			this->UpdateProgress((double) lineId / (double) numberOfCells);
		}

        // Get the data of the current fiber
        vtkCell * currentCell = input->GetCell(lineId);
        int numberOfFiberPoints = currentCell->GetNumberOfPoints();
        double inv_numberOfFiberPoints = 1.0/(double)numberOfFiberPoints;

        // Get cell containing fiber
		inputLines->GetNextCell(numberOfPoints, pointList);

		// Create an ID list for the output fiber
		vtkIdList * newFiberList = vtkIdList::New();

        // Previous point coordinates
        double* prev_p = (double*) malloc(3*sizeof(double));
        double* prev2_p = (double*) malloc(3*sizeof(double));
        double* prev3_p = (double*) malloc(3*sizeof(double));

		// Current point coordinates
		double p[3];

		double PI = 3.14159265358;

        // Loop through all points in the fiber
		for (int pointId = 0; pointId < numberOfFiberPoints; ++pointId)
		{
            // Get the point ID of the current fiber point
			vtkIdType currentPointId = currentCell->GetPointId(pointId);

			// Copy the point coordinates to the output
			inputPoints->GetPoint(currentPointId, p);

			if(pointId > 3)
            {
                vtkIdType newPointId = outputPoints->InsertNextPoint(p);
                newFiberList->InsertNextId(newPointId);

                for(int i = 0; i < numberOfScalarTypes; i++)
                {
                    // Get the scalar value
                    double scalar = inputPD->GetArray(i)->GetTuple1(currentPointId);

                    // Copy the scalar value to the output
                    outputScalarsList.at(i)->InsertNextTuple1(scalar);
                }

                //
                // Compute score for current point
                //

                double external_energy = 0.0;
                if(ps->useGlyphData)
                {
                    // Find the corresponding voxel
                    vtkIdType imagePointId = this->inputVolume->FindPoint(p[0], p[1], p[2]);
                    if (imagePointId == -1)
                    {
                        // outside the data, return zero data dependent score
                        external_energy = 0.0;
                    }
                    else
                    {
                        // compute difference vector
                        double* dp = Difference(p,prev_p);

                        // transform to spherical coordinates
                        double theta = acos(dp[2]);
                        double phi = atan(dp[1]/dp[0]);

                        // find nearest vector in glyph data angles list
                        double cost = 1e30;
                        int matchedId;
                        for(int i = 0; i<numberOfAngles; i++)
                        {
                            double * angles = anglesArray->GetTuple2(i);
                            double ncost = ((theta+PI) - (angles[0]+PI))*((theta+PI) - (angles[0]+PI)) + ((phi+PI) - (angles[1]+PI))*((phi+PI) - (angles[1]+PI));
                            if(ncost < cost)
                            {
                                cost = ncost;
                                matchedId = i;
                            }
                        }

                        // External energy
                        external_energy = radiiArray->GetComponent(imagePointId, matchedId);

                        // Normalize glyph data
                        if(ps->normalizeGlyphData)
                            external_energy /= 30000.0;

                        // Apply log
                        if(ps->applyLog)
                        {
                            if(external_energy > 0.0)
                                external_energy = log10(external_energy);
                        }


                        //double* matchedAngles = anglesArray->GetTuple2(matchedId);
                        //printf("p:%f %f %f, prev_p:%f %f %f, dp:%f %f %f \n", p[0], p[1], p[2], prev_p[0], prev_p[1], prev_p[2], dp[0], dp[1], dp[2]);
                        //printf("pointId: %d, theta:%f, phi:%f \n", pointId, theta, phi);
                        //printf("matched angles: theta:%f, phi:%f radius:%f\n", matchedAngles[0], matchedAngles[1],radius);
                    }
                }

                // Internal energy
                double internal_energy = 0.0;
                if(lambda != 0.0)
                {
                    // compute frenet frame and its derivatives
                    double* a1 = HalvedDifference(p,prev_p);
                    double* t0 = Normalize(a1);
                    double* prev_t0 = Normalize(HalvedDifference(prev_p, prev2_p));
                    double* prev2_t0 = Normalize(HalvedDifference(prev2_p, prev3_p));
                    double* t1 = HalvedDifference(t0,prev_t0);
                    double* n0 = Normalize(t1);
                    double* prev_n0 = Normalize(HalvedDifference(prev_t0,prev2_t0));
                    double* n1 = HalvedDifference(n0,prev_n0);
                    double* b0 = Normalize(Cross(t0,n0));
                    double* b1 = Cross(t0,n1);

                    // compute scalars
                    double curvature = Norm(t1);
                    double torsion =  Norm(b1);

                    // release memory
                    free(a1);
                    free(t0);
                    free(prev_t0);
                    free(prev2_t0);
                    free(t1);
                    free(n0);
                    free(prev_n0);
                    free(n1);
                    free(b0);
                    free(b1);

                    internal_energy = sqrt(curvature*curvature + muu*torsion + beta*inv_numberOfFiberPoints*inv_numberOfFiberPoints);
                }

                //printf("curvature:%f, torsion:%f \n", curvature, torsion);

                // Total score
                double score = external_energy + lambda * internal_energy;

                SMScalars->InsertNextTuple1(score);
            }

            // Set previous points
            if(prev2_p != NULL)
                memcpy(prev3_p,prev2_p,3*sizeof(double));
            if(prev_p != NULL)
                memcpy(prev2_p,prev_p,3*sizeof(double));
            memcpy(prev_p,p,3*sizeof(double));
		}

		// Add the new fiber to the output
		outputLines->InsertNextCell(newFiberList);

        // release memory
        free(prev_p);
        free(prev2_p);
        free(prev3_p);
	}

	// Normalize SM
	int numTuples = SMScalars->GetNumberOfTuples();
	if(ps->standardizeScalars)
	{
	    double mean = 0.0;
	    int nrTuples = numTuples;
	    for(vtkIdType i = 0; i < nrTuples; i++)
        {
            mean += SMScalars->GetTuple1(i);
        }
	    mean /= nrTuples;

	    double std_dev = 0.0;
	    for(vtkIdType i = 0; i < nrTuples; i++)
        {
            double val = SMScalars->GetTuple1(i);
            std_dev += (val - mean)*(val - mean);
        }
        std_dev = sqrt(std_dev/nrTuples);

        for(vtkIdType i = 0; i < numTuples; i++)
        {
            SMScalars->SetTuple1(i,(SMScalars->GetTuple1(i) - mean)/std_dev );
        }
	}

	// Add scalar arrays
	for(int i = 0; i < numberOfScalarTypes; i++)
    {
        vtkDoubleArray* outputScalars = outputScalarsList.at(i);
        output->GetPointData()->AddArray(outputScalars);
    }

    // Add SM scalar array
    output->GetPointData()->AddArray(SMScalars);
    output->GetPointData()->SetActiveScalars("SM");

	// Finalize the progress bar
	this->UpdateProgress(1.0);
}

} // namespace bmia
