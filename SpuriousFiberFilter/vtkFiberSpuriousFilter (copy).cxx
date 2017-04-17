
/** Includes */

#include "vtkFiberSpuriousFilter.h"
#include <ctime>


namespace bmia {


vtkStandardNewMacro(vtkFiberSpuriousFilter);


//-----------------------------[ Constructor ]-----------------------------\\

vtkFiberSpuriousFilter::vtkFiberSpuriousFilter()
{
	// Set default options
    this->inputVolume = NULL;
}


//------------------------------[ Destructor ]-----------------------------\\

vtkFiberSpuriousFilter::~vtkFiberSpuriousFilter()
{

}

//-------------------------------[ SetInputVolume ]-------------------------------\\

void vtkFiberSpuriousFilter::SetInputVolume(vtkImageData * image)
{
	// Store the pointer
	this->inputVolume = image;
}

//-------------------------------[ SetParameters ]-------------------------------\\

void vtkFiberSpuriousFilter::SetParameters(ParameterSettings* ps)
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
    double output = vec[0]*vec[0] + vec[1]*vec[1] + vec[2]*vec[2];
    free(vec);
    return output;
}

void PrintVector3(double* vec)
{
    if(vec == NULL)
        return;
    printf("%f, %f, %f \n", vec[0], vec[1], vec[2]);
}

void PrintVector2(double* vec)
{
    if(vec == NULL)
        return;
    printf("{%f, %f} \n", vec[0], vec[1]);
}

void PrintMatrix3x3(double* mat)
{
    printf("{%f, %f, %f \n", mat[0], mat[1], mat[2]);
    printf(" %f, %f, %f \n", mat[3], mat[4], mat[5]);
    printf(" %f, %f, %f} \n", mat[6], mat[7], mat[8]);
}


double* EulerAngles(double* input)
{
    double x = input[0];
    double y = input[1];
    double z = input[2];
    double* output = (double*)malloc(2*sizeof(double));
    if((x - 1)*(x - 1) < 1e-10)
    {
        output[0] = 1.5708;
        output[1] = 0;
    }
    else if((x + 1)*(x + 1) < 1e-10)
    {
        output[0] = -1.5708;
        output[1] = 0;
    }
    else
    {
        output[0] = acos(sqrt(y*y + z*z)*(z>=0?1:-1))*(x>0?1:-1);
        output[1] = -asin(y/sqrt(y*y + z*z)*(z>=0?1:-1));
    }
    return output;
}

double* R(double* input)
{
    double beta = input[0];
    double gamma = input[1];
    double* output = (double*)malloc(9*sizeof(double));
    output[0] = cos(beta);
    output[1] = 0;
    output[2] = sin(beta);
    output[3] = sin(beta)*sin(gamma);
    output[4] = cos(gamma);
    output[5] = -cos(beta)*sin(gamma);
    output[6] = -cos(gamma)*sin(beta);
    output[7] = sin(gamma);
    output[8] = cos(beta)*cos(gamma);

    free(input);

    return output;
}

double* Transpose3x3(double *src)
{
    double* dst = (double*)malloc(sizeof(double)*9);
    #pragma omp parallel
    for(int n = 0; n<9; n++)
    {
        int i = n/3;
        int j = n%3;
        dst[n] = src[3*j + i];
    }

    free(src);

    return dst;
}

double* Multiply(double* mat, double* vec)
{
    double* dst = (double*)malloc(sizeof(double)*3);
    dst[0] = mat[0]*vec[0] + mat[1]*vec[1]+mat[2]*vec[2];
    dst[1] = mat[3]*vec[0] + mat[4]*vec[1]+mat[5]*vec[2];
    dst[2] = mat[6]*vec[0] + mat[7]*vec[1]+mat[8]*vec[2];

    free(mat);

    return dst;
}

double* Subtract3(double* vec, double* vec2)
{
    double* dst = (double*)malloc(sizeof(double)*3);
    dst[0] = vec[0]-vec2[0];
    dst[1] = vec[1]-vec2[1];
    dst[2] = vec[2]-vec2[2];
    return dst;
}

double cot(double i) { return(1 / tan(i)); }

double kernel(double x, double y, double z, double b, double g, double D33, double D44, double t)
{
    return
    10.026513098524001*sqrt(D33*D44)*t*sqrt(D33*t)
    *
    (abs(b) < PI
    ?
    1.0/(100.53096491487338*D33*D44*t*t)*
        exp((-0.25*sqrt((1.*pow(-0.25*b*z +
                   (x*cos(0.5*b))/(1. - 0.041666666666666664*b*b),2))/(D33*D44)\
               + pow(b*b/D44 +
                pow(0.5*b*x + (0.5*z*cos(0.5*b))/
                    (1. - 0.041666666666666664*pow(b,2)),2)/D33,2)))/t)
    :
    1./(100.53096491487338*D33*D44*t*t)*
        exp((-0.25*sqrt((1.*pow(-0.25*b*z + 0.5*b*x*cot(0.5*b),2))/(D33*D44) +
              pow((b*b)/D44 + pow(0.5*b*x + 0.25*b*z*cot(0.5*b),2)/D33,2)))/t)
        )
    *
    (abs(g) < 0.3141592653589793
    ?
       1./(100.53096491487338*D33*D44*t*t)*
        exp((-0.25*sqrt((1.*pow(-0.25*g*z -
                   (1.*y*cos(0.5*g))/(1. - 0.041666666666666664*pow(g,2)),2))/
               (D33*D44) + pow((g*g)/D44 +
                pow(-0.5*g*y + (0.5*z*cos(0.5*g))/
                    (1. - 0.041666666666666664*pow(g,2)),2)/D33,2)))/t)
    :
       1./(100.53096491487338*D33*D44*pow(t,2))*
        exp((-0.25*sqrt((1.*pow(-0.25*g*z - 0.5*g*y*cot(0.5*g),2))/(D33*D44) +
              pow(pow(g,2)/D44 + pow(-0.5*g*y + 0.25*g*z*cot(0.5*g),2)/D33,2)))/t
          )
    );

}



double vtkFiberSpuriousFilter::k2(double* x, double* y, double* r, double* v)
{
    //printf("x:{%f,%f,%f}, y:{%f,%f,%f}, r:{%f,%f,%f}, v:{%f,%f,%f} \n", x[0], x[1], x[2], y[0], y[1], y[2], r[0], r[1], r[2], v[0], v[1], v[2]);

    //double* arg11 = R(EulerAngles(v));
    //PrintMatrix3x3(arg11);

    double* arg1 = Multiply(Transpose3x3(R(EulerAngles(v))),Subtract3(x,y));
    //PrintVector3(arg1);
//    PrintVector3(Subtract3(x,y));
//    PrintMatrix3x3(Transpose3x3(R(EulerAngles(v))));
//    PrintMatrix3x3(R(EulerAngles(v)));
//    PrintVector2(EulerAngles(v));

    //double* arg2 = EulerAngles(Multiply(Transpose3x3(R(EulerAngles(v))),r));
    double* arg2p = Multiply(Transpose3x3(R(EulerAngles(v))),r);
    double* arg2 = EulerAngles(arg2p);


//    PrintMatrix3x3(Transpose3x3(R(EulerAngles(v))));
//    PrintVector3(r);
//    PrintVector3(Multiply(Transpose3x3(R(EulerAngles(v))),r));
//    PrintVector2(arg2);

    double kernelval = kernel(arg1[0],arg1[1],arg1[2],arg2[0],arg2[1], this->ps->D33, this->ps->D44, this->ps->t);

    //free(arg1);
    //free(arg2);
    //free(arg2p);

    //printf("KERNEL VALUE: %f \n",kernelval);
    return kernelval;


    //double aap[3];
   //Transpose3x3(arg1,aap,3,3);


    //PrintMatrix3x3(aap);

//    double* arg1 = MatrixMult(Transpose(R(EulerAngles(v))),Subtract(x,y));
//    double* arg2 = EulerAngles(MatrixMult(Transpose(R(EulerAngles(v)))));
//    return kernel(arg1[0],arg1[1],arg1[2],arg2[0],arg2[1]);
}


void vtkFiberSpuriousFilter::Execute()
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

	vtkPolyData* inputCopy =  vtkPolyData::New();
	inputCopy->DeepCopy(input);

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
    //      PREPARE OUTPUT POLYDATA
    //

    // Prepare scalars list
    QList<vtkDoubleArray*> outputScalarsList;
    int numberOfScalarTypes = inputPD->GetNumberOfArrays();
    for(int i = 0; i < numberOfScalarTypes; i++)
    {
        vtkDoubleArray* outputScalars = vtkDoubleArray::New();
        outputScalars->SetName(inputPD->GetArray(i)->GetName());
        outputScalarsList.append(outputScalars);
    }

    // Add new scalar list for SM
    vtkDoubleArray* SMScalars = vtkDoubleArray::New();
    SMScalars->SetName("SpuriousFilter");
    SMScalars->SetNumberOfComponents(1);

    // Create temporary scalar list for tangents
    vtkDoubleArray* tangentScalars = vtkDoubleArray::New();
    tangentScalars->SetNumberOfComponents(3);

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
	vtkIdType numberOfPointsf2;
	vtkIdType * pointListf2;

    // Setup progress bar
    int numberOfCells = inputLines->GetNumberOfCells();
	int progressStep = 2;//numberOfCells / 25;
	progressStep += (progressStep == 0) ? 1 : 0;
	this->SetProgressText("Filtering fibers...");
	this->UpdateProgress(0.0);

	double minDistSquared = ps->minDist * ps->minDist;


	// save fibers in raw form


    // Compute all fiber tangents
	for (vtkIdType lineId = 0; lineId < numberOfCells; ++lineId)
	{
        // Get the data of the current fiber
        vtkCell * currentCell = input->GetCell(lineId);

        // Get cell containing fiber
		inputLines->GetNextCell(numberOfPoints, pointList);
		int numberOfFiberPoints = currentCell->GetNumberOfPoints();

        // Previous point coordinates
        double* prev_p = (double*) malloc(3*sizeof(double));

		// Current point coordinates
		double p[3];

        // Loop through all points in the fiber
		for (int pointId = 0; pointId < numberOfFiberPoints; ++pointId)
		{
            // Get the point ID of the current fiber point
			vtkIdType currentPointId = currentCell->GetPointId(pointId);

			// Copy the point coordinates to the output
			inputPoints->GetPoint(currentPointId, p);

			if(pointId > 1)
            {
                double* tangent = Normalize(Difference(p,prev_p));
                tangentScalars->InsertNextTuple3(tangent[0],tangent[1],tangent[2]);
                //printf("computed tangent: {%f,%f} \n", tangent[0],tangent[1]);
            }
            else
            {
                tangentScalars->InsertNextTuple3(0.0, 0.0, 0.0);
            }

            // Set previous point
            memcpy(prev_p,p,3*sizeof(double));
		}

        // release memory
        free(prev_p);
	}

	clock_t lastTime = clock();

	// Compute kernel for fibers
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
		//inputLines->GetCell(lineId,numberOfPoints, pointList);

		//printf("asdasdasdasd %d %d\n",numberOfFiberPoints,numberOfPoints);

		// Create an ID list for the output fiber
		vtkIdList * newFiberList = vtkIdList::New();

		// Current point coordinates
		double p[3];

        // Loop through all points in the fiber
		for (int pointId = 0; pointId < numberOfFiberPoints; ++pointId)
		{

            // Get the point ID of the current fiber point
			vtkIdType currentPointId = currentCell->GetPointId(pointId);

			// Copy the point coordinates to the output
			inputPoints->GetPoint(currentPointId, p);

            vtkIdType newPointId = outputPoints->InsertNextPoint(p);
            newFiberList->InsertNextId(newPointId);

            if(pointId > 1)
            {

                // copy the other scalars
//                for(int i = 0; i < numberOfScalarTypes; i++)
//                {
//                    // Get the scalar value
//                    double scalar = inputPD->GetArray(i)->GetTuple1(currentPointId);
//
//                    // Copy the scalar value to the output
//                    outputScalarsList.at(i)->InsertNextTuple1(scalar);
//                }

                //
                // Compute score for current point
                //

                // load f1 tangent
                double* tangent = tangentScalars->GetTuple3(currentPointId);
                //printf("loaded tangent: {%f,%f} \n", tangent[0],tangent[1]);

                // second loop around f2
                double score = 0;
                int totalNumPoint = 0;
                double pf2[3];
                for (vtkIdType lineIdf2 = 0; lineIdf2 < numberOfCells; ++lineIdf2)
                {
                    // skip equal fibers to prevent unnecessary calculations
                    if(lineId == lineIdf2)
                        continue;
//
                    vtkCell * currentCellf2 = inputCopy->GetCell(lineIdf2);
                    int numberOfFiberPointsf2 = currentCellf2->GetNumberOfPoints();

                    //inputLines->GetCell(lineIdf2,numberOfPointsf2, pointListf2);

                    for (int pointIdf2 = 0; pointIdf2 < numberOfFiberPointsf2; ++pointIdf2)
                    {
                        vtkIdType currentPointIdf2 = currentCellf2->GetPointId(pointIdf2);
                        inputPoints->GetPoint(currentPointIdf2, pf2);
                        if(pointIdf2 > 1)
                        {
                            double* tangentf2 = tangentScalars->GetTuple3(currentPointIdf2);

                            // check distance
                            double sval;
                            if(Norm(Subtract3(p,pf2)) + 10.0*Norm(Subtract3(tangent, tangentf2)) < minDistSquared)
                            {
                                // evaluate kernel
                                sval = k2(p,pf2,tangent,tangentf2);
                                //if(pointIdf2 == 10)
                                //    return;
                            }
                            else
                            {
                                // out of range, so output zero
                                sval = 0.0;
                            }

                            // add scoring value. check against NaN
                            //if(!(sval!=sval))
                                score += sval;

                            tangentf2 = NULL;
                        }
                    }
                }
                SMScalars->InsertNextTuple1(score);
                //printf("score for fiber point: %f \n", score);
            }
            else
                SMScalars->InsertNextTuple1(0.0);
		}

		// Add the new fiber to the output
		outputLines->InsertNextCell(newFiberList);

		// comptue remaining time
		clock_t newTime = clock();
		QString progtext = QString("Filtering fibers... (time remaining = %1 seconds)").arg(int(double(newTime - lastTime) / CLOCKS_PER_SEC * (numberOfCells - lineId)));
		this->SetProgressText(progtext.toLocal8Bit().data());
		lastTime = newTime;
	}

	// Add scalar arrays
	for(int i = 0; i < numberOfScalarTypes; i++)
    {
        vtkDoubleArray* outputScalars = outputScalarsList.at(i);
        output->GetPointData()->AddArray(outputScalars);
    }

    // Add SM scalar array
    output->GetPointData()->AddArray(SMScalars);
    output->GetPointData()->SetActiveScalars("SpuriousFilter");

	// Finalize the progress bar
	this->UpdateProgress(1.0);

	inputCopy->Delete();
}

} // namespace bmia
