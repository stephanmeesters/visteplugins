
/** Includes */

#include "vtkFiberROICutting.h"


namespace bmia {


vtkStandardNewMacro(vtkFiberROICutting);


//-----------------------------[ Constructor ]-----------------------------\\

vtkFiberROICutting::vtkFiberROICutting()
{
	// Set default options
    this->cutAtROIEnds = true;
}


//------------------------------[ Destructor ]-----------------------------\\

vtkFiberROICutting::~vtkFiberROICutting()
{

}


//-------------------------------[ Execute ]-------------------------------\\

void vtkFiberROICutting::Execute()
{
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

    // Get list of scalars
    QList<vtkDoubleArray*> outputScalarsList;
    int numberOfScalarTypes = inputPD->GetNumberOfArrays();
    for(int i = 0; i < numberOfScalarTypes; i++)
    {
        vtkDoubleArray* outputScalars = vtkDoubleArray::New();
        outputScalars->SetName(inputPD->GetArray(i)->GetName());
        outputScalarsList.append(outputScalars);
    }

	// Create a point set for the output
	vtkPoints * outputPoints = vtkPoints::New();
	output->SetPoints(outputPoints);
	outputPoints->Delete();

	// Create a line array for the output
	vtkCellArray * outputLines = vtkCellArray::New();
	output->SetLines(outputLines);
	outputLines->Delete();

	// Get image data from ROI data
	vtkImageData * roi = this->roiData->getVtkImageData();
    if (!roi)
	{
		vtkErrorMacro(<< "ROI data has not been set.");
		return;
	}

	// Convert ROI datatype to char
    vtkImageCast* imageCast = vtkImageCast::New();
    imageCast->SetInput(roi);
    imageCast->SetOutputScalarTypeToChar();
    imageCast->Update();
    vtkImageData * roiChar = imageCast->GetOutput();

    // Get the transformation matrix
    vtkObject* tfm;
    vtkMatrix4x4* transformationMatrix;
    if (roiData->getAttributes()->getAttribute("transformation matrix", tfm ))
    {
        transformationMatrix = vtkMatrix4x4::SafeDownCast(tfm);
        if (transformationMatrix == 0)
        {
            vtkErrorMacro(<< "not a valid transformation matrix");
            return;
        }
    }
    vtkMatrix4x4* transformationMatrixInv = vtkMatrix4x4::New();
    vtkMatrix4x4::Invert(transformationMatrix,transformationMatrixInv);

	// Number of points in the current fiber, and a list of its point IDs
	vtkIdType numberOfPoints;
	vtkIdType * pointList;

    // Setup progress bar
    int numberOfCells = inputLines->GetNumberOfCells();
	int progressStep = numberOfCells / 25;
	progressStep += (progressStep == 0) ? 1 : 0;
	this->SetProgressText("Cutting fibers using ROI...");
	this->UpdateProgress(0.0);

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

        // Get cell containing fiber
		inputLines->GetNextCell(numberOfPoints, pointList);

		// Evaluate if the fiber should be included in the output fibers
        double p[3];
        char prev_roi;
        bool excludeFiber = true;
        int beginPoint = -1, endPoint = -1, cuttedFiberLength = -1e30;
        for (int pointId = 0; pointId < numberOfFiberPoints; ++pointId)
		{
		    // Get the point ID of the current fiber point
			vtkIdType currentPointId = currentCell->GetPointId(pointId);

			// Get point coordinates
			inputPoints->GetPoint(currentPointId, p);

		    // Check if point is within ROI
			double vec[4] = {p[0],p[1],p[2],1};
			fiberMatrix->MultiplyPoint(vec,vec);
			transformationMatrixInv->MultiplyPoint(vec,vec);

            //printf("x:%d y:%d z:%d \n",(int)vec[0],(int)vec[1],(int)vec[2]);

            char* roiPixel = static_cast<char*>(roiChar->GetScalarPointer((int)vec[0],(int)vec[1],(int)vec[2]));

            if(pointId > 0)
            {
                if(this->cutAtROIEnds)
                {
                    // look for entrance or exit into ROI
                    char diff = roiPixel[0] - prev_roi;
                    //printf("pix1:%d pix2:%d diff:%d \n",roiPixel[0], prev_roi, diff);
                    if(diff != 0)
                    {
                        if(beginPoint == -1)
                            beginPoint = pointId;
                        else// if(pointId - beginPoint > cuttedFiberLength)
                        {
                            endPoint = pointId;
                            cuttedFiberLength = endPoint - beginPoint;
                        }
                    }
                }

                // include fiber if it enters ROI
                else if(roiPixel[0] == 1)
                {
                    excludeFiber = false;
                    break;
                }
            }

            // store previous ROI pixel
            prev_roi = roiPixel[0];
		}

		if(this->cutAtROIEnds)
		{
		    printf("beginpoint:%d endpoint:%d \n",beginPoint,endPoint);
            //if(!(beginPoint == -1 || endPoint == -1))
                excludeFiber = false;
		}
		else
		{
		    beginPoint = 0;
		    endPoint = numberOfFiberPoints;
		}

        if(!excludeFiber)
        {
            // Create an ID list for the output fiber
            vtkIdList * newFiberList = vtkIdList::New();

            // Current scalar value
            double scalar;

            // Loop through all points in the fiber
            for (int pointId = 0; pointId < numberOfFiberPoints; ++pointId)
            {
                // Get the point ID of the current fiber point
                vtkIdType currentPointId = currentCell->GetPointId(pointId);

                // Get point coordinates
                inputPoints->GetPoint(currentPointId, p);

                // Check if point is within ROI
                if(this->cutAtROIEnds)
                {
                    double vec[4] = {p[0],p[1],p[2],1};
                    fiberMatrix->MultiplyPoint(vec,vec);
                    transformationMatrixInv->MultiplyPoint(vec,vec);
                    char* roiPixel = static_cast<char*>(roiChar->GetScalarPointer((int)vec[0],(int)vec[1],(int)vec[2]));
                    if(roiPixel[0] == 0)
                        continue;
                }



                // Copy coordinates to output
                vtkIdType newPointId = outputPoints->InsertNextPoint(p);
                newFiberList->InsertNextId(newPointId);

                for(int i = 0; i < numberOfScalarTypes; i++)
                {
                    // Get the scalar value
                    double scalar = inputPD->GetArray(i)->GetTuple1(currentPointId);

                    // Copy the scalar value to the output
                    outputScalarsList.at(i)->InsertNextTuple1(scalar);
                }
            }

            // Add the new fiber to the output
            outputLines->InsertNextCell(newFiberList);
        }
	}

	// Add scalar arrays
	for(int i = 0; i < numberOfScalarTypes; i++)
    {
        vtkDoubleArray* outputScalars = outputScalarsList.at(i);
        output->GetPointData()->AddArray(outputScalars);
    }

	// Finalize the progress bar
	this->UpdateProgress(1.0);
}

} // namespace bmia
