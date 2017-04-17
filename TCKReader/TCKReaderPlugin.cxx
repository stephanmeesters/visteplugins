/**
 * Copyright (c) 2012, Biomedical Image Analysis Eindhoven (BMIA/e)
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 *
 *   - Redistributions of source code must retain the above copyright
 *     notice, this list of conditions and the following disclaimer.
 *
 *   - Redistributions in binary form must reproduce the above copyright
 *     notice, this list of conditions and the following disclaimer in
 *     the documentation and/or other materials provided with the
 *     distribution.
 *
 *   - Neither the name of Eindhoven University of Technology nor the
 *     names of its contributors may be used to endorse or promote
 *     products derived from this software without specific prior
 *     written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
 * FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
 * COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
 * INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
 * BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
 * CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
 * LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
 * ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 */

/*
 * TCKReaderPlugin.h
 *
 * 2013-05-16	Stephan Meesters
 * - First version
 *
 */

/** Includes */

#include "TCKReaderPlugin.h"



namespace bmia {


//-----------------------------[ Constructor ]-----------------------------\\

TCKReaderPlugin::TCKReaderPlugin() : plugin::Plugin("TCK Reader")
{

}


//------------------------------[ Destructor ]-----------------------------\\

TCKReaderPlugin::~TCKReaderPlugin()
{

}


//----------------------[ getSupportedFileExtensions ]---------------------\\

QStringList TCKReaderPlugin::getSupportedFileExtensions()
{
    QStringList list;
    list.push_back("tck");
    return list;
}


//---------------------[ getSupportedFileDescriptions ]--------------------\\

QStringList TCKReaderPlugin::getSupportedFileDescriptions()
{
	QStringList list;
	list.push_back("TCK files");
	return list;
}


//---------------------------[ loadDataFromFile ]--------------------------\\

void TCKReaderPlugin::loadDataFromFile(QString filename)
{
	// Print status message to the log
	this->core()->out()->logMessage("Trying to load data from file " + filename);

	// Create the Qt file handler
	ifstream TCKFile (filename.toUtf8().constData(), ios::in | ios::binary );

	// Try to open the input file
	if (TCKFile.fail())
	{
		this->core()->out()->logMessage("Could not open file " + filename + "!");
		return;
	}

	// temp variables
    int ii; double d; char c; short s;

    char endBuffer[] = "AAA";
    char end[] = "END";
    for(int i = 0; i<1000; i++)
    {
        TCKFile.read(reinterpret_cast<char*>(&c), sizeof(char));


        endBuffer[0] = endBuffer[1];
        endBuffer[1] = endBuffer[2];
        endBuffer[2] = c;

        if(!strcmp(endBuffer,end))
            break;
    }

	double inf=std::numeric_limits<double>::infinity();
    //double nan = std::numeric_limits<double>::
    QList< QList<float> > fibersList;

    int headerPos = TCKFile.tellg();
    int bugfixer = 0;
    int MAXTRIES = 200;
    while(bugfixer < MAXTRIES)
    {
        TCKFile.seekg(headerPos, TCKFile.beg);

        for(int i = 0; i<bugfixer; i++)
        {
            TCKFile.read(reinterpret_cast<char*>(&c), sizeof(char));
        }

        int i = 0, j = 0;
        fibersList.clear();

        int pos;
        float f;
        while(true)
        {
            pos = TCKFile.tellg();
            TCKFile.read(reinterpret_cast<char*>(&f), sizeof(float));
            if(abs(f) > 0.01)
                break;
        }
        TCKFile.seekg(pos, TCKFile.beg);

        int k = 0;
        bool abort;
        while(true)
        {
            //algo->UpdateProgress((double) j);

            QList<float> fiber;

            float f1,f2,f3;
            k = 0;
            abort = false;
            while(true)
            {
                TCKFile.read(reinterpret_cast<char*>(&f1), sizeof(float));
                TCKFile.read(reinterpret_cast<char*>(&f2), sizeof(float));
                TCKFile.read(reinterpret_cast<char*>(&f3), sizeof(float));

                if(f1 != f1 && f2!=f2 && f3!=f3)
                    break;

                if(f1 == inf && f2==inf && f3==inf)
                    break;

                if(TCKFile.eof())
                    break;

                f1 += 80;
                f1 *= 0.5;

                f2 += 120;
                f2 *= 0.5;

                f3 += 60;
                f3 *= 0.5;

                fiber.append(f1);
                fiber.append(f2);
                fiber.append(f3);

                k++;

                if(f1 > 10000 || f2 > 10000 || f3 > 10000)
                {
                    abort = true;
                    break;
                }
            }

            if(abort)
            {
                break;
            }

            fibersList.append(fiber);

            j++;

            if(TCKFile.eof())
                break;

            if(f1 == inf && f2==inf && f3==inf)
                break;
        }

        if(abort)
        {
            printf("Bugfix attempt: %d\n",bugfixer);
            bugfixer++;
        }
        else
            break;

    }

    if(bugfixer == MAXTRIES)
    {
        this->core()->out()->logMessage("Error loading " + filename + "!");
        return;
    }

    //
    //  Transformation
    //

    vtkMatrix4x4 * mat = vtkMatrix4x4::New();
    mat->Identity();

    vtkTransform* transform = vtkTransform::New();
    transform->SetMatrix(mat);
    transform->Translate(-80,-120,-60);
    transform->Scale(2,2,2);
    mat = transform->GetMatrix();


    // close .TCK file
    TCKFile.close();

     // Create polydata
    vtkPolyData* output = vtkPolyData::New();

	// Create a point set for the output
	vtkPoints * outputPoints = vtkPoints::New();
	output->SetPoints(outputPoints);
	outputPoints->Delete();

    // Cell array holding the pathways.
    // Each pathway is a single cell (line).
    // Each cell holds the id values to points from the vtkPoints list
    vtkCellArray * outputLines = vtkCellArray::New();
    output->SetLines(outputLines);
    outputLines->Delete();

    // Loop over pathways
    int counter = 0;
    int numPathways = fibersList.length();
    for(int i = 0; i < numPathways; i++)
    {
        QList<float> fiber = fibersList.at(i);

        int numberOfFiberPoints = fiber.length()/3;

        // Create a cell representing a fiber
        outputLines->InsertNextCell(numberOfFiberPoints);

        // Loop over points in the pathway
        for(int j = 0; j<numberOfFiberPoints; j++)
        {
            outputPoints->InsertNextPoint(fiber[j*3],fiber[j*3+1],fiber[j*3+2]);
            outputLines->InsertCellPoint(counter + j);
        }

        counter += numberOfFiberPoints;
    }

    //
    //  Save to dataset
    //

    // Short name of the data set
	QString shortName = filename;

	// Find the last slash
	int lastSlash = filename.lastIndexOf("/");

	// If the filename does not contain a slash, try to find a backslash
	if (lastSlash == -1)
	{
		lastSlash = filename.lastIndexOf("\\");
	}

	// Throw away everything up to and including the last slash
	if (lastSlash != -1)
	{
		shortName = shortName.right(shortName.length() - lastSlash - 1);
	}

	// Find the last dot in the remainder of the filename
	int lastPoint = shortName.lastIndexOf(".");

	// Throw away everything after and including the last dot
	if (lastPoint != -1)
	{
		shortName = shortName.left(lastPoint);
	}

	// Create a new data set for the transfer function
	data::DataSet * ds = new data::DataSet(shortName, "fibers", output);

	// Fibers should be visible, and the visualization pipeline should be updated
	ds->getAttributes()->addAttribute("isVisible", 1.0);
	ds->getAttributes()->addAttribute("updatePipeline", 1.0);

	// Copy the transformation matrix to the output
	ds->getAttributes()->addAttribute("transformation matrix", vtkObject::SafeDownCast(mat));

    // Add the data set to the manager
	this->core()->data()->addDataSet(ds);
}

} // namespace bmia


Q_EXPORT_PLUGIN2(libTCKReaderPlugin, bmia::TCKReaderPlugin)
