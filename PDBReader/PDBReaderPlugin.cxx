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
 * PDBReaderPlugin.h
 *
 * 2013-05-16	Stephan Meesters
 * - First version
 *
 */

/** Includes */

#include "PDBReaderPlugin.h"



namespace bmia {


//-----------------------------[ Constructor ]-----------------------------\\

PDBReaderPlugin::PDBReaderPlugin() : plugin::Plugin("PDB Reader")
{

}


//------------------------------[ Destructor ]-----------------------------\\

PDBReaderPlugin::~PDBReaderPlugin()
{

}


//----------------------[ getSupportedFileExtensions ]---------------------\\

QStringList PDBReaderPlugin::getSupportedFileExtensions()
{
    QStringList list;
    list.push_back("pdb");
    return list;
}


//---------------------[ getSupportedFileDescriptions ]--------------------\\

QStringList PDBReaderPlugin::getSupportedFileDescriptions()
{
	QStringList list;
	list.push_back("PDB files");
	return list;
}


//---------------------------[ loadDataFromFile ]--------------------------\\

void PDBReaderPlugin::loadDataFromFile(QString filename)
{
	// Print status message to the log
	this->core()->out()->logMessage("Trying to load data from file " + filename);

	// Create the Qt file handler
	ifstream pdbFile (filename.toUtf8().constData(), ios::in | ios::binary );

	// Try to open the input file
	if (pdbFile.fail())
	{
		this->core()->out()->logMessage("Could not open file " + filename + "!");
		return;
	}

    // temp variables
    int ii; double d; char c; short s;

    // load header size
    unsigned int headersize;
    pdbFile.read(reinterpret_cast<char*>(&headersize), sizeof(unsigned int));

    // load transformation matrix
    QMatrix4x4 matrix;
    for(int i =0; i<4; i++)
    {
        for(int j = 0; j<4; j++)
        {
            pdbFile.read(reinterpret_cast<char*>(&d), sizeof(double));
            matrix(i,j) = d;
        }
    }

    // load number of stats
    int numstats;
    pdbFile.read(reinterpret_cast<char*>(&numstats), sizeof(int));
    std::cout << "Number of stats: " << numstats << std::endl;

    // load statheaders
    StatHeader* statheaders = new StatHeader[numstats];
    for(int i = 0; i<numstats; i++)
    {
        pdbFile.read(reinterpret_cast<char*>(&ii), sizeof(int));
        statheaders[i].i1 = ii;

        pdbFile.read(reinterpret_cast<char*>(&ii), sizeof(int));
        statheaders[i].i2 = ii;

        pdbFile.read(reinterpret_cast<char*>(&ii), sizeof(int));
        statheaders[i].i3 = ii;

        for(int j = 0; j<255; j++)
        {
            pdbFile.read(&c, sizeof(char));
            statheaders[i].c1[j] = c;
        }

        for(int j = 0; j<255; j++)
        {
            pdbFile.read(&c, sizeof(char));
            statheaders[i].c2[j] = c;
        }

        pdbFile.read(reinterpret_cast<char*>(&ii), sizeof(int));
        statheaders[i].i4 = ii;

        pdbFile.read(reinterpret_cast<char*>(&s), sizeof(short)); // dummy

        std::cout << "i1 " << statheaders[i].i1 << std::endl;
        std::cout << "i2 " << statheaders[i].i2 << std::endl;
        std::cout << "i3 " << statheaders[i].i3 << std::endl;
        std::cout << "c1 " << statheaders[i].c1 << std::endl;
        std::cout << "c2 " << statheaders[i].c2 << std::endl;
        std::cout << "i4 " << statheaders[i].i4 << std::endl;
    }

    // count number of point stats
    int numPointStats = 0;
    for(int i = 0; i<numstats; i++)
    {
        if(statheaders[i].i2 == 1)
            numPointStats++;
    }
    bool isScored = numPointStats > 0;
    if(!isScored)
        this->core()->out()->logMessage("PDB file contains no scoring values.");

    // number of algo's
    int numAlgos;
    pdbFile.read(reinterpret_cast<char*>(&numAlgos), sizeof(int));

    // load algoheaders
    AlgoHeader* algoheaders = new AlgoHeader[numAlgos];
    for(int i = 0; i<numAlgos; i++)
    {
        for(int j = 0; j<255; j++)
        {
            pdbFile.read(&c, sizeof(char));
            algoheaders[i].c1[j] = c;
        }

        for(int j = 0; j<255; j++)
        {
            pdbFile.read(&c, sizeof(char));
            algoheaders[i].c2[j] = c;
        }

        pdbFile.read(reinterpret_cast<char*>(&ii), sizeof(int));
        algoheaders[i].i1 = ii;

        pdbFile.read(reinterpret_cast<char*>(&s), sizeof(short)); // dummy
    }

    // version number
    int versionNumber;
    pdbFile.read(reinterpret_cast<char*>(&versionNumber), sizeof(int));
    std::cout << "PDB version number: " << versionNumber << std::endl;

    // skip to end of header
    pdbFile.seekg(headersize, pdbFile.beg);

    // number of pathways
    int numPathways;
    pdbFile.read(reinterpret_cast<char*>(&numPathways), sizeof(int));

	// create progress bar
	vtkAlgorithm* algo = vtkAlgorithm::New();
	algo->SetProgressText("Loading pathways ...");
	algo->UpdateProgress(0.0);
    this->core()->out()->createProgressBarForAlgorithm(algo, "PDB reader");

    // load pathways
    Pathway* pathways = new Pathway[numPathways];
    int totalNumberOfPoints = 0;
    for(int i = 0; i<numPathways; i++)
    {
        // update progress bar
        algo->UpdateProgress((double) i / (double) numPathways);

        pdbFile.read(reinterpret_cast<char*>(&ii), sizeof(int));
        pathways[i].headerSize = ii;

        int pos = pdbFile.tellg();

        pdbFile.read(reinterpret_cast<char*>(&ii), sizeof(int));
        pathways[i].numPoints = ii;

        totalNumberOfPoints += pathways[i].numPoints;

        pdbFile.read(reinterpret_cast<char*>(&ii), sizeof(int));
        pathways[i].algoInt = ii;

        pdbFile.read(reinterpret_cast<char*>(&ii), sizeof(int));
        pathways[i].seedPointIndex = ii;

        pathways[i].pathStats = new double[numstats];
        for(int j = 0; j<numstats; j++)
        {
            pdbFile.read(reinterpret_cast<char*>(&d), sizeof(double));
            pathways[i].pathStats[j] = d;
        }

        // skip to end of header of pathway
        pdbFile.seekg(pathways[i].headerSize + pos, pdbFile.beg);

        pathways[i].points = new double[pathways[i].numPoints*3]; // 3D points
        for(int j = 0; j<pathways[i].numPoints; j++)
        {
            for(int k =0; k<3; k++)
            {
                pdbFile.read(reinterpret_cast<char*>(&d), sizeof(double));
                pathways[i].points[j*3+k] = d;
            }
        }

        if(isScored)
        {
            pathways[i].pointStats = new double[pathways[i].numPoints * numstats]; // scalar score per point

            for(int k =0; k<numstats; k++)
            {
                // stat value per point
                if(statheaders[k].i2 == 1)
                {
                    for(int j = 0; j<pathways[i].numPoints; j++)
                    {
                        pdbFile.read(reinterpret_cast<char*>(&d), sizeof(double));
                        pathways[i].pointStats[k*pathways[i].numPoints + j] = d;
                    }
                }

                // stat value per fiber (single value)
                else
                {
                    for(int j = 0; j<pathways[i].numPoints; j++)
                    {
                        pathways[i].pointStats[k*pathways[i].numPoints + j] = pathways[i].pathStats[k];
                    }
                }
            }
        }
    }

    // close .pdb file
    pdbFile.close();

    // notify load complete
    this->core()->out()->logMessage("PDB processed: " + filename);

    //
    //  Transformation
    //
    for(int i = 0; i<numPathways; i++)
    {
        for(int j = 0; j<pathways[i].numPoints; j++)
        {
            // translation
            pathways[i].points[j*3+0] += 80;
            pathways[i].points[j*3+1] += 120;
            pathways[i].points[j*3+2] += 60;

            // scaling
            pathways[i].points[j*3+0] *= 0.5;
            pathways[i].points[j*3+1] *= 0.5;
            pathways[i].points[j*3+2] *= 0.5;
        }
    }

    vtkMatrix4x4 * mat = vtkMatrix4x4::New();
    mat->Identity();

    vtkTransform* transform = vtkTransform::New();
    transform->SetMatrix(mat);
    transform->Translate(-80,-120,-60);
    transform->Scale(2,2,2);
    mat = transform->GetMatrix();

    //
    //  VTK conversion
    //

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

    // Array holding the ConTrack scoring values
    QList<vtkDoubleArray*> scoringList;
    if(isScored)
    {
        for(int k = 0; k<numstats; k++)
        {
            vtkDoubleArray* scoring = vtkDoubleArray::New();
            scoring->SetName(statheaders[k].c1);
            scoring->SetNumberOfTuples(totalNumberOfPoints);
            scoringList.append(scoring);
        }
    }

    // Loop over pathways
    int counter = 0;
    for(int i = 0; i < numPathways; i++)
    {
        int numberOfFiberPoints = pathways[i].numPoints;

        // Create a cell representing a fiber
        outputLines->InsertNextCell(numberOfFiberPoints);

        // Loop over points in the pathway
        for(int j = 0; j<numberOfFiberPoints; j++)
        {
            outputPoints->InsertNextPoint(pathways[i].points[j*3],pathways[i].points[j*3+1],pathways[i].points[j*3+2]);
            outputLines->InsertCellPoint(counter + j);
            if(isScored)
            {
                for(int k = 0; k<numstats; k++)
                {
                    vtkDoubleArray* scoring = scoringList.at(k);
                    scoring->SetTuple1(counter + j, pathways[i].pointStats[k*numberOfFiberPoints + j]);
                }
            }
        }

        counter += numberOfFiberPoints;
    }

    // Set active scalars to ConTrack score
    if(isScored)
    {
        for(int k = 0; k<numstats; k++)
        {
            vtkDoubleArray* scoring = scoringList.at(k);
            output->GetPointData()->AddArray(scoring);
        }
        if(numstats >= 2)
            output->GetPointData()->SetActiveScalars(statheaders[2].c1); // Default: "Scoring"
        else
            output->GetPointData()->SetActiveScalars(statheaders[0].c1);
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

	// remove progress bar
    this->core()->out()->deleteProgressBarForAlgorithm(algo);

    // clean up
    algo->Delete();
    transform->Delete();
}

} // namespace bmia


Q_EXPORT_PLUGIN2(libPDBReaderPlugin, bmia::PDBReaderPlugin)
