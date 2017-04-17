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

#ifndef bmia_PDBReaderPlugin_h
#define bmia_PDBReaderPlugin_h


/** Includes - Main Header */

#include "DTITool.h"

/** Includes - VTK */
#include <vtkPolyData.h>
#include <vtkCellArray.h>
#include <vtkPoints.h>
#include <vtkDoubleArray.h>
#include <vtkPointData.h>
#include <vtkTransform.h>
#include <vtkAlgorithm.h>

/** Includes - Qt */

#include <QFile>
#include <QByteArray>
#include <QMatrix4x4>

/** Includes - C++ */

//#include <iostream>
#include <fstream>
#include <assert.h>


namespace bmia {


/** A plugin for reading PDB fiber files.
*/

/** Data structs **/

struct StatHeader
{
    int i1;             /* is this statistic good as a luminance encoding for paths? (not used) */
    int i2;             /* is this statistic stored per point? */
    int i3;             /* does it show up in the stat panel by default? (not used) */
    char c1[255];       /* aggregate name */
    char c2[255];       /* local name */
    int i4;             /* unique ID */
};

struct AlgoHeader
{
    char c1[255];
    char c2[255];
    int i1;
};

struct Pathway
{
    int headerSize;
    int numPoints;
    int algoInt;
    int seedPointIndex;
    double* pathStats;
    double* points;         // {{x_1,y_1,z_1},{x_2,y_2,z_2},...,{x_n,y_n,z_n}} -- n = number of points
    double* pointStats;     // {{s1_1,s1_2,...s1_n},{s2_1,s2_2,...s2_n},...,{sm_1,sm_2,...sm_n}} -- m = number of point stats
};

/** Class **/

class PDBReaderPlugin : public plugin::Plugin, public data::Reader
{
    Q_OBJECT
    Q_INTERFACES(bmia::plugin::Plugin)
    Q_INTERFACES(bmia::data::Reader)

	public:

		/** Constructor */

		PDBReaderPlugin();

		/** Destructor */

		~PDBReaderPlugin();

		/** Returns the list of file extensions supported by this reader plugin. */

		QStringList getSupportedFileExtensions();

		/** Returns a list containing short descriptions of the supported file
			types. The number of descriptions and their order should match those
			of the list returned by "getSupportedFileExtensions". */

		QStringList getSupportedFileDescriptions();

		/** Load transfer function data from the given file and make it available
			to the data manager.
			@param filename Name if the desired transfer function file. */

		void loadDataFromFile(QString filename);

};


} // namespace bmia


#endif // bmia_PDBReaderPlugin_h
