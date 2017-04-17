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

#ifndef vtkInteractorStyleTrackballPositionPicker_h
#define vtkInteractorStyleTrackballPositionPicker_h

/** Includes - Standard C++ */
#include <iostream>

/** Includes - VTK */
#include <vtkSmartPointer.h>
#include <vtkCellArray.h>
#include <vtkCellData.h>
#include <vtkDoubleArray.h>
#include <vtkPoints.h>
#include <vtkPolyLine.h>
#include <vtkPolyData.h>
#include <vtkPolyDataMapper.h>
#include <vtkActor.h>
#include <vtkRenderWindow.h>
#include <vtkRenderer.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkRendererCollection.h>

#include <vtkPointPicker.h>
#include <vtkCellPicker.h>
#include <vtkWorldPointPicker.h>
#include <vtkPropPicker.h>
#include <vtkPropCollection.h>
#include <vtkInteractorStyleTrackballCamera.h>
#include <vtkObjectFactory.h>
#include <vtkProperty.h>

/** Includes - Qt */
#include <QList>

#include "data/DataSet.h"

#include "IsosurfaceVisualization.h"
#include <vtkImageSlice.h>
#include <vtkProp3D.h>

#include <vtkTexture.h>

// Forward declariation
class IsosurfaceVisualization;

namespace bmia
{

/**
This class allows the user to
*/
class vtkInteractorStyleTrackballPositionPicker : public vtkInteractorStyleTrackballCamera
{

	public:
		/**
		Constructor
		*/
		vtkInteractorStyleTrackballPositionPicker();

		/**
		Destructor
		*/
		~vtkInteractorStyleTrackballPositionPicker();

		/**
		Constructor Call
		*/
		static vtkInteractorStyleTrackballPositionPicker* New();

		/**
		VTK "Type" macros.
		*/
		vtkTypeMacro(vtkInteractorStyleTrackballPositionPicker, vtkInteractorStyleTrackballCamera);

		/**
		Set up renderer
		*/
		void SetRenderProcess(vtkRenderer * renderer);

		/**
		Set parent class
		*/
		void SetParentClass(void* parentClass);

		/**
		Rewrite the OnLeftButtonDown function inherited
		*/
		virtual void OnLeftButtonDown();

		/**
		The vtkRenderer
		*/
		vtkRenderer			* Renderer;

		/**
		Reference to parent class for passing through position data
		*/
        IsosurfaceVisualization * parentClass;
};

}//namespace bmia
#endif
