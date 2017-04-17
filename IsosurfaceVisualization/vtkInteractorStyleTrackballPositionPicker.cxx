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

/** Includes */

#include "vtkInteractorStyleTrackballPositionPicker.h"

namespace bmia
{


//-----------------------------[ Constructor ]-----------------------------\\

vtkInteractorStyleTrackballPositionPicker::vtkInteractorStyleTrackballPositionPicker()
{
	// Initial setting for parameters
}


//------------------------------[ Destructor ]-----------------------------\\

vtkInteractorStyleTrackballPositionPicker::~vtkInteractorStyleTrackballPositionPicker()
{
	// Destroy

}


vtkStandardNewMacro(vtkInteractorStyleTrackballPositionPicker);


//--------------------------[ SetRenderProcess ]---------------------------\\

void vtkInteractorStyleTrackballPositionPicker::SetRenderProcess(vtkRenderer * renderer)
{
	// Get the vtkRenderer
	this->Renderer = renderer;
}


//---------------------------[ SetParentClass ]-----------------------\\

void vtkInteractorStyleTrackballPositionPicker::SetParentClass(void* parentClass)
{
	// Get the parent class
	this->parentClass = static_cast<IsosurfaceVisualization*>(parentClass);
}


//-----------------------------[ OnLeftButtonDown ]------------------------\\

void vtkInteractorStyleTrackballPositionPicker::OnLeftButtonDown()
{
    // Forward events
	vtkInteractorStyleTrackballCamera::OnLeftButtonDown();

	// Get the cellpicker from the interactor
	vtkCellPicker *CellPicker = vtkCellPicker::SafeDownCast(Interactor->GetPicker());
	CellPicker->PickFromListOff();

	// Perform pick action on the click event position.Return 1 is successfully picked, otherwise 0.
	int Picked = CellPicker->Pick(this->Interactor->GetEventPosition()[0],
								   this->Interactor->GetEventPosition()[1],
								   0,  // always zero
								   this->Renderer);

	// Perform if pick action successed.
	if(Picked)
	{
        /*double pickPos[3];
        CellPicker->GetPickPosition( pickPos );
        double xp = pickPos[0];
        double yp = pickPos[1];
        double zp = pickPos[2];*/

        /*vtkTexture* tex = CellPicker->GetTexture();

        if(tex == NULL)
            std::cout << "no texture!" << std::endl;
        else
            tex->Print(std::cout);*/

        vtkPoints* points = CellPicker->GetPickedPositions();

        vtkIdType VTKId = 0;

        for(int i = 0; i<points->GetNumberOfPoints(); i++)
        {
            double pickPos[3];
            points->GetPoint(VTKId,pickPos);

            std::cout << "-----------------------" << std::endl;

            for(int k = 0; k<3; k++)
                std::cout << "picked " << pickPos[k] << std::endl;

            this->parentClass->setClippingPlanesPosition(pickPos);
        }


	}
}

} // namespace bmia
