#ifndef bmia_ConnectivityMeasurePlugin_vtkFiberSelectAnterior_h
#define bmia_ConnectivityMeasurePlugin_vtkFiberSelectAnterior_h

/** Includes - Main Header */

#include "DTITool.h"

/** Includes - STD */

#include <algorithm>
#include <vector>

/** Includes - VTK */

#include <vtkPolyDataToPolyDataFilter.h>
#include <vtkPolyData.h>
#include <vtkObjectFactory.h>
#include <vtkPointData.h>
#include <vtkDoubleArray.h>
#include <vtkPoints.h>
#include <vtkObjectFactory.h>
#include <vtkDataArray.h>
#include <vtkCellArray.h>
#include <vtkCell.h>
#include <vtkIdList.h>
#include <vtkImageCast.h>
#include <vtkObject.h>
#include <vtkMatrix4x4.h>

/** Includes - Qt */

#include <QMap>

namespace bmia {

/** This class is used to
*/

class vtkFiberSelectAnterior : public vtkPolyDataToPolyDataFilter
{
	public:

		/** Constructor Call */

		static vtkFiberSelectAnterior * New();

		/** VTK Macro */

		vtkTypeMacro(vtkFiberSelectAnterior, vtkPolyDataToPolyDataFilter);

        void SetFiberTransformationMatrix(vtkMatrix4x4* m)
        {
            this->fiberMatrix = m;
        }

        void SetNumberOfAnteriorFibers(int b)
        {
            this->numberOfAnteriorFibers = b;
        }

	protected:

		/** Main entry point of the filter. */

		virtual void Execute();

		/** Constructor. */

		vtkFiberSelectAnterior();

		/** Destructor. */

		~vtkFiberSelectAnterior();

        /** Fiber transformation matrix */
        vtkMatrix4x4* fiberMatrix;

        /** Number of anterior fibers */
        int numberOfAnteriorFibers;

        /** Selected scalar value type */
		int scalarType;

}; // class vtkFiberSelectAnterior


} // namespace bmia


#endif // bmia_ConnectivityMeasurePlugin_vtkFiberSelectAnterior_h
