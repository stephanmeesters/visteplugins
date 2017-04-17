#ifndef bmia_ConnectivityMeasurePlugin_vtkFiberROICutting_h
#define bmia_ConnectivityMeasurePlugin_vtkFiberROICutting_h

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

class vtkFiberROICutting : public vtkPolyDataToPolyDataFilter
{
	public:

		/** Constructor Call */

		static vtkFiberROICutting * New();

		/** VTK Macro */

		vtkTypeMacro(vtkFiberROICutting, vtkPolyDataToPolyDataFilter);

        void SetROIData(data::DataSet * d)
        {
            this->roiData = d;
        }

        void SetFiberTransformationMatrix(vtkMatrix4x4* m)
        {
            this->fiberMatrix = m;
        }

        void SetCutFibersAtROIEnds(bool b)
        {
            this->cutAtROIEnds = b;
        }

	protected:

		/** Main entry point of the filter. */

		virtual void Execute();

		/** Constructor. */

		vtkFiberROICutting();

		/** Destructor. */

		~vtkFiberROICutting();

		/** ROI data */
        data::DataSet * roiData;

        /** Fiber transformation matrix */
        vtkMatrix4x4* fiberMatrix;

        /** Option: cut at ROI ends */
        bool cutAtROIEnds;

        /** Selected scalar value type */

		int scalarType;

}; // class vtkFiberROICutting


} // namespace bmia


#endif // bmia_ConnectivityMeasurePlugin_vtkFiberROICutting_h
