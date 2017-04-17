#ifndef bmia_ConnectivityMeasurePlugin_vtkFiberScoringMeasuresFilter_h
#define bmia_ConnectivityMeasurePlugin_vtkFiberScoringMeasuresFilter_h


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
#include <vtkImageData.h>

/** Includes - Qt */

#include <QMap>

/** Includes - Custom Files */

#include "ScoringMeasuresTypes.h"

namespace bmia {

/** This class is used to
*/

class vtkFiberScoringMeasuresFilter : public vtkPolyDataToPolyDataFilter
{
	public:

		/** Constructor Call */

		static vtkFiberScoringMeasuresFilter * New();

		/** VTK Macro */

		vtkTypeMacro(vtkFiberScoringMeasuresFilter, vtkPolyDataToPolyDataFilter);

		/** Set a new input volume.
			@param image		New discrete sphere function volume. */

		void SetInputVolume(vtkImageData * image);

		/** Set parameter settings
			@param ps		Parameter settings struct. */

		void SetParameters(ParameterSettings* ps);

        /** Method for curve generation **/

		enum TypeOfCurve
        {
            CURVE_TYPE_GEODESIC = 0,
            CURVE_TYPE_ELASTICA
        };


	protected:

		/** Main entry point of the filter. */

		virtual void Execute();

		/** Constructor. */

		vtkFiberScoringMeasuresFilter();

		/** Destructor. */

		~vtkFiberScoringMeasuresFilter();

		/** Input volume, containing an array defining the spherical directions,
			an array with the radius per direction per voxel, and (optionally)
			a triangles array defining the topology of the glyphs. */

		vtkImageData * inputVolume;

		ParameterSettings* ps;

}; // class vtkFiberScoringMeasuresFilter


} // namespace bmia


#endif // bmia_ConnectivityMeasurePlugin_vtkFiberScoringMeasuresFilter_h
