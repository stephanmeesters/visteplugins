#ifndef bmia_IsosurfaceVisualization_h
#define bmia_IsosurfaceVisualization_h

/** Includes - Main Header */

#include "DTITool.h"

/** Includes - Qt */

#include "ui_IsosurfaceVisualization.h"

/** Includes - STL */

#include <algorithm>

/** Includes - VTK */

#include <vtkPropAssembly.h>
#include <vtkActor.h>
#include <vtkSphereSource.h>
#include <vtkPolyDataMapper.h>
#include <vtkSmartPointer.h>
#include <vtkImageThreshold.h>
#include <vtkImageGaussianSmooth.h>
#include <vtkMarchingCubes.h>
#include <vtkLODActor.h>
#include <vtkProperty.h>
#include <vtkMatrix4x4.h>
#include <vtkWindowedSincPolyDataFilter.h>
#include <vtkDoubleArray.h>
#include <vtkPointData.h>
#include <vtkClipPolyData.h>
#include <vtkPlane.h>
#include <vtkDecimatePro.h>
#include <vtkSmoothPolyDataFilter.h>
#include <vtkPlaneSource.h>
#include <vtkBox.h>
#include <vtkCubeSource.h>
#include <vtkSuperquadricSource.h>
#include <vtkProp3D.h>
#include <vtkOrientationMarkerWidget.h>
#include <vtkLookupTable.h>
#include <vtkImageActor.h>
#include <vtkImageMapToColors.h>
#include <vtkImageResliceMapper.h>
#include <vtkImageSlice.h>
#include <vtkImageProperty.h>
#include <vtkTransform.h>
#include <vtkImageMapToColors.h>
#include <vtkExtractPolyDataGeometry.h>
#include <vtkPolyDataConnectivityFilter.h>
#include <vtkColorTransferFunction.h>
#include <vtkPiecewiseFunction.h>
#include <vtkDiscretizableColorTransferFunction.h>
#include <vtkVersion.h>
#include <vtkDepthSortPolyData.h>
#include <vtkLineSource.h>
#include <vtkObject.h>
#include <vtkPolyDataWriter.h>
#include <vtkPolyData.h>
#include <vtkPoints.h>
#include <vtkSelectEnclosedPoints.h>
#include <vtkDiskSource.h>
#include <vtkPolyDataMapper2D.h>
#include <vtkActor2D.h>
#include <vtkProperty2D.h>
#include <vtkColorSeries.h>
#include <vtkThreshold.h>
#include <vtkPolyDataNormals.h>
#include <vtkPointData.h>
#include <vtkContextView.h>
#include <vtkContextScene.h>
#include <vtkChart.h>
#include <vtkChartXY.h>
#include <vtkPlot.h>
#include <vtkTable.h>
#include <vtkQtTableView.h>
#include <vtkTimerLog.h>
#include <vtkFloatArray.h>
#include <vtkAxis.h>
#include <vtkScalarBarActor.h>
#include <vtkColorSeries.h>
#include <vtkSmoothPolyDataFilter.h>
#include <vtkContextScene.h>
#include <vtkContextItem.h>
#include <vtkContextMouseEvent.h>
#include <vtkXMLMaterial.h>
#include <vtkOpenGLExtensionManager.h>
#include <vtkgl.h>

/** Includes - Custom */

#include "gui/MetaCanvas/vtkSubCanvas.h"
#include "gui/MetaCanvas/vtkMedicalCanvas.h"
#include "vtkInteractorStyleTrackballPositionPicker.h"
#include "Helpers/vtkImageSliceActor.h"
#include "vtkCurvaturesShapeIndex.h"

class vtkInteractorStyleTrackballPositionPicker;

namespace Ui
{
    class IsosurfaceVisualizationForm;
}

namespace bmia
{

/** ModelInfo struct holding all the settings of a generated mesh */
typedef struct
{
    // raw data
    data::DataSet* ds;
    double scalarRange[2];
    int dataDimensions[3];

    // polydata DataSet
    data::DataSet* ds_poly;
    bool bDsPolyAdded;

    // mesh and mesh settings
    vtkPolyData* polydata;
	vtkPolyDataMapper* polyDataMapper;
    vtkLODActor* prop;
    double minimumThreshold;
    double maximumThreshold;
    double smoothing;
    double alpha;
    QString colorString;
    double color[3];
    bool visible;
    int renderStyle;
    double reduction;
    double specular;
    double global_bounds[6];
    bool selectLargestComponent;

    // clipping cube
    vtkBox* clippingCube;
    double bounds[6];
    double clippingValues[3];
    bool flipped[3];
    vtkExtractPolyDataGeometry* extractPolyFunc;
    bool invertClipping;
    bool alignPlanesToPick;

    // orthogonal planes
    QList<vtkImageSlice*> orthogonalPlanesModels;
    QList<vtkPlane*> orthogonalPlanes;
    bool planesActivated[3];
    vtkLookupTable* defaultLUT;

    // base layer
    data::DataSet* dsBaseLayer;
    int baseLayerLUTIndex;

    // overlay
    data::DataSet* dsOverlay;
    int overlayIndex;
    int overlayLUTIndex;

    // fiber selection
    int selectedFiberData;

	// curvature
	int curvatureLUTIndex;
	int curvatureType;
	bool usingCurvature;

} ModelInfo;

class IsosurfaceVisualization :  public plugin::AdvancedPlugin,
                    public plugin::Visualization,
                    public plugin::GUI,
                    public data::Consumer
{
    Q_OBJECT
    Q_INTERFACES(bmia::plugin::Plugin)
    Q_INTERFACES(bmia::plugin::AdvancedPlugin)
    Q_INTERFACES(bmia::plugin::GUI)
    Q_INTERFACES(bmia::data::Consumer)
    Q_INTERFACES(bmia::plugin::Visualization)

public:

    QString getPluginVersion()
    {
        return "1.0.0";
    }

    IsosurfaceVisualization();
    ~IsosurfaceVisualization();

    virtual void init();

    /**
     * Return the VTK actor that renders the cone.
     */
    vtkProp * getVtkProp();

    /**
     * Return the widget that is shown in the GUI
     */
    QWidget * getGUI();

    /**
     * This function is called when a new data set becomes available.
     *
     * @param ds The new data set that was added.
     */
    void dataSetAdded(data::DataSet * d);

    /**
     * This function is called when an already available data set was changed.
     *
     * @param ds The data set that has been updated.
     */
    void dataSetChanged(data::DataSet * d);

    /**
     * This function is called when a data set that was available has been removed.
     *
     * @param ds The data set that was removed from the pool of data sets.
     */
    void dataSetRemoved(data::DataSet * d);

    /**
     * Set the clipping plane position
     *
     * @param pos Array holding the position coordinates.
     */
    void setClippingPlanesPosition(double* pos);

protected slots:

    void comboBoxDataChanged();
    void comboBoxOverlayChanged();
    void comboBoxOverlayLUTChanged();
    void comboBoxBaseLayerLUTChanged();
    void checkBoxVisibleChanged(bool checked);
    void comboBoxStyleChanged();
    void buttonUpdateClicked();
    void inputMaximumThresholdChanged(double value);
    void inputMinimumThresholdChanged(double value);
    void inputSmoothingChanged(double value);
    void inputColorChanged();
    void inputAlphaChanged(double value);
    void inputReductionChanged(double value);
    void inputSpecularChanged(double value);
	void comboBoxCurvatureLUTChanged();
	void comboBoxCurvatureTypeChanged();
    void horizontalSliderXChanged(int value);
    void horizontalSliderYChanged(int value);
    void horizontalSliderZChanged(int value);
    void checkBoxXChanged(bool checked);
    void checkBoxYChanged(bool checked);
    void checkBoxZChanged(bool checked);
    void checkBoxFlipXChanged(bool checked);
    void checkBoxFlipYChanged(bool checked);
    void checkBoxFlipZChanged(bool checked);
    void buttonSaveMeshClicked();
    void spinXChanged(int value);
    void spinYChanged(int value);
    void spinZChanged(int value);
    void checkBoxInvertClippingChanged(bool checked);
    void checkBoxLargestComponentChanged(bool checked);
    void checkBoxAlignPlanesToPickChanged(bool checked);

private:

    /** If plugin inherits from plugin::GUI */
    QWidget * widget;

    /** QT form */
    Ui::IsosurfaceVisualizationForm * form;

    /** Holds the scalar volume data sets */
    QList<data::DataSet *> datasets;

    /** Holds the transfer function data sets */
    QList<data::DataSet *> tf_datasets;

    /** Transfer functions for orthogonal planes */
    QList<vtkColorTransferFunction*> colorTransferFunctions;

    /** Color loopup tables */
    QList<vtkLookupTable*> lookUpTables;

    /** The collection of all the actors that this plugin can render.
    		This is the object that will be returned by getVtkProp().  */
    vtkPropAssembly * assembly;

    /** Collection of all the ModelInfo structs */
    QList<ModelInfo*> modelInfoList;

    /** ModelInfo struct of current model */
    ModelInfo* current_modelInfo;

    /** ModelInfo struct of base layer plane */
    ModelInfo* baseLayer_modelInfo;

    /** ModelInfo struct of overlay plane */
    ModelInfo* overlay_modelInfo;

    /** Reconstruct mesh on call to updateRenderingModels()  */
    bool bModelDirty;

    double clickedPoint[3];

    int findInputDataSet(data::DataSet * ds);

    /** Connect GUI controls to their respective "SLOT" functions. */
    void connectAll();

    /** Disconnect all GUI controls */
    void disconnectAll();

    /** Disconnect all GUI controls */
    void updateRenderingModels();

    /**
     * Create a model info struct for the new mesh.
     *
     * @param ds The new data set that was added.
     */
    void createModelInfo(data::DataSet * d);

    /**
     * Create a new lookup table based on the loaded transfer function.
     *
     * @param ds The new data set that was added.
     */
    void createLookupTable(data::DataSet * d, int index = -1);

    /**
     *  Update the mesh after adjusting the clipping plane
     *
     *  @param      direction   0:x 1:y 2:z
     *  @param      value       Clipping plane position for that direction
     */
    void updateClippingPlaneSlider(int direction, int value, bool render = true);

    /**
     *  Update the mesh after enabling/disabling the clipping plane
     *
     *  @param      direction   0:x 1:y 2:z
     *  @param      checked     Enabled/disabled
     */
    void updateClippingPlaneEnabled(int direction, bool checked);

    /**
     *  Disable GUI elements to prevent propagation
     */
    void BlockSignals();

    /**
     *  Enable GUI elements
     */
    void AllowSignals();

    /** Position picker */
    void setupClippingPlanesPicker();
    vtkInteractorStyleTrackballPositionPicker * styleTrackballPP;

    // 2d points in medical planes
    QList<vtkActor*> pointer2DList;

    // temporary scalar bar graphic
    vtkSmartPointer<vtkScalarBarActor> scalarBar;

    // dataset holding selected point settings
	data::DataSet* settings;
};

}

#endif  // bmia_IsosurfaceVisualization_h
