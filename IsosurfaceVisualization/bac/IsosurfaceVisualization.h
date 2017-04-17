#ifndef bmia_IsosurfaceVisualization_h
#define bmia_IsosurfaceVisualization_h

#include "DTITool.h"
#include "ui_IsosurfaceVisualization.h"

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
#include "gui/MetaCanvas/vtkSubCanvas.h"
#include <vtkOrientationMarkerWidget.h>

#include "vtkLookupTable.h"
#include "vtkImageActor.h"
#include "vtkImageMapToColors.h"
#include "vtkImageResliceMapper.h"
#include "vtkImageSlice.h"
#include "vtkImageProperty.h"
#include "vtkTransform.h"
#include "vtkImageMapToColors.h"
#include "vtkExtractPolyDataGeometry.h"
#include "vtkPolyDataConnectivityFilter.h"
#include "vtkColorTransferFunction.h"
#include "vtkPiecewiseFunction.h"
#include "vtkDiscretizableColorTransferFunction.h"
#include "vtkVersion.h"
#include "vtkDepthSortPolyData.h"
#include "vtkLineSource.h"

// picker
#include "vtkInteractorStyleTrackballPositionPicker.h"

// temporary
#include "vtkDebugLeaks.h"

#include "vtkPolyDataWriter.h"
#include "vtkPolyData.h"

#include "vtkPoints.h"
#include "vtkSelectEnclosedPoints.h"

#include "gui/MetaCanvas/vtkMedicalCanvas.h"
#include "vtkDiskSource.h"
#include "vtkPolyDataMapper2D.h"
#include "vtkActor2D.h"
#include "vtkProperty2D.h"

#include "Helpers/vtkImageSliceActor.h"

// labels
#include "vtkStringArray.h"
#include "vtkVertexGlyphFilter.h"
#include "vtkPointSetToLabelHierarchy.h"
#include "vtkLabelPlacementMapper.h"

// curvature
#include "vtkCurvaturesShapeIndex.h"
#include <vtkColorSeries.h>
#include <vtkThreshold.h>
#include <algorithm>
#include <vtkPolyDataNormals.h>

#include <vtkPointData.h>


// qt window test
#include "vtkContextView.h"
#include "vtkContextScene.h"
#include "vtkChart.h"
#include "vtkChartXY.h"
#include "vtkPlot.h"
#include "vtkTable.h"
#include "vtkQtTableView.h"
#include "vtkTimerLog.h"
#include "vtkFloatArray.h"
#include "vtkAxis.h"

#include "vtkScalarBarActor.h"
#include <vtkColorSeries.h>

#include <vtkSmoothPolyDataFilter.h>

#include <vtkContextScene.h>
#include <vtkContextItem.h>
#include <vtkContextMouseEvent.h>

#include <vtkXMLMaterial.h>
#include <vtkOpenGLExtensionManager.h>
#include <vtkgl.h>

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

/** MeasuredPoint struct holding coordinates and 3D actors of a measured point */
typedef struct
{
    bool set;
    double x, y, z;
    vtkActor* sphere;
    vtkActor2D* label;

} MeasuredPoint;

/** Generic vector3 struct */
typedef struct
{
    double x, y, z;
} Vec3;

/** Holding pointdata and location of the most anterior point of a fiber */
typedef struct
{
    QList< Vec3* > data;
    QList< double > scalarData;
    int anteriorPointIndex;
    double userPointRefinement;
} FiberData;

/** Holding fiberdata */
typedef struct
{
    data::DataSet* ds;
    QList<FiberData*> selectedLines;
    int userSelectedLine;
} SortedFibers;

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

    void buttonSetPointAClicked();
    void buttonSetPointBClicked();

    void buttonSetLineColorClicked();
    void buttonSaveMeasurementClicked();

	void lineEditNamePointAChanged(QString value);
	void lineEditNamePointBChanged(QString value);

	void comboBoxFiberDataChanged();

	void fiberSelectUpdate(int value);
	void fiberRefinementUpdate(int value);
	void fiberRefinementUpdate(double value);

private:

    QWidget secondWindow;

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

    /** Fiber information structs */
    QList<SortedFibers*> sortedFibersList;

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
    QList<MeasuredPoint*> measuredPointList;
    vtkActor* measuredLine;
	vtkActor2D* measuredLabels;
	vtkPoints* measuredLabelPoints;
	vtkStringArray* measuredLabelStrings;

    QList<vtkActor*> depthElectrodeBlobs; // temporary
    QColor currentElectrodesColor; // temporary

    void setMeasuredPoint(int id);
    void calculateDistance();

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

	vtkActor2D* GenerateLabels(vtkSmartPointer<vtkPoints> points, vtkSmartPointer<vtkStringArray> labels);

    /** Position picker */
    void setupClippingPlanesPicker();
    vtkInteractorStyleTrackballPositionPicker * styleTrackballPP;

    QList<vtkActor*> pointer2DList;

    void processFiberAnteriorSorting(SortedFibers* sortedFibers);

    void fiberPointSelect();

    vtkSmartPointer<vtkScalarBarActor> scalarBar; //temp
};

}

#endif  // bmia_IsosurfaceVisualization_h
