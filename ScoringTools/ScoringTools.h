#ifndef bmia_ScoringTools_h
#define bmia_ScoringTools_h

/** Includes - Main Header */

#include "DTITool.h"

/** Includes - Qt */

#include "ui_ScoringTools.h"

/** Includes - VTK */

#include <vtkPropAssembly.h>
#include <vtkProp3D.h>
#include <vtkActor.h>
#include <vtkActor2D.h>
#include <vtkPolyData.h>
#include <vtkStringArray.h>
#include <vtkVertexGlyphFilter.h>
#include <vtkPointSetToLabelHierarchy.h>
#include <vtkLabelPlacementMapper.h>
#include <vtkPoints.h>
#include <vtkSmartPointer.h>
#include <vtkPointData.h>
#include <vtkPolyDataMapper.h>
#include <vtkDoubleArray.h>
#include <vtkContextView.h>
#include <vtkContextScene.h>
#include <vtkChartXY.h>
#include <vtkTable.h>
#include <vtkIntArray.h>
#include <vtkFloatArray.h>
#include <vtkPlot.h>
#include <vtkAxis.h>

/** Includes - Custom filters */
#include "vtkFiberSelectionFilter.h"
#include "vtkFiberROICutting.h"

/** Includes - Common types */
#include "ScoringTypes.h"

namespace Ui
{
    class ScoringToolsForm;
}

namespace bmia
{

//
//  STRUCTS
//

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
    data::DataSet* ds_processed;

    QList<FiberData*> selectedLines;
    int userSelectedLine;
    int selectedScalarType;
    int numberOfScalarTypes;

    bool processed;
    bool hasScalars;

    QList<ThresholdSettings*> scalarThresholdSettings;

    double lengthOfFiberRange[2];
    double lengthOfFiberSetting[2];

    int prunePercentage;

    QString outputFiberDataName;

} SortedFibers;

//
//  CLASS
//

class ScoringTools :  public plugin::AdvancedPlugin,
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

    ScoringTools();
    ~ScoringTools();

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

protected slots:

    void fibersComboChanged(int index);
    void scalarTypeComboChanged(int index);
    void averageValueMinSliderChanged(int value);
    void averageValueMinSpinBoxChanged(double value);
    void averageValueMaxSliderChanged(int value);
    void averageValueMaxSpinBoxChanged(double value);
    void minkowskiAverageValueMinSliderChanged(int value);
    void minkowskiAverageValueMinSpinBoxChanged(double value);
    void minkowskiAverageValueMaxSliderChanged(int value);
    void minkowskiAverageValueMaxSpinBoxChanged(double value);
    void updateButtonClicked();
    void outputLineEditChanged(QString text);
    void setActiveScalarsButtonClicked();
    void fiberLengthSliderChanged(int value);
    void fiberLengthSpinBoxChanged(double value);
    void globalMinimumSliderChanged(int value);
    void globalMinimumSpinBoxChanged(double value);
    void globalMaximumSliderChanged(int value);
    void globalMaximumSpinBoxChanged(double value);
    void displayHistogramButtonClicked();
    void displayOutputHistogramButtonClicked();
    void minkowskiOrderSpinBoxChanged(int value);
    void prunePercentageSpinBoxChanged(int value);

private:

    /** If plugin inherits from plugin::GUI */
    QWidget * widget;

    /** QT form */
    Ui::ScoringToolsForm * form;

    /** The collection of all the actors that this plugin can render.
    		This is the object that will be returned by getVtkProp().  */
    vtkPropAssembly * assembly;

    //
    //  Qt communication
    //

    /** Connect GUI controls to their respective "SLOT" functions. */
    void connectAll();

    /** Disconnect all GUI controls */
    void disconnectAll();

    //
    //  Plugin base
    //

    enum HistogramType
    {
        HISTOGRAM_INPUT = 0,
        HISTOGRAM_OUTPUT
    };

    /** Fiber information structs */
    QList<SortedFibers*> sortedFibersList;

    /** Find a fiber data struct in the sortedFibersList **/
    int FindInputDataSet(data::DataSet * ds);

    void SelectFiberDataSet(int index);

    int selectedFiberDataset;

    void SelectScalarType(int index);
    void ComputeFibers();

    void BlockSignals();
    void AllowSignals();

    void UpdateGUI();
    SortedFibers* GetSortedFibers();
    ThresholdSettings* GetThresholdSettings();

    void EnableGUI();
    void DisableGUI();
    void SetActiveScalars();

    void ComputeFiberLengthRange();

    //QWidget* histogramWindow;
    void ShowHistogram(HistogramType histType);

    /** ROI scalar data list */
    QList<data::DataSet *> roiDataSets;
};

}

#endif  // bmia_ScoringTools_h
