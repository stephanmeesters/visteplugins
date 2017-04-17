#ifndef bmia_ScoringMeasures_h
#define bmia_ScoringMeasures_h

/** Includes - Main Header */

#include "DTITool.h"

/** Includes - Qt */

#include "ui_ScoringMeasures.h"

/** Includes - VTK */

#include <vtkPropAssembly.h>
#include <vtkProp3D.h>
#include <vtkActor.h>
#include <vtkActor2D.h>
#include <vtkImageData.h>
#include <vtkPointData.h>
#include <vtkDecimatePro.h>
#include <vtkSmoothPolyDataFilter.h>
#include <vtkSplineFilter.h>

/** Includes - Custom */
#include "vtkFiberScoringMeasuresFilter.h"
#include "ScoringMeasuresTypes.h"

namespace Ui
{
    class ScoringMeasuresForm;
}

namespace bmia
{

//
//  STRUCTS
//

/** Holding fiberdata */
typedef struct
{
    data::DataSet* ds;
    data::DataSet* ds_processed;

    bool processed;

    QString outputFiberDataName;
    int selectedGlyphData;

    ParameterSettings* ps;

    vtkPolyData* preprocessedPolyData;

} SortedFibers;

//
//  CLASS
//

class ScoringMeasures :  public plugin::AdvancedPlugin,
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

    ScoringMeasures();
    ~ScoringMeasures();

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
    void glyphDataComboChanged(int index);
    void updateButtonClicked();
    //void usedInScoringCheckBoxChanged(bool checked);
    void lambdaSliderChanged(int value);
    void lambdaSpinBoxChanged(double value);
    void betaSliderChanged(int value);
    void betaSpinBoxChanged(double value);
    void muuSliderChanged(int value);
    void muuSpinBoxChanged(double value);
    void standardizeScalarsCheckBoxChanged(bool checked);
    void normalizeGlyphDataCheckBoxChanged(bool checked);
    void applyLogarithmCheckBoxChanged(bool checked);

private:

    /** If plugin inherits from plugin::GUI */
    QWidget * widget;

    /** QT form */
    Ui::ScoringMeasuresForm * form;

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

    /** Fiber information structs */
    QList<SortedFibers*> sortedFibersList;

    /** Glyph data list */
    QList<data::DataSet *> glyphDataSets;

    /** Find a fiber data struct in the sortedFibersList **/
    int FindInputDataSet(data::DataSet * ds);

    int selectedFiberDataset;

    void BlockSignals();
    void AllowSignals();
    void UpdateGUI();
    void EnableGUI();
    void DisableGUI();

    void SelectFiberDataSet(int index);
    void SelectGlyphDataSet(int index);
    void ComputeScore();

    SortedFibers* GetSortedFibers();
    ParameterSettings* GetParameterSettings();
};

}

#endif  // bmia_ScoringMeasures_h
