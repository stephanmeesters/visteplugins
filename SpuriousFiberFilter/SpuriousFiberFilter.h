#ifndef bmia_SpuriousFiberFilter_h
#define bmia_SpuriousFiberFilter_h

/** Includes - Main Header */

#include "DTITool.h"

/** Includes - Qt */

#include "ui_SpuriousFiberFilter.h"

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
#include <vtkPolyDataWriter.h>
#include <vtkCellArray.h>

/** Includes - Custom */
#include "vtkFiberSpuriousFilter.h"
#include "SpuriousFiberFilterTypes.h"
#include "vtkFiberSelectAnterior.h"

#include <QString>
#include <QFileInfo>
#include <stdlib.h>

namespace Ui
{
    class SpuriousFiberFilterForm;
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
    ParameterSettings* ps;
    vtkPolyData* preprocessedPolyData;

} SortedFibers;

//
//  CLASS
//

class SpuriousFiberFilter :  public plugin::AdvancedPlugin,
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

    SpuriousFiberFilter();
    ~SpuriousFiberFilter();

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
    void updateButtonClicked();
    void D33SliderChanged(int value);
    void D33SpinBoxChanged(double value);
    void D44SliderChanged(int value);
    void D44SpinBoxChanged(double value);
    void TimeSliderChanged(int value);
    void TimeSpinBoxChanged(double value);
    void WindowSizeSliderChanged(int value);
    void WindowSizeSpinBoxChanged(double value);
    void EpsilonSliderChanged(int value);
    void EpsilonSpinBoxChanged(double value);
    void applyKernelInBothDirsCheckBoxChanged(bool checked);
    void minDistSliderChanged(int value);
    void minDistSpinBoxChanged(double value);
    void applySubsamplingCheckBoxChanged(bool checked);
    void subsampleSpinBoxChanged(double value);
    void processBatchButtonClicked();
    void applySelectAnteriorCheckBoxChanged(bool checked);
    void anteriorSelectSpinBoxChanged(int value);
    void cutoffSliderChanged(int value);
    void cutoffSpinBoxChanged(double value);
    void outputLineEditChanged(QString value);

private:

    /** If plugin inherits from plugin::GUI */
    QWidget * widget;

    /** QT form */
    Ui::SpuriousFiberFilterForm * form;

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

    /** Find a fiber data struct in the sortedFibersList **/
    int FindInputDataSet(data::DataSet * ds);

    int selectedFiberDataset;

    void BlockSignals();
    void AllowSignals();
    void UpdateGUI();
    void EnableGUI();
    void DisableGUI();

    void SelectFiberDataSet(int index);
    void ComputeScore();

    SortedFibers* GetSortedFibers();
    ParameterSettings* GetParameterSettings();
};

}

#endif  // bmia_SpuriousFiberFilter_h
