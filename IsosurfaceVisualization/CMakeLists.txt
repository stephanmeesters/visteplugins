OPTION(PLUGINS_BUILD_ISOSURFACEVISUALIZATION "Build Isosurface Visualization Plugin" OFF)

if(PLUGINS_BUILD_ISOSURFACEVISUALIZATION)

	set( BMIA_ISOSURFACEVISUALIZATION_PLUGIN )

	FIND_PACKAGE(VTK REQUIRED)
	INCLUDE(${VTK_USE_FILE})
	FIND_PACKAGE(Qt4 REQUIRED)
	INCLUDE(${QT_USE_FILE})

	IF(NOT VTK_USE_CHARTS)
		MESSAGE(FATAL_ERROR "Isosurface Visualization Plugin requires VTK_USE_CHARTS.")
	ENDIF(NOT VTK_USE_CHARTS)

	SET (SRCS
		IsosurfaceVisualization.cxx
		IsosurfaceVisualization.h
		vtkInteractorStyleTrackballPositionPicker.cxx
		vtkInteractorStyleTrackballPositionPicker.h
		vtkCurvaturesShapeIndex.cxx
		vtkCurvaturesShapeIndex.h
	)

	include_directories(${CMAKE_CURRENT_BINARY_DIR})
	include_directories(${BMIA_LIBS_SOURCE_DIR})

	QT4_WRAP_UI(UiHeaders IsosurfaceVisualization.ui)
	SET_SOURCE_FILES_PROPERTIES(IsosurfaceVisualization.cxx PROPERTIES OBJECT_DEPENDS "${UiHeaders}")

	QT4_WRAP_CPP(MOC_SRCS IsosurfaceVisualization.h)
	ADD_LIBRARY(IsosurfaceVisualization SHARED ${SRCS} ${MOC_SRCS})

	TARGET_LINK_LIBRARIES(IsosurfaceVisualization bmia_core bmia_plugin bmia_data vtkHybrid vtkCharts)

endif(PLUGINS_BUILD_ISOSURFACEVISUALIZATION)
