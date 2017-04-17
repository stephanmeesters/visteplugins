OPTION(PLUGINS_BUILD_DISTANCEMEASURES "Build Distance Measures Plugin" OFF)

if(PLUGINS_BUILD_DISTANCEMEASURES)

	set( BMIA_DISTANCEMEASURES_PLUGIN )

	FIND_PACKAGE(VTK REQUIRED)
	INCLUDE(${VTK_USE_FILE})
	FIND_PACKAGE(Qt4 REQUIRED)
	INCLUDE(${QT_USE_FILE})

	IF(NOT VTK_USE_CHARTS)
		MESSAGE(FATAL_ERROR "Distance Measures Plugin requires VTK_USE_CHARTS.")
	ENDIF(NOT VTK_USE_CHARTS)

	SET (SRCS
		DistanceMeasures.cxx
		DistanceMeasures.h
	)

	include_directories(${CMAKE_CURRENT_BINARY_DIR})
	include_directories(${BMIA_LIBS_SOURCE_DIR})

	QT4_WRAP_UI(UiHeaders DistanceMeasures.ui)
	SET_SOURCE_FILES_PROPERTIES(DistanceMeasures.cxx PROPERTIES OBJECT_DEPENDS "${UiHeaders}")

	QT4_WRAP_CPP(MOC_SRCS DistanceMeasures.h)
	ADD_LIBRARY(DistanceMeasures SHARED ${SRCS} ${MOC_SRCS})

	TARGET_LINK_LIBRARIES(DistanceMeasures bmia_core bmia_plugin bmia_data vtkHybrid vtkCharts)

endif(PLUGINS_BUILD_DISTANCEMEASURES)
