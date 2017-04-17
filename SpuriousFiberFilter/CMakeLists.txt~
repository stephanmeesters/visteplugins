OPTION(PLUGINS_BUILD_SCORINGMEASURES "Build Scoring Measures Plugin" OFF)

if(PLUGINS_BUILD_SCORINGMEASURES)

	set( BMIA_SCORINGMEASURES_PLUGIN )

	FIND_PACKAGE(VTK REQUIRED)
	INCLUDE(${VTK_USE_FILE})
	FIND_PACKAGE(Qt4 REQUIRED)
	INCLUDE(${QT_USE_FILE})

	SET (SRCS
		ScoringMeasures.cxx
		ScoringMeasures.h
		ScoringMeasuresTypes.h
		vtkFiberScoringMeasuresFilter.cxx
		vtkFiberScoringMeasuresFilter.h
	)

	include_directories(${CMAKE_CURRENT_BINARY_DIR})
	include_directories(${BMIA_LIBS_SOURCE_DIR})

	QT4_WRAP_UI(UiHeaders ScoringMeasures.ui)
	SET_SOURCE_FILES_PROPERTIES(ScoringMeasures.cxx PROPERTIES OBJECT_DEPENDS "${UiHeaders}")

	QT4_WRAP_CPP(MOC_SRCS ScoringMeasures.h)
	ADD_LIBRARY(ScoringMeasures SHARED ${SRCS} ${MOC_SRCS})

	TARGET_LINK_LIBRARIES(ScoringMeasures bmia_core bmia_plugin bmia_data vtkHybrid vtkCharts)

endif(PLUGINS_BUILD_SCORINGMEASURES)
