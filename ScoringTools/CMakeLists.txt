OPTION(PLUGINS_BUILD_SCORINGTOOLS "Build Scoring Tools Plugin" OFF)

if(PLUGINS_BUILD_SCORINGTOOLS)

	set( BMIA_SCORINGTOOLS_PLUGIN )

	FIND_PACKAGE(VTK REQUIRED)
	INCLUDE(${VTK_USE_FILE})
	FIND_PACKAGE(Qt4 REQUIRED)
	INCLUDE(${QT_USE_FILE})

	SET (SRCS
		ScoringTools.cxx
		ScoringTools.h
		vtkFiberSelectionFilter.cxx
		vtkFiberSelectionFilter.h
		ScoringTypes.h
		vtkFiberROICutting.cxx
		vtkFiberROICutting.h
	)

	include_directories(${CMAKE_CURRENT_BINARY_DIR})
	include_directories(${BMIA_LIBS_SOURCE_DIR})

	QT4_WRAP_UI(UiHeaders ScoringTools.ui)
	SET_SOURCE_FILES_PROPERTIES(ScoringTools.cxx PROPERTIES OBJECT_DEPENDS "${UiHeaders}")

	QT4_WRAP_CPP(MOC_SRCS ScoringTools.h)
	ADD_LIBRARY(ScoringTools SHARED ${SRCS} ${MOC_SRCS})

	TARGET_LINK_LIBRARIES(ScoringTools bmia_core bmia_plugin bmia_data vtkHybrid vtkCharts)

endif(PLUGINS_BUILD_SCORINGTOOLS)
