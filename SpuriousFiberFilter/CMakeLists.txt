OPTION(PLUGINS_BUILD_SPURIOUSFIBERFILTER "Build Spurious Fiber Filter Plugin" ON)

if(PLUGINS_BUILD_SPURIOUSFIBERFILTER)

	set( BMIA_SPURIOUSFIBERFILTER_PLUGIN )

	FIND_PACKAGE(VTK REQUIRED)
	INCLUDE(${VTK_USE_FILE})
	FIND_PACKAGE(Qt4 REQUIRED)
	INCLUDE(${QT_USE_FILE})

	SET (SRCS
		SpuriousFiberFilter.cxx
		SpuriousFiberFilter.h
		SpuriousFiberFilterTypes.h
		vtkFiberSpuriousFilter.cxx
		vtkFiberSpuriousFilter.h
		vtkFiberSelectAnterior.cxx
		vtkFiberSelectAnterior.h
	)

	include_directories(${CMAKE_CURRENT_BINARY_DIR})
	include_directories(${BMIA_LIBS_SOURCE_DIR})

	QT4_WRAP_UI(UiHeaders SpuriousFiberFilter.ui)
	SET_SOURCE_FILES_PROPERTIES(SpuriousFiberFilter.cxx PROPERTIES OBJECT_DEPENDS "${UiHeaders}")

	QT4_WRAP_CPP(MOC_SRCS SpuriousFiberFilter.h)
	ADD_LIBRARY(SpuriousFiberFilter SHARED ${SRCS} ${MOC_SRCS})

	TARGET_LINK_LIBRARIES(SpuriousFiberFilter bmia_core bmia_plugin bmia_data vtkHybrid vtkCharts)

endif(PLUGINS_BUILD_SPURIOUSFIBERFILTER)
