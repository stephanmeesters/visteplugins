OPTION(PLUGINS_BUILD_CAMERATOOLS "Build Camera Tools Plugin" OFF)

if(PLUGINS_BUILD_CAMERATOOLS)

	set( BMIA_CAMERATOOLS_PLUGIN )

	FIND_PACKAGE(VTK REQUIRED)
	INCLUDE(${VTK_USE_FILE})

	FIND_PACKAGE(Qt4 REQUIRED)
	INCLUDE(${QT_USE_FILE})

	SET (SRCS
		CameraTools.cxx
		CameraTools.h
	)

	qt4_wrap_cpp( MOC_SRCS CameraTools.h )

	add_library( CameraToolsPlugin SHARED ${SRCS} ${MOC_SRCS} )
	target_link_libraries( CameraToolsPlugin bmia_plugin bmia_data vtkVolumeRendering vtkIO vtkWidgets vtkGraphics vtkRendering )

endif(PLUGINS_BUILD_CAMERATOOLS)