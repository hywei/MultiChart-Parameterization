file(GLOB HEADERS *.h)
file(GLOB SOURCES *.cpp)

find_package(OpenGL)
find_package(Qt4 COMPONENTS QtCore QtGui QtOpenGL REQUIRED)

include(${QT_USE_FILE})

include_directories( $ENV{QTDIR}/include
                     $ENV{QTDIR}/include/QtCore
                     $ENV{QTDIR}/include/QtGui
                     $ENV{QTDIR}/include/QtOpenGl
                     ${Boost_INCLUDE_DIR}
                     ${PROJECT_SOURCE_DIR}/include
                   )

link_directories( $ENV{QTDIR}/lib
                  ${PROJECT_SOURCE_DIR}/lib
                )

QT4_WRAP_CPP(MOCFILES ${HEADERS})

LINK_DIRECTORIES(/usr/local/lib/hj)
add_executable( MainWindow ${HEADERS} ${SOURCES} ${MOCFILES})

set ( DEPENDENCIES ui
                   opengl
                   param
                   meshmodel
                   numerical
                   common
                   graphite
                   alglib
                    )

add_dependencies (
  alglib
  graphite
  common
  numerical
  meshmodel
  param
  opengl
  ui
  ) 
                                      

if(WIN32)
	if(MSVC)
		set ( OPENGL_LIBRARIES opengl32.lib glu32.lib glaux.lib)
        set ( NUMERIC_LIBRARIES cblas.lib lapack.lib linalg.lib
              sparseRelease.lib sparse.lib)
        set ( QT_USE_LIBRARIES qtmain QtCore4 QtGui4 QtOpenGL4)
	endif(MSVC)
else ()
	set ( OPENGL_LIBRARIES libGL.so libGLU.so)
    set ( NUMERIC_LIBRARIES blas lapack sparse)
    set ( HJ_3RD_LIBRARIES hj-sparse-solver hj-sparse-util)
    set ( QT_USE_LIBRARIES QtCore QtGui QtOpenGL )
endif ()


target_link_libraries( MainWindow
                       ${OPENGL_LIBRARIES}
                       ${QT_USE_LIBRARIES}
                       ${NUMERIC_LIBRARIES}
                       ${HJ_3RD_LIBRARIES}
                       ${DEPENDENCIES}
                       quad_dof
                     )
                  
        
                  
