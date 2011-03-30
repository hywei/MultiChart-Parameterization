find_package(OpenGL)
find_package(Qt4 COMPONENTS QtCore QtGui QtOpenGL REQUIRED)

include(${QT_USE_FILE})


include_directories( $ENV{QTDIR}/include
                     $ENV{QTDIR}/include/QtOpenGL
                     ${Boost_INCLUDE_DIR}
                   )

file(GLOB HEADERS *.h)
file(GLOB SOURCES *.cpp)

add_library(ui STATIC ${HEADERS} ${SOURCES})
