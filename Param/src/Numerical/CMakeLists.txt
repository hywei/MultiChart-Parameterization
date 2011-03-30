file(GLOB HEADERS *.h)
file(GLOB SOURCES *.cpp)

include_directories( ${PROJECT_SOURCE_DIR}/include)
include_directories( ${PROJECT_SOURCE_DIR}/include/hj_3rd)
include_directories( ${PROJECT_SOURCE_DIR}/src/Graphite)

add_library(numerical STATIC ${HEADERS} ${SOURCES})
