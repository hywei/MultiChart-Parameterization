file(GLOB HEADERS *.h)
file(GLOB SOURCES *.cpp)

add_library(common STATIC ${HEADERS} ${SOURCES})
