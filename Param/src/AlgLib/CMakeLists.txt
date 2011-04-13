#file(GLOB HEADERS *.h)
#file(GLOB SOURCES *.cpp)

set( HEADERS
  ap.h
  alglibinternal.h
  alglibmisc.h
  linalg.h
  solvers.h
  optimization.h
)

set( SOURCES
  ap.cpp
  alglibinternal.cpp
  alglibmisc.cpp
  linalg.cpp
  solvers.cpp
  optimization.cpp
)

add_library(alglib STATIC ${HEADERS} ${SOURCES})