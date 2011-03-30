include_directories( ${PROJECT_SOURCE_DIR}/src/Graphite )

set ( include_files ../Graphite/)
#message("${include_files}")

set ( HEADERS
      OGF/basic/common/common.h
      OGF/basic/common/linkage.h
      OGF/basic/debug/assert.h
      OGF/basic/debug/logger.h
      OGF/basic/os/environment.h
      OGF/basic/types/counted.h
      OGF/basic/types/smart_pointer.h
      OGF/basic/types/types.h
      OGF/math/common/common.h
      OGF/math/numeric/blas.h
      OGF/math/numeric/lapack.h
      OGF/math/symbolic/polynomial.h
      OGF/math/symbolic/stencil.h
      OGF/math/symbolic/symbolic.h
      )
      
set ( SOURCES
      OGF/basic/common/
      OGF/basic/debug/
      OGF/basic/os/
      OGF/basic/types/
      OGF/math/numeric/blas.cpp
      OGF/math/numeric/lapack.cpp
      OGF/math/symbolic/polynomial.cpp
      OGF/math/symbolic/stencil.cpp
      OGF/math/symbolic/symbolic.cpp
     )

#file(GLOB HEADERS OGF/basic/common/*.h
#                  OGF/basic/debug/*.h
#                  OGF/basic/os/*.h
#                  OGF/basic/types/*.h
#                  OGF/math/common/*.h
#                  OGF/math/numeric/*.h
#                  OGF/math/symbolic/*.h)

# message("${HEADERS}")
                 
                  
#file(GLOB SOURCES OGF/basic/common/*.cpp
#                  OGF/basic/debug/*.cpp
#                  OGF/basic/os/*.cpp
#                  OGF/basic/types/*.cpp
#                  OGF/math/common/*.cpp
#                  OGF/math/numeric/*.cpp
#                  OGF/math/symbolic/*.cpp)

add_library(graphite STATIC ${HEADERS} ${SOURCES})
