enable_language(Fortran)                                    
SET(CMAKE_Fortran_FLAGS "-fno-underscoring")

include_directories( ${Boost_INCLUDE_DIR}
                     ${PROJECT_SOURCE_DIR}/include
                     ${PROJECT_SOURCE_DIR}/include/hj_3rd
                     /usr/include/libmpalglib)

set ( HEADERS Parameterization.h
              Barycentric.h
              ParamPatch.h
              ChartCreator.h
              TransFunctor.h
              ParamDrawer.h
              TriDistortion.h
              Parameter.h
              ChartOptimization.h
              quad_dof.h
              )

set ( SOURCES Parameterization.cc
              Barycentric.cc
              TransFunctor.cc
              ChartCreator.cc
              ParamDrawer.cc
              TriDistortion.cc
              Parameter.cc
              ChartOptimization.cc
              )              

add_library(quad_dof SHARED  quad_dof-1.f)
            
add_library(param STATIC ${HEADERS} ${SOURCES} )
