include_directories( ${Boost_INCLUDE_DIR}
                     ${PROJECT_SOURCE_DIR}/include
                     ${PROJECT_SOURCE_DIR}/include/hj_3rd)

set ( HEADERS Parameterization.h
              Barycentric.h
              ParamPatch.h
              ChartCreator.h
              TransFunctor.h
              ParamDrawer.h
              TriDistortion.h
              Parameter.h
              )

set ( SOURCES Parameterization.cc
              Barycentric.cc
              TransFunctor.cc
              ChartCreator.cc
              ParamDrawer.cc
              TriDistortion.cc
              Parameter.cc
              )
              

                  

add_library(param STATIC ${HEADERS} ${SOURCES})
