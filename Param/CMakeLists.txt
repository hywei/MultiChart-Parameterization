cmake_minimum_required(VERSION 2.8)

if(WIN32)
        set(CMAKE_FIND_LIBRARY_PREFIXES "")
        set(CMAKE_FIND_LIBRARY_SUFFIXES ".lib")
else()
        set(CMAKE_FIND_LIBRARY_PREFIXES "lib")
        set(CMAKE_FIND_LIBRARY_SUFFIXES ".so")
endif()

#find_package( OpenGL REQUIRED)
#find_package( Boost REQUIRED)

set( Boost_INCLUDE_DIR "/usr/include/")

#set(CMAKE_BUILD_TYPE Debug)
#set(CMAKE_CXX_COMPILER /usr/bin/g++)


project(Parameter)


message("${BOOST_INCLUDE_DIR}")

message("${CMAKE_BUILD_TYPE}")
message("${CMAKE_CXX_COMPILER}")

if(WIN32)
        set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/bin/win32)
else()
        set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/bin/linux)
endif()
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/lib)

add_subdirectory(src/Graphite)
add_subdirectory(src/Common)
add_subdirectory(src/Numerical)
add_subdirectory(src/ModelMesh)
add_subdirectory(src/MainWindow)
add_subdirectory(src/OpenGL)
add_subdirectory(src/Param)
add_subdirectory(src/UI)

     
        
       
