# CMake entry point
cmake_minimum_required (VERSION 2.6)

set( component_SOURCES  
	 Distance.h
	 Distance.cpp
	 Initialization.h
	 Initialization.cpp
	 IOHandler.h
	 IOHandler.cpp
	 Metric.h
	 Metric.cpp
	 PreComputing.h
	 PreComputing.cpp
	 Silhouette.h
	 Silhouette.cpp
   ) # Add the source-files for the component here

include(CheckCXXCompilerFlag)
CHECK_CXX_COMPILER_FLAG("-std=c++11" COMPILER_SUPPORTS_CXX11)
CHECK_CXX_COMPILER_FLAG("-std=c++0x" COMPILER_SUPPORTS_CXX0X)
if(COMPILER_SUPPORTS_CXX11)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
elseif(COMPILER_SUPPORTS_CXX0X)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x")
else()
        message(STATUS "The compiler ${CMAKE_CXX_COMPILER} has no C++11 support. Please use a different C++ compiler.")
endif()

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp")

include_directories(
    /usr/include/
    /usr/include/eigen3
    /usr/local/include/
    /usr/local/include/vtk-6.3/
	./src
    .
)

file( GLOB component_SOURCES *.cpp )

add_library( common ${component_SOURCES} )