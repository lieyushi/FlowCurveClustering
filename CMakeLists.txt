# CMake entry point
cmake_minimum_required (VERSION 2.6)
project (MyClustering)

add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/src/Common)

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

# add "-O2" for library eigen optimization for matrix computation.
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp -O2")

include_directories(
    /usr/include/
    /usr/include/eigen3
    /usr/local/include/
    /usr/local/include/vtk-6.3/
    src/Common/
	  ./src
    .
)

set(ALL_LIBS
	common
)

add_definitions(
#	-DTW_STATIC
#	-DTW_NO_LIB_PRAGMA
#	-DTW_NO_DIRECT3D
#	-DGLEW_STATIC
#	-D_CRT_SECURE_NO_WARNINGS
)

###############################################
# generate kmeans executable
set(KMeans_SRC
	${CMAKE_CURRENT_SOURCE_DIR}/src/KMeans/PCA_Cluster.h
    ${CMAKE_CURRENT_SOURCE_DIR}/src/KMeans/PCA_Cluster.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/src/KMeans/main.cpp
)
add_executable(kmeans ${KMeans_SRC})
target_link_libraries(kmeans ${ALL_LIBS})

###############################################


###############################################
# generate query executable
set(Query_SRC
	${CMAKE_CURRENT_SOURCE_DIR}/src/Query/Query.h
    ${CMAKE_CURRENT_SOURCE_DIR}/src/Query/Query.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/src/Query/main.cpp
)
add_executable(query ${Query_SRC})
target_link_libraries(query ${ALL_LIBS})
###############################################


###############################################
# generate birch executable
set( Birch_SRC
   )
file( GLOB Birch_SRC ${CMAKE_CURRENT_SOURCE_DIR}/src/Birch/*.cpp )  
add_executable(birch ${Birch_SRC})
target_link_libraries(birch ${ALL_LIBS})
###############################################


###############################################
# generate birch executable
set( KMediods_SRC
   )
file( GLOB KMediods_SRC ${CMAKE_CURRENT_SOURCE_DIR}/src/KMedoids/*.cpp )  
add_executable(kmedoids ${KMediods_SRC})
target_link_libraries(kmedoids ${ALL_LIBS})
###############################################


###############################################
# generate dbscan executable
set( DBSCAN_SRC
   )
file( GLOB DBSCAN_SRC ${CMAKE_CURRENT_SOURCE_DIR}/src/DBSCAN/*.cpp )  
add_executable(dbscan ${DBSCAN_SRC})
target_link_libraries(dbscan ${ALL_LIBS})
###############################################


###############################################
# generate dbscan executable
set( OPTICS_SRC
   )
file( GLOB OPTICS_SRC ${CMAKE_CURRENT_SOURCE_DIR}/src/OPTICS/*.cpp )  
add_executable(optics ${OPTICS_SRC})
target_link_libraries(optics ${ALL_LIBS})
###############################################


###############################################
# generate ahc_dist executable
set( AHC_DIST_SRC
   )
file( GLOB AHC_DIST_SRC ${CMAKE_CURRENT_SOURCE_DIR}/src/AHC_DIST/*.cpp )  
add_executable(ahc_dist ${AHC_DIST_SRC})
target_link_libraries(ahc_dist ${ALL_LIBS})
###############################################


###############################################
# generate ahc executable
set( AHC_SRC
   )
file( GLOB AHC_SRC ${CMAKE_CURRENT_SOURCE_DIR}/src/AHC/*.cpp )  
add_executable(ahc ${AHC_SRC})
target_link_libraries(ahc ${ALL_LIBS})
###############################################


###############################################
# generate ahc executable
set( SC_SRC
   )
file( GLOB SC_SRC ${CMAKE_CURRENT_SOURCE_DIR}/src/SpectralClustering/*.cpp )  
add_executable(sc ${SC_SRC})
target_link_libraries(sc ${ALL_LIBS})
###############################################


###############################################
# generate affinityPropagation executable
set( AP_SRC
   )
file( GLOB AP_SRC ${CMAKE_CURRENT_SOURCE_DIR}/src/AffinityPropagation/*.cpp )  
add_executable(affinityPropagation ${AP_SRC})
target_link_libraries(affinityPropagation ${ALL_LIBS})
###############################################


###############################################
# generate readClustering executable
set( RC_SRC
   )
file( GLOB RC_SRC ${CMAKE_CURRENT_SOURCE_DIR}/src/ReadClustering/*.cpp )  
add_executable(readClustering ${RC_SRC})
target_link_libraries(readClustering ${ALL_LIBS})
###############################################