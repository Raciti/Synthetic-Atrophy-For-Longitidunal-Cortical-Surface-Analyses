cmake_minimum_required(VERSION 2.8)
 
PROJECT(parcellateSingleSurface)
 
find_package(VTK REQUIRED)
include(${VTK_USE_FILE})

find_package(ITK REQUIRED)
include(${ITK_USE_FILE})

add_executable(parcellateSingleSurface parcellateSingleSurface.cxx)
 
target_link_libraries(parcellateSingleSurface ${VTK_LIBRARIES} ${ITK_LIBRARIES})
