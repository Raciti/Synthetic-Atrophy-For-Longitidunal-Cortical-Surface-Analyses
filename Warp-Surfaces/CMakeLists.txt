cmake_minimum_required(VERSION 2.8)

PROJECT(warpSurfs)

Find_package(VTK REQUIRED)
include(${VTK_USE_FILE})

Find_package(ITK REQUIRED)
include(${ITK_USE_FILE})

add_executable(warpSurfs warpSurfaces.cxx)

target_link_libraries(warpSurfs ${ITK_LIBRARIES} ${VTK_LIBRARIES})
