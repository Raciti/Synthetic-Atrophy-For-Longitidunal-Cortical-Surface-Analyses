cmake_minimum_required(VERSION 2.8)

project(distanceBetweenSurfaces)

find_package(ITK REQUIRED)
include(${ITK_USE_FILE})

find_package(VTK REQUIRED)
include(${VTK_USE_FILE})

add_executable(distanceBetweenSurfaces distanceBetweenSurfaces.cxx)

target_link_libraries(distanceBetweenSurfaces ${ITK_LIBRARIES} ${VTK_LIBRARIES})
