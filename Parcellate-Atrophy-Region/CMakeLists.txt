cmake_minimum_required(VERSION 2.8)

project(parcellateAtrophyRegion)

find_package(ITK REQUIRED)
include(${ITK_USE_FILE})

find_package(VTK REQUIRED)
include(${VTK_USE_FILE})

add_executable(parcellateAtrophyRegion parcellateAtrophyRegion.cxx)

target_link_libraries(parcellateAtrophyRegion ${ITK_LIBRARIES} ${VTK_LIBRARIES})
