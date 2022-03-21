cmake_minimum_required(VERSION 2.8)

project(performErosion)

find_package(ITK REQUIRED)
include(${ITK_USE_FILE})

add_executable(performErosion performErosion.cxx)

target_link_libraries(performErosion ${ITK_LIBRARIES})
