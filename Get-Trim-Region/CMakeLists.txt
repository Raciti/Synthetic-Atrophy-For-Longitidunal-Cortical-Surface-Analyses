cmake_minimum_required(VERSION 2.8)

project(getTrimRegion)

find_package(ITK REQUIRED)
include(${ITK_USE_FILE})

add_executable(getTrimRegion getTrimRegion.cxx)

target_link_libraries(getTrimRegion ${ITK_LIBRARIES})
