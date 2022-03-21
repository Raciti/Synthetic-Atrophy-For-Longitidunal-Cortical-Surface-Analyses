cmake_minimum_required(VERSION 2.8)

project(maskField)

find_package(ITK REQUIRED)
include(${ITK_USE_FILE})

find_package(VTK REQUIRED)
include(${VTK_USE_FILE})

add_executable(maskField gradientMask_surfacebased.cxx)

target_link_libraries(maskField ${ITK_LIBRARIES} ${VTK_LIBRARIES})
