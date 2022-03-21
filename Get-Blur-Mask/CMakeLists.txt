cmake_minimum_required(VERSION 2.8)

project(getBlurMask)

find_package(ITK REQUIRED)
include(${ITK_USE_FILE})

find_package(VTK REQUIRED)
include(${VTK_USE_FILE})

add_executable(getBlurMask getBlurMask.cxx)

target_link_libraries(getBlurMask ${ITK_LIBRARIES} ${VTK_LIBRARIES})
