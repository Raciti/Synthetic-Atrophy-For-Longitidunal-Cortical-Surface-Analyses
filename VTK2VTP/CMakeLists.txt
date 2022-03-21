cmake_minimum_required(VERSION 2.8)

PROJECT(VTK2VTP)

Find_package(VTK REQUIRED)
include(${VTK_USE_FILE})

add_executable(VTK2VTP vtk2vtp.cxx)

target_link_libraries(VTK2VTP ${VTK_LIBRARIES})
