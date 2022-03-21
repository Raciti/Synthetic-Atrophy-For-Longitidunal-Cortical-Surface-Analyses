cmake_minimum_required(VERSION 2.8)

PROJECT(VTP2VTK)

Find_package(VTK REQUIRED)
include(${VTK_USE_FILE})

add_executable(VTP2VTK vtp2vtk.cxx)

target_link_libraries(VTP2VTK ${VTK_LIBRARIES})
