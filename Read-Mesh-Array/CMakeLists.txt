cmake_minimum_required(VERSION 3.12)

PROJECT(readMeshArray)

Find_package(VTK REQUIRED)
include(${VTK_USE_FILE})

add_executable(readMeshArray readMeshArray.cxx)

target_link_libraries(readMeshArray ${VTK_LIBRARIES})
