cmake_minimum_required(VERSION 2.8)

PROJECT(copyScalars)

Find_package(VTK REQUIRED)
include(${VTK_USE_FILE})

add_executable(copyScalars copyScalars.cxx)

target_link_libraries(copyScalars ${VTK_LIBRARIES})
