cmake_minimum_required(VERSION 3.6)
project(nurbs_library)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -lboost_iostreams -lboost_system -lboost_filesystem")
set(SOURCE_FILES main.cpp NurbsCurve.cpp NurbsCurve.h gnuplot-iostream.h)
add_executable(nurbs_library ${SOURCE_FILES})