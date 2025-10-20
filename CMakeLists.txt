cmake_minimum_required(VERSION 3.0)

project(SeaMeshProject)

set(CMAKE_CXX_STANDARD 14)

find_package(VTK REQUIRED)
find_package(Python COMPONENTS Interpreter Development)
find_package(Boost 1.83.0 REQUIRED)
find_package(Eigen3 3.4 REQUIRED NO_MODULE)
find_package(GDAL REQUIRED)

message (STATUS "VTK_VERSION    : ${VTK_VERSION}")
message (STATUS "Python_VERSION : ${Python_VERSION}")
message (STATUS "Boost_VERSION  : ${Boost_VERSION}")
message (STATUS "Eigen3_VERSION : ${Eigen3_VERSION}")

add_subdirectory(src)
