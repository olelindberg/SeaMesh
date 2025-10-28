cmake_minimum_required(VERSION 3.0)

project(SeaMeshProject LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(PROJ_DIR "/usr/lib/x86_64-linux-gnu") # root install folder
set(CMAKE_PREFIX_PATH "${PROJ_DIR}")

find_package(VTK REQUIRED)
find_package(Python COMPONENTS Interpreter Development)
find_package(Boost 1.83.0 REQUIRED)
find_package(Eigen3 3.4 REQUIRED NO_MODULE)
find_package(GDAL REQUIRED)
find_package(PROJ REQUIRED CONFIG)

message (STATUS "VTK_VERSION    : ${VTK_VERSION}")
message (STATUS "Python_VERSION : ${Python_VERSION}")
message (STATUS "Boost_VERSION  : ${Boost_VERSION}")
message (STATUS "Eigen3_VERSION : ${Eigen3_VERSION}")
message (STATUS "PROJ_VERSION   : ${PROJ_VERSION}")

add_subdirectory(src)
