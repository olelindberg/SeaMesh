add_subdirectory(test)

list (APPEND CMAKE_MODULE_PATH "${ShipMesh_SOURCE_DIR}/cmake-modules")

set(SOURCE_FILES 
main.cpp
#GeoTiffReader.cpp
VtkUtil.cpp
QuadTreeBuilder.cpp
QuadTreeSearch.cpp
QuadTree.cpp
Profiler.cpp
QuadTreePointsInserter.cpp
QuadTreeSeaBuilder.cpp
QuadTreeBuilderBase.cpp
TreesUtil.cpp
MortonUtil.cpp
)

add_executable(SeaMesh ${SOURCE_FILES})

target_link_libraries(SeaMesh PUBLIC
${VTK_LIBRARIES}
/usr/local/shapelib-1.5.0/lib/libshp.so
Eigen3::Eigen
GDAL::GDAL
)
target_include_directories(SeaMesh PUBLIC
/usr/local/shapelib-1.5.0/include  
${GDAL_INCLUDE_DIR}
${Boost_INCLUDE_DIR}
)
