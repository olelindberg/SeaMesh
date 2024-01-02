add_executable(test02_boost_polygon_voronoi main.cpp)

target_link_libraries(test02_boost_polygon_voronoi PUBLIC
${VTK_LIBRARIES}
Eigen3::Eigen
)

target_include_directories(test02_boost_polygon_voronoi PUBLIC
${Boost_INCLUDE_DIR}
)
