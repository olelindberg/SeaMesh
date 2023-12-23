set(SOURCE_FILES 
main.cpp
)

add_executable(test01_boost_rtree ${SOURCE_FILES})

target_link_libraries(test01_boost_rtree PUBLIC
)

target_include_directories(test01_boost_rtree PUBLIC
)
