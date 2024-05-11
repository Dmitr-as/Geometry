cmake_minimum_required(VERSION 3.7)

project(geometry LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

enable_testing()

find_package(Boost REQUIRED)
find_package(Boost COMPONENTS unit_test_framework REQUIRED)

FILE( GLOB GET_HEAD_FILES *.h )

add_executable(tests_base_geometry
    ${GET_HEAD_FILES}
    tests_base_geometry.cpp
)

target_link_libraries(tests_base_geometry
    ${Boost_LIBRARIES}
)

add_subdirectory( test )
