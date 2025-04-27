cmake_minimum_required(VERSION 3.7)

project(geometry LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

FILE( GLOB GET_HEAD_FILES *.h )

# add_executable(tests_base_geometry
#     ${GET_HEAD_FILES}
#     tests.cpp
# )

enable_testing()

add_subdirectory( test )
