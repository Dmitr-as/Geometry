include(CTest)

FILE( GLOB TEST_FILES RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} *.cpp)

SET(GCC_COVERAGE_COMPILE_FLAGS "-fprofile-arcs -ftest-coverage")
set(CMAKE_CXX_FLAGS "-O0 -coverage -fprofile-arcs -ftest-coverage") # set compiler flags
SET(GCC_COVERAGE_LINK_FLAGS "-lgcov")

# find required tools
find_program(GCOV lcov)
find_program(LCOV lcov)
find_program(GENHTML genhtml)

if(LCOV AND GENHTML)
    # add coverage target
    add_custom_target(coverage
        COMMAND ctest
        COMMAND ${LCOV} --directory . --capture --output-file coverage.info
        COMMAND ${GENHTML} --demangle-cpp -o coverage coverage.info # generate report
        WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
    )
# отчёт о покрытии в консоли
# lcov --directory . --capture --output-file coverage.info
# genhtml --demangle-cpp -o coverage coverage.info
endif()

## Create the gcov-clean target. This cleans the build as well as generated
## .gcda and .gcno files.
add_custom_target(covclear
    COMMAND ${CMAKE_MAKE_PROGRAM} clean
    COMMAND rm -f ${OBJECT_DIR}/*.gcno
    COMMAND rm -f ${OBJECT_DIR}/*.gcda
    COMMAND rm -R coverage
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
)

foreach( fl ${TEST_FILES} )
    string( REPLACE ".cpp" "" file ${fl} )
    set(test_name "test_${file}")
    add_executable( ${test_name} ${fl} )
    add_test(NAME ${test_name} COMMAND ${test_name})
endforeach(fl)

