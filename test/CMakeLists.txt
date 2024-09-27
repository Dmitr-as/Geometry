include(CTest)
find_package(Qt6 REQUIRED COMPONENTS Test)

FILE( GLOB TEST_FILES RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} *.cpp)

#set(CMAKE_CXX_FLAGS "-g -O0 -Wall -coverage -fprofile-arcs -ftest-coverage")

#set(OBJECT_DIR ${CMAKE_BINARY_DIR}/CMakeFiles/RunAdder.dir)
#message("-- Object files will be output to: ${OBJECT_DIR}")

# Create the gcov target. Run coverage tests with 'make gcov'
#add_custom_target(gcov
#    COMMAND mkdir -p coverage
#    COMMAND ${CMAKE_MAKE_PROGRAM} test
#    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
#    )
#add_custom_command(TARGET gcov
#    COMMAND echo "=================== GCOV ===================="
#    COMMAND gcov -b ${CMAKE_SOURCE_DIR}/*.cpp -o ${OBJECT_DIR}
#    | grep -A 5 "Adder.cpp" > CoverageSummary.tmp
#    COMMAND cat CoverageSummary.tmp
#    COMMAND echo "-- Coverage files have been output to ${CMAKE_BINARY_DIR}/coverage"
#    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/coverage
#    )
#add_dependencies(gcov RunAdder)
## Make sure to clean up the coverage folder
#set_property(DIRECTORY APPEND PROPERTY ADDITIONAL_MAKE_CLEAN_FILES coverage)

## Create the gcov-clean target. This cleans the build as well as generated
## .gcda and .gcno files.
#add_custom_target(scrub
#COMMAND ${CMAKE_MAKE_PROGRAM} clean
#COMMAND rm -f ${OBJECT_DIR}/*.gcno
#COMMAND rm -f ${OBJECT_DIR}/*.gcda
#WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
#)

set(CMAKE_CXX_FLAGS "-O0 -coverage -fprofile-arcs -ftest-coverage") # set compiler flags

# find required tools
find_program(LCOV lcov)
find_program(GENHTML genhtml)

if(LCOV AND GENHTML)
  # add coverage target
  add_custom_target(coverage
    # gather data
    COMMAND ${LCOV} --directory . --capture --output-file coverage.info
    # generate report
    COMMAND ${GENHTML} --demangle-cpp -o coverage coverage.info
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR})
# отчёт о покрытии
# lcov --directory . --capture --output-file coverage.info
# genhtml --demangle-cpp -o coverage coverage.info
endif()

foreach( fl ${TEST_FILES} )
    string( REPLACE ".cpp" "" file ${fl} )
    set(test_name "test_${file}")
    add_executable( ${test_name} ${fl} )

    target_link_libraries(${test_name}
        PRIVATE
            Qt::Test
    )

    add_test(NAME ${test_name} COMMAND ${test_name})

endforeach(fl)

