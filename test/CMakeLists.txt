FILE( GLOB TEST_FILES RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} *.cpp)

foreach( fl ${TEST_FILES} )
    string( REPLACE ".cpp" "" file ${fl} )
    set(test_name "test_${file}")
    add_executable( ${test_name} ${fl} )

    add_test(NAME ${test_name} COMMAND ${test_name})

endforeach(fl)

