#[[
add_executable(b_base_geometry
    btest.cpp
)

target_link_libraries(b_base_geometry
    ${Boost_LIBRARIES}
)

add_test(NAME Test_base_geometry COMMAND b_base_geometry)



add_executable(c_base_geometry
    ctest.cpp
)

target_link_libraries(c_base_geometry
    ${Boost_LIBRARIES}
)

add_test(NAME Test_cbase_geometry COMMAND c_base_geometry)
]]#

FILE( GLOB TEST_FILES RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} *.cpp)

foreach( fl ${TEST_FILES} )
    string( REPLACE ".cpp" "" file ${fl} )
    set(test_name "test_${file}")
    add_executable( ${test_name} ${fl} )

    target_link_libraries( ${test_name} PRIVATE ${Boost_LIBRARIES} )
    #[[
    target_include_directories( ${test_name}
        PRIVATE
            "$<BUILD_INTERFACE:${CMAKE_BINARY_DIR}/src>"
            "$<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>"
    )
    ]]#
    add_test(NAME ${test_name} COMMAND ${test_name})

endforeach(fl)
