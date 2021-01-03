cmake_minimum_required(VERSION 3.11)

if(${CMAKE_SOURCE_DIR} STREQUAL ${CMAKE_BINARY_DIR})
    message(FATAL_ERROR "DO NOT BUILD in-tree.")
endif()


project(lost)

file(GLOB_RECURSE SOURCES CONFIGURE_DEPENDS "src/*.cpp")

add_executable(${CMAKE_PROJECT_NAME} ${SOURCES})

set_property(TARGET ${CMAKE_PROJECT_NAME} PROPERTY CXX_STANDARD 11)
set_property(TARGET ${CMAKE_PROJECT_NAME} PROPERTY CXX_STANDARD_REQUIRED ON)

target_link_libraries(${CMAKE_PROJECT_NAME} cairo)

if (UNIX)
    execute_process(
            COMMAND ./download-bsc.sh
            OUTPUT ./bright-star-catalog.tsv
            WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
    )
else()

endif()

enable_testing()
add_subdirectory(test)
