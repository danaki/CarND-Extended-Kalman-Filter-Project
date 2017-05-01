# Sorry, man, builtin Ubuntu in Windows 10 has cmake 2.8.12.2
cmake_minimum_required (VERSION 2.8)

set(PROJECT_NAME_STR ExtendedKF)
project(${PROJECT_NAME_STR} C CXX)

find_package(Threads REQUIRED)

set(CMAKE_VERBOSE_MAKEFILE ON)
set(CMAKE_BUILD_STRIP FALSE)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g2")
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -g2")

if(CMAKE_COMPILER_IS_GNUCXX)
    add_definitions(-Wall -ansi -Wno-deprecated -pthread -g)
endif()

add_definitions(-std=c++0x)

set(sources
    src/FusionEKF.cpp
    src/kalman_filter.cpp
    src/main.cpp
    src/tools.cpp)

set(COMMON_INCLUDES src)

add_executable(ExtendedKF ${sources})

set(EXT_PROJECTS_DIR ${PROJECT_SOURCE_DIR}/ext)


#-------------------
# Module source
#-------------------
set(PROJECT_LIB_NAME ${PROJECT_NAME_STR}_impl)

include_directories(${COMMON_INCLUDES})
file(GLOB SRC_FILES ${PROJECT_SOURCE_DIR}/src/*.cpp)
add_library(${PROJECT_LIB_NAME} ${SRC_FILES})

add_subdirectory(${EXT_PROJECTS_DIR}/gtest)


#-------------------
# Test
#-------------------
add_definitions(${MSVC_COMPILER_DEFS})
enable_testing()
set(PROJECT_TEST_NAME ${PROJECT_NAME_STR}_test)
include_directories(${GTEST_INCLUDE_DIRS} ${COMMON_INCLUDES})

file(GLOB TEST_SRC_FILES ${PROJECT_SOURCE_DIR}/test/*.cpp)
add_executable(${PROJECT_TEST_NAME} ${TEST_SRC_FILES})
add_dependencies(${PROJECT_TEST_NAME} googletest)


if(NOT WIN32 OR MINGW)
    target_link_libraries(${PROJECT_TEST_NAME}
        ${GTEST_LIBS_DIR}/libgtest.a
        ${GTEST_LIBS_DIR}/libgtest_main.a
        )
else()
    target_link_libraries(${PROJECT_TEST_NAME}
        debug ${GTEST_LIBS_DIR}/DebugLibs/${CMAKE_FIND_LIBRARY_PREFIXES}gtest${CMAKE_FIND_LIBRARY_SUFFIXES}
        optimized ${GTEST_LIBS_DIR}/ReleaseLibs/${CMAKE_FIND_LIBRARY_PREFIXES}gtest${CMAKE_FIND_LIBRARY_SUFFIXES}
        )
    target_link_libraries(${PROJECT_TEST_NAME}
        debug ${GTEST_LIBS_DIR}/DebugLibs/${CMAKE_FIND_LIBRARY_PREFIXES}gtest_main${CMAKE_FIND_LIBRARY_SUFFIXES}
        optimized ${GTEST_LIBS_DIR}/ReleaseLibs/${CMAKE_FIND_LIBRARY_PREFIXES}gtest_main${CMAKE_FIND_LIBRARY_SUFFIXES}
        )
endif()

target_link_libraries(${PROJECT_TEST_NAME} ${CMAKE_THREAD_LIBS_INIT} ${PROJECT_LIB_NAME})

add_test(test1 ${PROJECT_TEST_NAME})
