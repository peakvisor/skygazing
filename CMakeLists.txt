cmake_minimum_required(VERSION 3.1)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -O3 -ffast-math -Wall  -Wno-missing-braces")
set(CMAKE_CXX_STANDARD 17)

project(skygazing_tests)
add_executable(skygazing_tests test/skygazing_test.cpp)
target_include_directories(skygazing_tests PUBLIC lib include)

if (SKYGAZING_AAPLUS_SRC)
    project(skygazing_tests_aaplus)
    add_executable(skygazing_tests_aaplus test/skygazing_test_aaplus.cpp)
    target_include_directories(skygazing_tests_aaplus PUBLIC ${SKYGAZING_AAPLUS_SRC} lib include)
    target_link_libraries(skygazing_tests_aaplus ${SKYGAZING_AAPLUS_SRC}/cmake-build-release/lib/libaaplus.dylib)
endif()
