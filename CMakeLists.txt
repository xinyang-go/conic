cmake_minimum_required(VERSION 3.16)
project(conic)

set(CMAKE_CXX_STANDARD 17)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Werror=return-type")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -Wall")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3 -flto -march=native")

find_package(Eigen3 3.3 REQUIRED)

add_library(conic src/conic.cpp)
target_include_directories(conic PUBLIC include ${EIGEN3_INCLUDE_DIR})

option(BUILD_EXAMPLE "build example" ON)
if(BUILD_EXAMPLE)
    find_package(OpenCV 4 REQUIRED)
    find_package(Ceres 2 REQUIRED)

    add_executable(example_static example/example_static.cpp)
    target_link_libraries(example_static conic ${OpenCV_LIBS})
    set_target_properties(example_static PROPERTIES
            RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/example")
            
    add_executable(example_dynamic example/example_dynamic.cpp)
    target_link_libraries(example_dynamic conic ${OpenCV_LIBS} ${CERES_LIBRARIES})
    set_target_properties(example_dynamic PROPERTIES
            RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/example")
endif()
