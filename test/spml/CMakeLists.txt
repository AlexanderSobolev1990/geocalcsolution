cmake_minimum_required(VERSION 3.7)
project(test_spml LANGUAGES CXX)
get_filename_component(GEOCALCSOLUTION_DIR ${CMAKE_CURRENT_SOURCE_DIR}/../../../ ABSOLUTE) # Путь к корневой директории решения (solution)
message(STATUS "CMake version: ${CMAKE_VERSION}, Project: ${PROJECT_NAME}, GEOCALCSOLUTION_DIR: ${GEOCALCSOLUTION_DIR}")

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/bin/${CMAKE_BUILD_TYPE}/tests) # Директрия для тестов

enable_testing(true)

set(CMAKE_INCLUDE_CURRENT_DIR ON)
#set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
#set(BOOST_TEST_LOG_LEVEL all)

find_package (Boost COMPONENTS unit_test_framework REQUIRED)
#find_package(Armadillo REQUIRED)
#include_directories(/usr/include/python3.8)

set(CMAKE_INCLUDE_CURRENT_DIR ON)
#set(CTEST_OUTPUT_ON_FAILURE ON)
#-----------------------------------------------------------------------------------------------------------------------
# geodesy
add_executable(test_spml_geodesy test_spml_geodesy.cpp)
add_test(NAME test_spml_geodesy COMMAND test_spml_geodesy)
target_link_libraries(test_spml_geodesy spml ${Boost_LIBRARIES})
#-----------------------------------------------------------------------------------------------------------------------

