cmake_minimum_required(VERSION 3.7)
project(spml LANGUAGES CXX) # SPML (Special Program Modules Library) - СБ ПМ (Cпециальная Библиотека Программных Модулей)
get_filename_component(GEOCALCSOLUTION_DIR ${CMAKE_CURRENT_SOURCE_DIR}/../../ ABSOLUTE) # Путь к корневой директории решения (solution)
message(STATUS "CMake version: ${CMAKE_VERSION}, Project: ${PROJECT_NAME}, GEOCALCSOLUTION_DIR: ${GEOCALCSOLUTION_DIR}")

set(CMAKE_INCLUDE_CURRENT_DIR ON)
#set(CMAKE_CXX_STANDARD 11)
#set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(HEADERS
    include/spml.h
    include/consts.h
    include/convert.h
    include/compare.h    
    include/geodesy.h
    include/units.h
    )

set(SOURCES
    src/spml.cpp
    src/convert.cpp
    src/geodesy.cpp
    )

add_library(${PROJECT_NAME} STATIC ${HEADERS} ${SOURCES}) # Статическая библиотека

target_include_directories(${PROJECT_NAME}
    PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}/include        
)

find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)
find_package(Armadillo REQUIRED)
find_package(Boost REQUIRED COMPONENTS system program_options)

include_directories(${ARMADILLO_INCLUDE_DIRS})
target_link_libraries(${PROJECT_NAME}
    PUBLIC
        ${ARMADILLO_LIBRARIES}
        pthread
        Boost::system
        Boost::program_options
)

## Включить замеры времени сборки
#set_property(GLOBAL PROPERTY RULE_LAUNCH_COMPILE "${CMAKE_COMMAND} -E time")
#set_property(GLOBAL PROPERTY RULE_LAUNCH_LINK "${CMAKE_COMMAND} -E time")
