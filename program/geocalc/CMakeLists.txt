cmake_minimum_required(VERSION 3.7)
project(geocalc LANGUAGES CXX) # geodetic calculator
get_filename_component(GEOCALCSOLUTION_DIR ${CMAKE_CURRENT_SOURCE_DIR}/../../ ABSOLUTE) # Путь к корневой директории решения (solution)
message(STATUS "CMake version: ${CMAKE_VERSION}, Project: ${PROJECT_NAME}, GEOCALCSOLUTION_DIR: ${GEOCALCSOLUTION_DIR}")

set(CMAKE_INCLUDE_CURRENT_DIR ON)
#SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pthread")

set(HEADERS
)

set(SOURCES
    src/main_geocalc.cpp
)

add_executable(${PROJECT_NAME} ${HEADERS} ${SOURCES})

target_include_directories(${PROJECT_NAME}
    PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}/include
        ${GEOCALCSOLUTION_DIR}/program/spml/include
)

target_link_libraries(${PROJECT_NAME}
    PUBLIC
        spml
)
#-----------------------------------------------------------------------------------------------------------------------
# Установим пакет geocalc
include(CPack)
install(TARGETS geocalc
        DESTINATION bin
        COMPONENT geocalc) # Установка самого пакета
#install(FILES "${GEOCALCSOLUTION_DIR}/desktop/coordcalc.desktop"
#        DESTINATION share/applications
#        COMPONENT coordcalc) # Установка ярлыка пакета
#install(FILES "${GEOCALCSOLUTION_DIR}/desktop/coordcalc.desktop"
#        DESTINATION share/applications
#        COMPONENT coordcalc) # Установка ярлыка пакета
#install(FILES "${GEOCALCSOLUTION_DIR}/icons/logo_coordcalc.gif"
#        DESTINATION share/icons
#        COMPONENT coordcalc) # Установка иконки пакета
cpack_add_component_group(geocalc)
cpack_add_component(geocalc
                    DISPLAY_NAME geocalc
                    DESCRIPTION "geocalc program"
                    GROUP geocalc
                    INSTALL_TYPES Full)
#-----------------------------------------------------------------------------------------------------------------------
# Включить замеры времени сборки
#set_property(GLOBAL PROPERTY RULE_LAUNCH_COMPILE "${CMAKE_COMMAND} -E time")
#set_property(GLOBAL PROPERTY RULE_LAUNCH_LINK "${CMAKE_COMMAND} -E time")
