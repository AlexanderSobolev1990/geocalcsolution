cmake_minimum_required(VERSION 3.7)
project(geocalcsolution) # GEOCALCSOLUTION (Geodetic Calculator Solution) - Геодезический калькулятор
message(STATUS "CMake version: ${CMAKE_VERSION}, Project: ${PROJECT_NAME}")

if(NOT CMAKE_BUILD_TYPE) # Если не указан тип сборки, собрать в режиме Release
  set(CMAKE_BUILD_TYPE Release)
endif()
message(WARNING "Build mode: ${CMAKE_BUILD_TYPE}")

# Библиотеки - отдельно в директорию lib
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/bin/${CMAKE_BUILD_TYPE}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/bin/${CMAKE_BUILD_TYPE}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/bin/${CMAKE_BUILD_TYPE})

if(CMAKE_BUILD_TYPE MATCHES "^[Rr]elease") # Актуально для Linux
    option(BUILD_DOCUMENTATION "Build project documentation (Requires Doxygen)" ON) # Опционально построение документации: ON|OFF
endif()

set(CMAKE_CXX_FLAGS_RELEASE "-D__DEBUG__")

add_compile_options(
#    -Wall
#    -fprofile-arcs
#    -ftest-coverage
    "$<$<CONFIG:RELEASE>:-O3>"
    "$<$<CONFIG:DEBUG>:-O0>"
)

#add_link_options(
#    -lgcov
#    --coverage
#    )

add_subdirectory(program/spml) # SPML (Special Program Modules Library) - СБ ПМ (Cпециальная Библиотека Программных Модулей)
add_subdirectory(program/geocalc) # GEOCALC (Geodetic Calculator) - Геодезический калькулятор

# Тесты
#enable_testing()
#add_subdirectory(test/spml)
#add_subdirectory(test/spml/geodesy)

# Генерация документации
if(BUILD_DOCUMENTATION)
    find_package(Doxygen
        REQUIRED dot) # dot - требуется Graphviz dot утилита для построения графиков в документации
    if(DOXYGEN_FOUND)
        add_subdirectory(doc) # В поддиректории проекта /doc должен быть CMakeLists.txt с настройками Doxygen
    else()
        message(WARNING "Doxygen need to be installed to generate the doxygen documentation")
    endif()
endif()

# Настройки создания пакета при помощи CPack
set(CPACK_GENERATOR "DEB")
set(CPACK_DEB_COMPONENT_INSTALL ON)

#set(CPACK_DEBIAN_GEOCALC_PACKAGE_DEPENDS "libarmadillo-dev,liblapack-dev,libblas-dev,libboost-dev,libboost-test-dev")
set(CPACK_COMPONENT_GEOCALC_DESCRIPTION "Geodetic Calculator")

set(CPACK_PACKAGE_VERSION_MAJOR 22) # Год
set(CPACK_PACKAGE_VERSION_MINOR 12) # Месяц
set(CPACK_PACKAGE_VERSION_PATCH 21) # День-Версия

set(CPACK_SET_DESTDIR true) # Если выставлено - будет инсталлироваться в /usr/local/bin, иначе в /usr/bin
#set(CMAKE_INSTALL_PREFIX /opt)
#set(CPACK_PACKAGING_INSTALL_PREFIX /opt)

include(CPack)
