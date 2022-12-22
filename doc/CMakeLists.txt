cmake_minimum_required(VERSION 3.7)
set(DOXYGEN_TARGET_NAME project_documentation) # Вторая часть имени файла с настройками Doxygen
get_filename_component(GEOCALCSOLUTION_DIR ${CMAKE_CURRENT_SOURCE_DIR}/../ ABSOLUTE) # Путь к корневой директории решения (solution)

set(DOXYGEN_INPUT_DIRECTORY # Откуда брать исходники для документации
    ${GEOCALCSOLUTION_DIR}/program/spml
    ${GEOCALCSOLUTION_DIR}/program/geocalc
)
set(DOXYGEN_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}) # Куда положить сгенерированную документацию

if(${CMAKE_VERSION} VERSION_GREATER "3.8")
#if(CMAKE_VERSION VERSION_LESS 3.9)
    # Версия CMake >= 3.9 умеет сама создавать файл настроек doxygen, что мы и делаем:
    message(WARNING "Doxygen do not need config file, options are taken from CMakeLists.txt")

    # Укажем ключи Doxygen
    set(DOXYGEN_DOXYFILE_ENCODING UTF-8)
    set(DOXYGEN_OUTPUT_LANGUAGE RUSSIAN)
    set(DOXYGEN_PROJECT_NAME "Геодезический калькулятор / Geodetic calculator")
    set(DOXYGEN_PROJECT_BRIEF "Программная документация / Program documentation")
    set(DOXYGEN_EXTRACT_ALL YES) #NO) #
    set(DOXYGEN_EXTRACT_PRIVATE YES) #NO) #
    set(DOXYGEN_EXTRACT_STATIC YES) #NO) #
    set(DOXYGEN_ENABLE_PREPROCESSING YES)
    set(DOXYGEN_MACRO_EXPANSION YES)
    set(DOXYGEN_EXPAND_ONLY_PREDEF YES)
#    set(DOXYGEN_EXTRACT_PRIVATE YES) # NO) #
    set(DOXYGEN_USE_MDFILE_AS_MAINPAGE "${GEOCALCSOLUTION_DIR}/README.md")
    set(DOXYGEN_GENERATE_HTML YES)
    set(DOXYGEN_DOT_GRAPH_MAX_NODES 4)
    set(DOXYGEN_FULL_PATH_NAMES NO)

    set(DOXYGEN_GENERATE_LATEX YES) # NO) #
    set(DOXYGEN_USE_PDFLATEX YES)
    set(DOXYGEN_PDF_HYPERLINKS YES) # NO) #
    set(DOXYGEN_LATEX_SOURCE_CODE YES) # NO) #
    set(DOXYGEN_LATEX_HIDE_INDICES)

    set(DOXYGEN_GENERATE_TESTLIST YES)
    set(DOXYGEN_GENERATE_TODOLIST YES)
#    set(DOXYGEN_CALL_GRAPH NO)
#    set(DOXYGEN_CALLER_GRAPH NO)
#    set(DOXYGEN_CLASS_DIAGRAMS NO)
#    set(DOXYGEN_CLASS_GRAPH NO)
#    set(DOXYGEN_COLLABORATION_GRAPH NO)
#    set(DOXYGEN_GROUP_GRAPHS NO)
#    set(DOXYGEN_INCLUDE_GRAPH NO)
#    set(DOXYGEN_INCLUDE_BY_GRAPH NO)
#    set(DOXYGEN_GRAPHICAL_HIERARCHY NO)
#    set(DOXYGEN_DIRECTORY_GRAPH NO)
#    set(DOXYGEN_HAVE_DOT NO)

    set(DOXYGEN_GENERATE_TREEVIEW YES)
#    set(DOXYGEN_SOURCE_BROWSER NO) # Не включать исходники! *.cpp
    set(DOXYGEN_SOURCE_BROWSER YES) # Включать исходники! *.cpp
    set(DOXYGEN_HTML_TIMESTAMP YES)

    # Создадим файл настроек Doxygen
    doxygen_add_docs(
        ${DOXYGEN_TARGET_NAME}
        ${GEOCALCSOLUTION_DIR}/README.md
        ${DOXYGEN_INPUT_DIRECTORY}
        ALL
        WORKING_DIRECTORY ${DOXYGEN_OUTPUT_DIRECTORY}
        COMMENT "Generating doxygen project documentation"
    )
else()
    # Версия CMake <= 3.9 НЕ умеет сама создавать файл настроек doxygen, поэтому файл должен быть заранее подготовлен!
    set(DOXYGEN_SETTINGS_FILE ${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile.${DOXYGEN_TARGET_NAME}) # CMake дает имя файлу настроек Doxyfile.{user_defined}

    if(EXISTS ${DOXYGEN_SETTINGS_FILE}) # Если файл настроек есть
        add_custom_target(
            ${DOXYGEN_TARGET_NAME}
            ALL
            COMMAND ${DOXYGEN_EXECUTABLE} ${DOXYGEN_SETTINGS_FILE} # doxygen-у надо передать имя файла с настройками
            WORKING_DIRECTORY ${DOXYGEN_OUTPUT_DIRECTORY}
            COMMENT "Generating doxygen project documentation"
        )
    else()
        message(WARNING "Doxygen config file was not found, project documenation creation skipped!")
    endif()
endif()
