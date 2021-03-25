cmake_minimum_required(VERSION 3.5)

project(lab03 LANGUAGES CXX)

set(CMAKE_INCLUDE_CURRENT_DIR ON)

set(CMAKE_AUTOUIC ON)
set(CMAKE_AUTOMOC ON)
set(CMAKE_AUTORCC ON)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_compile_options(
	-Werror

	-Wall
	-Wextra
	-Wpedantic

	-Wcast-align
	-Wcast-qual
	-Wconversion
	-Wctor-dtor-privacy
	-Wenum-compare
	-Wfloat-equal
	-Wnon-virtual-dtor
	-Wold-style-cast
	-Woverloaded-virtual
	-Wredundant-decls
	-Wsign-conversion
	-Wsign-promo
)

set(
	HEADERS
	combinatorics.hpp
	experiment.hpp
	mainwindow.hpp
	random.hpp
	simulate.hpp
	statistics.hpp
)

set(
	SOURCES
	combinatorics.cpp
	experiment.cpp
	mainwindow.cpp
	random.cpp
	simulate.cpp
	statistics.cpp
)

set(
	FORMS
	mainwindow.ui
)

find_package(Qt5 COMPONENTS Widgets REQUIRED)
find_package(Qt5 COMPONENTS PrintSupport REQUIRED)

add_executable(${PROJECT_NAME} ${HEADERS} ${SOURCES} ${FORMS} main.cpp)

target_link_libraries(${PROJECT_NAME} PRIVATE Qt5::Widgets)
target_link_libraries(${PROJECT_NAME} PRIVATE Qt5::PrintSupport)
