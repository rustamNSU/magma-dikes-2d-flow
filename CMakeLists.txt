cmake_minimum_required (VERSION 3.12)
project("dike-2d-flow")

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

set(CMAKE_INSTALL_PREFIX ${CMAKE_CURRENT_SOURCE_DIR}/Installation)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

find_package (Eigen3 3.4 REQUIRED NO_MODULE)
if(Eigen3_FOUND)
    message("Eigen3 found: ${Eigen3_DIR}")
else()
	message(FATAL_ERROR "Eigen3 was not found!")
endif()

# Download and build external projects
include(FetchContent)
include("${PROJECT_SOURCE_DIR}/cmake/fetchHighFive.cmake")
add_subdirectory("${PROJECT_SOURCE_DIR}/3rdparty")
add_subdirectory("${PROJECT_SOURCE_DIR}/src")