cmake_minimum_required (VERSION 3.12)
project("dike-2d-flow")

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

set(CMAKE_INSTALL_PREFIX ${CMAKE_CURRENT_SOURCE_DIR}/Installation)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

find_package(Eigen3 3.4 REQUIRED)
# find_package(suitesparse CONFIG REQUIRED)
find_package(UMFPACK CONFIG REQUIRED)
find_package(HighFive REQUIRED)
find_package(nlohmann_json REQUIRED)
find_package(spdlog REQUIRED)
find_package(fmt REQUIRED)
option(USE_OMP "with openmp" OFF)
option(BUILD_TESTS "build tests" OFF)


if(USE_OMP)
	find_package(OpenMP REQUIRED)
	set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    add_compile_definitions(USE_OMP)
endif()
set(HIGHFIVE_USE_BOOST Off)
set(HIGHFIVE_USE_EIGEN On)


add_subdirectory(${PROJECT_SOURCE_DIR}/src)
# if(BUILD_TESTS)
	add_subdirectory(${PROJECT_SOURCE_DIR}/tests)
# endif()