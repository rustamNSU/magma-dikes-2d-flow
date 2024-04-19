#!/bin/bash
mkdir -p ./build
cd ./build

cmake -B .\
      -S ..\
      -D CMAKE_BUILD_TYPE=Debug\
      -D CMAKE_TOOLCHAIN_FILE=~/vcpkg/scripts/buildsystems/vcpkg.cmake
cmake --build . --parallel
