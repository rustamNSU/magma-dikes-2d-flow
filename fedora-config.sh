#!/bin/bash
mkdir -p ./build
cd ./build

cmake -B .\
      -S ..\
      -D CMAKE_BUILD_TYPE=RelWithDebInfo\
      -D CMAKE_TOOLCHAIN_FILE=~/vcpkg/scripts/buildsystems/vcpkg.cmake
cmake --build . --parallel
