#!/bin/bash
mkdir -p ./debug
cd ./debug

cmake -B .\
      -S ..\
      -D CMAKE_BUILD_TYPE=Debug\
      -D CMAKE_TOOLCHAIN_FILE=/Volumes/SamsungSSD/vcpkg/scripts/buildsystems/vcpkg.cmake
cmake --build . --parallel
