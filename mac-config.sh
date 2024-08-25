#!/bin/bash
mkdir -p ./build
cd ./build

cmake -B .\
      -S ..\
      -D CMAKE_BUILD_TYPE=Release\
      -D CMAKE_TOOLCHAIN_FILE=/Volumes/SamsungSSD/vcpkg/scripts/buildsystems/vcpkg.cmake
cmake --build . --parallel
