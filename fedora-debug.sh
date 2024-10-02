#!/bin/bash
mkdir -p ./debug
cd ./debug

cmake -B .\
      -S ..\
      -D CMAKE_BUILD_TYPE=Debug \
      -D CMAKE_TOOLCHAIN_FILE=~/vcpkg/scripts/buildsystems/vcpkg.cmake \
      -D BUILD_TESTS=ON
cmake --build . --parallel -j4

# -DCMAKE_BUILD_TYPE=Debug \
# -DCMAKE_C_FLAGS_DEBUG="-g -O0" \
# -DCMAKE_CXX_FLAGS_DEBUG="-g -O0" \