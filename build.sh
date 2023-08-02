#!/bin/bash
mkdir -p ./build
cd ./build

cmake -D'CMAKE_BUILD_TYPE=Debug'\
      -D'Eigen3_DIR="/opt/homebrew/include/eigen3"' ..
cmake --build . --parallel
