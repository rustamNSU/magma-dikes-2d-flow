#!/bin/bash
mkdir -p ./release
cd ./release

cmake -D'CMAKE_BUILD_TYPE=Release'\
      -D'Eigen3_DIR="/opt/homebrew/include/eigen3"' ..
cmake --build . --parallel
