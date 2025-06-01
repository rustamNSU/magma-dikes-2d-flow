# A quasi-2D model of dike propagation

## Building source code (Linux, macOS)

1. Clone this repository.

2. Install [vcpkg](https://github.com/microsoft/vcpkg) to manage third-party dependencies:
   - Follow the official guide to bootstrap vcpkg.
   - Install required packages:
     ```bash
     ./vcpkg install eigen3 suitesparse-umfpack highfive nlohmann-json spdlog fmt
     ```

3. Configure and build the project:
   - This script uses `CMAKE_TOOLCHAIN_FILE` to inform CMake of the vcpkg toolchain.
   - You may need to modify the value of `-D CMAKE_TOOLCHAIN_FILE=~/vcpkg/scripts/buildsystems/vcpkg.cmake` in the script to match the actual location of your `vcpkg` installation.
   - For Linux (tested on Fedora) and macOS:
     ```bash
     sh fedora-config.sh
     ```

## Citing this work

If you use this model in your research, please cite the following manuscript:

```bibtex
@unpublished{abdullin2025quasi2d,
  author       = {Abdullin, Rustam and Melnik, Oleg and Rust, Alison},
  title        = {A quasi-2D model of dike propagation},
  year         = {2025},
  note         = {Unpublished manuscript, in preparation}
}