set(HIGHFIVE_USE_BOOST Off)
set(HIGHFIVE_USE_EIGEN On)

FetchContent_Declare(
    highfive
    URL https://github.com/BlueBrain/HighFive/archive/refs/tags/v2.7.1.tar.gz)
FetchContent_MakeAvailable(highfive)