add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/utils)
add_executable(runner 
    main.cpp 
    Mesh.cpp
    Elasticity.cpp
    DikeData.cpp
    InputData.cpp
    Schedule.cpp
    DikeModel2d.cpp
    Utils.cpp
    MagmaState.cpp
    TimestepController.cpp
    ReservoirData.cpp
)
target_link_libraries(runner
    Eigen3::Eigen
    nlohmann_json::nlohmann_json
    HighFive
    spdlog::spdlog
    utils
)
target_include_directories(runner PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
install(TARGETS runner DESTINATION bin)