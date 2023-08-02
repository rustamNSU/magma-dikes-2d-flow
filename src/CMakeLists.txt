add_executable(
    runner 
    main.cpp 
    Mesh.cpp
    Elasticity.cpp
    DikeData.cpp
    InputData.cpp
    Schedule.cpp
    MassBalance.cpp
    Utils.cpp
    MagmaState.cpp
    Writers.cpp
)
target_link_libraries(
    runner
    Eigen3::Eigen
    nlohmann_json::nlohmann_json
    HighFive
)
target_include_directories(runner PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
install(TARGETS runner DESTINATION bin)