add_executable(
    runner 
    main.cpp 
    Mesh.cpp
    Elasticity.cpp
    DikeData.cpp
    InputData.cpp
    Schedule.cpp
    MassBalance.cpp
)
target_link_libraries(
    runner
    Eigen3::Eigen
    nlohmann_json::nlohmann_json
)
install(TARGETS runner DESTINATION bin)