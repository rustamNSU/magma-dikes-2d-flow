add_executable(
    runner 
    main.cpp 
    Mesh.cpp
    Elasticity.cpp
    DikeData.cpp
)
target_link_libraries(
    runner
    Eigen3::Eigen
)
install(TARGETS runner DESTINATION bin)