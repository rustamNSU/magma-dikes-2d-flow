add_library(utils ViscosityModels.cpp)
target_include_directories(utils PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
install(TARGETS utils DESTINATION lib)