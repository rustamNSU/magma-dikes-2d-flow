add_library(utils INTERFACE)
target_include_directories(utils SYSTEM INTERFACE
    ${CMAKE_CURRENT_SOURCE_DIR}
)