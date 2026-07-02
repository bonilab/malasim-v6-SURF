# Collect all .cpp and .h files in the current directory
file(GLOB DXG_GENERATOR_SOURCES
        "${CMAKE_CURRENT_SOURCE_DIR}/*.cpp"
        "${CMAKE_CURRENT_SOURCE_DIR}/*.h"
)

# Define the DxGGenerator executable
add_executable(DxGGenerator
        ${DXG_GENERATOR_SOURCES}
)

find_package(fmt CONFIG REQUIRED)
find_package(GSL REQUIRED)
find_package(yaml-cpp CONFIG REQUIRED)
find_package(spdlog REQUIRED)
find_package(date CONFIG REQUIRED)
find_package(CLI11 CONFIG REQUIRED)
find_package(unofficial-sqlite3 CONFIG REQUIRED)

# Ensure MalaSimCore is built first
add_dependencies(DxGGenerator MalaSimCore)

# Link MalaSimCore to DxGGenerator
target_link_libraries(DxGGenerator  PUBLIC MalaSimCore
)

# Include the directory containing MalaSimCore headers
target_include_directories(DxGGenerator PRIVATE ${CMAKE_SOURCE_DIR}/src)
