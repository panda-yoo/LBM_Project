cmake_minimum_required(VERSION 3.10)
project(MyPalabosExample LANGUAGES CXX)

# === Compiler setup ===
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_BUILD_TYPE Release)

# === Find MPI ===
find_package(MPI REQUIRED)
include_directories(${MPI_INCLUDE_PATH})

# === Palabos include + library paths ===
include_directories(
    /usr/local/include/palabos/src
    /usr/local/include/palabos/externalLibraries
    ${CMAKE_CURRENT_SOURCE_DIR}/include
)

link_directories(/usr/local/lib)

# === Collect sources ===
file(GLOB_RECURSE SOURCES "src/*.cpp")

# === Create executable ===
add_executable(example ${SOURCES} main.cpp)

# === Ensure we use MPI compile + link options ===
target_compile_options(example PRIVATE ${MPI_CXX_COMPILE_FLAGS})
target_link_libraries(example PRIVATE palabos ${MPI_CXX_LIBRARIES})
