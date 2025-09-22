# CMakeLists.txt
# Set the minimum required version of CMake
cmake_minimum_required(VERSION 3.10)

# Define the project name and specify the language
project(PalabosSim CXX)

# --- Project-wide Settings ---
# Set the C++ standard to C++17
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Add compiler options
add_compile_options(-Wno-narrowing)

# --- Palabos Library Definition ---
# Define the location of the Palabos source code
set(PALABOS_ROOT /opt/palabos)

# Find all Palabos source files.
# WARNING: Using file(GLOB...) can cause issues where new files aren't detected.
file(GLOB_RECURSE PALABOS_SRC
    "${PALABOS_ROOT}/src/*.cpp"
    "${PALABOS_ROOT}/externalLibraries/*.cpp"
)


# 1. Create a static library target for Palabos
add_library(palabos_lib STATIC ${PALABOS_SRC})

# Make sure our new library can find all necessary Palabos headers
target_include_directories(palabos_lib PUBLIC
    "${PALABOS_ROOT}/src"
    "${PALABOS_ROOT}/externalLibraries"
    "${PALABOS_ROOT}/include"
)

# --- Your Executable ---
# Create your executable from your personal source file
add_executable(PalabosSim palabos_test.cpp src/obstacle.cpp)

# Add your local 'include' directory to this target's include path
target_include_directories(PalabosSim PRIVATE include) # <-- THIS IS THE NEW LINE

# 2. Link your executable against the Palabos library and pthreads
target_link_libraries(PalabosSim PRIVATE palabos_lib pthread)


set_target_properties(PalabosSim PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}/data"
)




# cmake_minimum_required(VERSION 3.10)
# project(PalabosSim)

# add_compile_options(-Wno-narrowing)
# # Set the C++ standard
# set(CMAKE_CXX_STANDARD 17)
# set(CMAKE_CXX_STANDARD_REQUIRED ON)

# # Define the location of Palabos inside the container
# set(PALABOS_ROOT /opt/palabos)

# # Add include directories for both main source AND external libraries
# include_directories(${PALABOS_ROOT}/src)
# include_directories(${PALABOS_ROOT}/externalLibraries)
# include_directories(${PALABOS_ROOT}/include)

# target_include_directories("${PALABOS_ROOT}/include")
# # Find all Palabos source files (.cpp) from both directories
# file(GLOB_RECURSE PALABOS_SRC
#     "${PALABOS_ROOT}/src/*.cpp"
#     "${PALABOS_ROOT}/externalLibraries/*.cpp"
    
# )
# add_library(palabos_lib STATIC ${PALABOS_SRC})
# # Define our executable using your source file from the error message
# add_executable(PalabosSim src/obstacle.cpp )

# # Palabos uses threads, so we need to link against the pthreads library
# target_link_libraries(PalabosSim pthread)