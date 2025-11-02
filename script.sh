
#!/bin/bash

# --- Configuration ---
# Use the first script argument as the program name, or default to 'my_example'.
EXECUTABLE_NAME=${1:-my_example}
BUILD_DIR="build"

echo "--> Cleaning up old build directory..."
rm -rf "$BUILD_DIR"

echo "--> Creating build directory and configuring with CMake..."
# Create the build directory, navigate into it, and run cmake.
# The '&&' operator ensures the script will stop if a command fails.


# mkdir "$BUILD_DIR" && cd "$BUILD_DIR" && cmake ..
mkdir "$BUILD_DIR" && cd "$BUILD_DIR" && cmake -DCMAKE_CXX_FLAGS="-fsanitize=address -g -O1" ..


# Check if cmake succeeded before trying to build.
if [ $? -ne 0 ]; then
    echo "❌ CMake configuration failed. Aborting."
    exit 1
fi

echo "--> Compiling the project with make..."
make -j4 

# Check if the executable was created successfully, then run it.
if [ -f "$EXECUTABLE_NAME" ]; then
    echo "✅ Build successful!"
    echo "--- Running '$EXECUTABLE_NAME' ---"
    # mpirun -np 1 ./"$EXECUTABLE_NAME"
    ./"$EXECUTABLE_NAME"

    echo "--- Program finished. ---"
else
    echo "❌ Build failed. Executable '$EXECUTABLE_NAME' not found."
    exit 1
fi
