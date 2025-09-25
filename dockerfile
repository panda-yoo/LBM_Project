# Use a specific version of Ubuntu for reproducibility
FROM ubuntu:22.04

# Avoid interactive prompts during installation
ENV DEBIAN_FRONTEND=noninteractive

# Set up the Python virtual environment path
ENV VENV_PATH=/opt/venv

# Add the venv's bin directory to the system's PATH.
# This makes commands like 'python' and 'pip' use the venv's executables by default.
ENV PATH="$VENV_PATH/bin:$PATH"

# Update package lists and install dependencies in a single layer
RUN apt-get update && apt-get install -y \
    build-essential \
    cmake \
    python3 \
    python3-pip \
    python3-venv \
    && rm -rf /var/lib/apt/lists/*

# Create the virtual environment
RUN python3 -m venv $VENV_PATH

# Copy ONLY the requirements file first to leverage Docker's layer cache
COPY requirements.txt .

# Install Python packages into the virtual environment
# We don't need to 'activate' it because we added it to the PATH
RUN pip install --no-cache-dir -r requirements.txt

# Set the working directory for the project
WORKDIR /LBM_Project

# Copy all the project files
COPY palabos /opt/palabos
COPY src ./src
COPY include ./include
COPY analysis.py .
COPY CMakeLists.txt .
COPY palabos_test.cpp .