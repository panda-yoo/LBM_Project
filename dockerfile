# Use the latest version of Ubuntu as the base
FROM ubuntu:latest

# Avoid interactive prompts during installation
ENV DEBIAN_FRONTEND=noninteractive

# Update package lists and install dependencies (git is no longer needed for this)
RUN apt-get update && apt-get install -y \
    build-essential \
    cmake \
    python3

# Copy the local Palabos source code into the /opt/ directory in the image
COPY palabos /opt/palabos

# Set the working directory for when the container runs
WORKDIR /work