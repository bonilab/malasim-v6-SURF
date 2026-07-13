#!/usr/bin/env bash

set -euo pipefail

# The script must be stored in the project root directory.
PROJECT_ROOT="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
BUILD_DIR="$PROJECT_ROOT/build"
BIN_DIR="$BUILD_DIR/bin"

cd "$PROJECT_ROOT"

# Safety check to avoid cleaning the wrong directory.
if [[ ! -f "$PROJECT_ROOT/Makefile" ]]; then
    echo "Error: Makefile not found in:"
    echo "  $PROJECT_ROOT"
    echo "Place this script in the project root directory."
    exit 1
fi

echo "Project root: $PROJECT_ROOT"

# Clone vcpkg only if it does not already exist.
if [[ ! -d "$PROJECT_ROOT/vcpkg/.git" ]]; then
    echo "Cloning vcpkg..."
    git clone https://github.com/microsoft/vcpkg.git "$PROJECT_ROOT/vcpkg"
else
    echo "Using existing vcpkg directory."
fi

export VCPKG_ROOT="$PROJECT_ROOT/vcpkg"

make setup-vcpkg

echo "Cleaning previous build files..."

if [[ -d "$BUILD_DIR" ]]; then
    # Protect against build/bin unexpectedly being a file or symbolic link.
    if [[ -e "$BIN_DIR" && ! -d "$BIN_DIR" ]]; then
        echo "Error: $BIN_DIR exists but is not a directory."
        exit 1
    fi

    # Remove everything directly inside build/ except build/bin/.
    find "$BUILD_DIR" \
        -mindepth 1 \
        -maxdepth 1 \
        ! -name "bin" \
        -exec rm -rf -- {} +

    # Preserve all other files and folders in build/bin/, but remove
    # executables that will be rebuilt.
    rm -f -- \
        "$BIN_DIR/DxGGenerator" \
        "$BIN_DIR/MalaSim" \
        "$BIN_DIR/malasim_test"
else
    mkdir -p "$BIN_DIR"
fi

# Ensure build/bin exists even on the first build.
mkdir -p "$BIN_DIR"

echo "Preserved contents of build/bin/:"
find "$BIN_DIR" -mindepth 1 -maxdepth 1 -printf '  %f\n' 2>/dev/null || true

echo "Generating Release build..."

make generate \
    BUILD_TYPE=Release \
    BUILD_TESTS=ON \
    ENABLE_COVERAGE=OFF

echo "Building project..."

make build BUILD_TYPE=Release

echo "Running tests..."

make test

echo "Build completed successfully."