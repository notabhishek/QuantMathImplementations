# Create a build directory (recommended)
mkdir cmake-build-release 
cd cmake-build-release

# Generate the build system
cmake .. -DCMAKE_CXX_COMPILER=clang++ -DCMAKE_BUILD_TYPE=Release

# Compile the project
make