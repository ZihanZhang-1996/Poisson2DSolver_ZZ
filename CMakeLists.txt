# CMakeLists.txt (root)

cmake_minimum_required(VERSION 3.16)
project(pde2d LANGUAGES CXX)

# ---- Modern CMake: Boost policy ----
cmake_policy(SET CMP0167 NEW)   # Suppress FindBoost deprecation warning

# ---- Options ----
option(PDE2D_BUILD_SHARED "Build pde2d as a shared library" ON)

# ---- Dependencies ----
# Eigen (header-only)
find_package(Eigen3 REQUIRED)

# Modern Boost (find via config, imported targets)
find_package(Boost 1.70 REQUIRED COMPONENTS math_c99)

# ---- Library: pde2d ----
set(PDE2D_SOURCES
  src/pde2d/Mesh2D.cpp
  src/pde2d/problem_syb.cpp
  src/pde2d/fem_2D.cpp
  src/pde2d/solver_newton.cpp
  src/pde2d/io.cpp
  # src/pde2d/quadrature_basic.cpp
  src/pde2d/quadrature_boost.cpp
)

if(PDE2D_BUILD_SHARED)
  add_library(pde2d SHARED ${PDE2D_SOURCES})
else()
  add_library(pde2d STATIC ${PDE2D_SOURCES})
endif()

add_library(pde2d::pde2d ALIAS pde2d)

target_compile_features(pde2d PUBLIC cxx_std_17)


target_include_directories(pde2d
  PUBLIC
    ${PROJECT_SOURCE_DIR}/include
    ${PROJECT_SOURCE_DIR}/third_party    
)

target_link_libraries(pde2d
  PUBLIC
    Eigen3::Eigen
    Boost::math_c99
)

# (Optional) nice warnings
if(MSVC)
  target_compile_options(pde2d PRIVATE /W4 /permissive-)
else()
  target_compile_options(pde2d PRIVATE -Wall -Wextra -Wpedantic)
endif()

# ---- Executable: client ----
add_executable(pde2d_client
  src/client/main.cpp
)

target_link_libraries(pde2d_client
  PRIVATE
    pde2d::pde2d
)

target_compile_features(pde2d_client PRIVATE cxx_std_17)

# Create "input" folder in the build directory
file(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/input)
file(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/output)

# Copy example file into it
configure_file(
  ${PROJECT_SOURCE_DIR}/data/example_input.json
  ${CMAKE_BINARY_DIR}/input/example_input.json
  COPYONLY
)


# ---- Tests ----
enable_testing()

add_executable(test_basic
  test/test_basic.cpp   
)


target_link_libraries(test_basic
  PRIVATE
    pde2d::pde2d
)

target_compile_features(test_basic PRIVATE cxx_std_17)

add_test(NAME run_test_basic COMMAND test_basic)