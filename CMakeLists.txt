cmake_minimum_required(VERSION 3.16)
project(MacQC LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE Release)
endif()

# Keep floating-point arithmetic identical to the reference (no fused multiply-add).
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
  add_compile_options(-ffp-contract=off -Wall -Wextra)
endif()

add_library(macqc STATIC
  backend/src/linalg.cpp
  backend/src/integrals.cpp
  backend/src/report.cpp
  backend/src/scf.cpp
  backend/src/fci.cpp
  backend/src/calculation.cpp
  backend/src/api.cpp)
target_include_directories(macqc PUBLIC backend/include)
target_include_directories(macqc PRIVATE third_party)

add_executable(macqc_cli backend/cli/main.cpp)
target_link_libraries(macqc_cli PRIVATE macqc)
set_target_properties(macqc_cli PROPERTIES OUTPUT_NAME macqc)
