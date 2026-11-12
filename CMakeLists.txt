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

# ΨLab for the terminal (Linux, macOS): plain-text input file -> plain-text output file.
add_library(psilab_text STATIC backend/cli/psilab_text.cpp)
target_include_directories(psilab_text PUBLIC backend/cli)
target_link_libraries(psilab_text PUBLIC macqc)

add_executable(psilab backend/cli/psilab.cpp)
target_link_libraries(psilab PRIVATE psilab_text)

enable_testing()
add_executable(macqc_tests backend/tests/test_macqc.cpp)
target_link_libraries(macqc_tests PRIVATE macqc)
target_include_directories(macqc_tests PRIVATE third_party)
target_compile_definitions(macqc_tests PRIVATE
  MACQC_TEST_DATA_DIR="${CMAKE_CURRENT_SOURCE_DIR}/backend/tests/data")
add_test(NAME macqc_tests COMMAND macqc_tests)
add_test(NAME cli_smoke COMMAND macqc_cli --print 1)
add_test(NAME cli_json COMMAND macqc_cli --json --n 2)
add_test(NAME cli_invalid COMMAND macqc_cli --R -1)
set_tests_properties(cli_invalid PROPERTIES WILL_FAIL TRUE)

add_executable(psilab_tests backend/tests/test_psilab.cpp)
target_link_libraries(psilab_tests PRIVATE psilab_text)
add_test(NAME psilab_tests COMMAND psilab_tests)
add_test(NAME psilab_example
  COMMAND psilab ${CMAKE_CURRENT_SOURCE_DIR}/examples/heh.inp ${CMAKE_CURRENT_BINARY_DIR}/heh.out)
add_test(NAME psilab_missing_input COMMAND psilab does-not-exist.inp)
set_tests_properties(psilab_missing_input PROPERTIES WILL_FAIL TRUE)
