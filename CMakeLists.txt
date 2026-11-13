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

find_package(Threads REQUIRED)

add_library(macqc STATIC
  backend/src/linalg.cpp
  backend/src/integrals.cpp
  backend/src/report.cpp
  backend/src/scf.cpp
  backend/src/fci.cpp
  backend/src/molecule.cpp
  backend/src/basis.cpp
  backend/src/gto_integrals.cpp
  backend/src/rhf.cpp
  backend/src/molecular.cpp
  backend/src/gemm.cpp
  backend/src/mo_integrals.cpp
  backend/src/ci.cpp
  backend/src/cc_general.cpp
  backend/src/cc_closed_shell.cpp
  backend/src/cc_estimate.cpp
  backend/src/job.cpp
  backend/src/calculation.cpp
  backend/src/api.cpp)
target_include_directories(macqc PUBLIC backend/include)
target_include_directories(macqc PRIVATE third_party)
target_link_libraries(macqc PUBLIC Threads::Threads)
# Basis-set files (Gaussian94, from basissetexchange.org): the source tree is the default
# location for builds that are not installed; installed programs use <prefix>/share/psilab/basis.
target_compile_definitions(macqc PRIVATE MACQC_BASIS_DIR="${CMAKE_CURRENT_SOURCE_DIR}/basis")

add_executable(macqc_cli backend/cli/main.cpp)
target_link_libraries(macqc_cli PRIVATE macqc)
set_target_properties(macqc_cli PROPERTIES OUTPUT_NAME macqc)

# ΨLab for the terminal (Linux, macOS): plain-text input file -> plain-text output file.
add_library(psilab_text STATIC backend/cli/psilab_text.cpp)
target_include_directories(psilab_text PUBLIC backend/cli)
target_link_libraries(psilab_text PUBLIC macqc)

add_executable(psilab backend/cli/psilab.cpp)
target_link_libraries(psilab PRIVATE psilab_text)
# "cmake --install build" puts psilab into <prefix>/bin (default /usr/local/bin), so it runs as "psilab".
install(TARGETS psilab RUNTIME DESTINATION bin)
install(DIRECTORY basis/ DESTINATION share/psilab/basis FILES_MATCHING PATTERN "*.gbs" PATTERN "*.md")

enable_testing()
add_executable(macqc_tests backend/tests/test_macqc.cpp)
target_link_libraries(macqc_tests PRIVATE macqc)
target_include_directories(macqc_tests PRIVATE third_party)
target_compile_definitions(macqc_tests PRIVATE
  MACQC_TEST_DATA_DIR="${CMAKE_CURRENT_SOURCE_DIR}/backend/tests/data")
add_test(NAME macqc_tests COMMAND macqc_tests)
# ΨLab 3.0: basis sets, integrals, RHF (against PySCF and published HF limits), API.
add_executable(macqc_rhf_tests backend/tests/test_rhf.cpp)
target_link_libraries(macqc_rhf_tests PRIVATE macqc)
target_include_directories(macqc_rhf_tests PRIVATE third_party)
target_compile_definitions(macqc_rhf_tests PRIVATE
  MACQC_TEST_DATA_DIR="${CMAKE_CURRENT_SOURCE_DIR}/backend/tests/data")
add_test(NAME macqc_rhf_tests COMMAND macqc_rhf_tests)
set_tests_properties(macqc_rhf_tests PROPERTIES
  ENVIRONMENT "PSILAB_BASIS_DIR=${CMAKE_CURRENT_SOURCE_DIR}/basis" TIMEOUT 3600)
add_executable(macqc_fci_tests backend/tests/test_fci.cpp)
target_link_libraries(macqc_fci_tests PRIVATE macqc)
target_include_directories(macqc_fci_tests PRIVATE third_party)
target_compile_definitions(macqc_fci_tests PRIVATE
  MACQC_TEST_DATA_DIR="${CMAKE_CURRENT_SOURCE_DIR}/backend/tests/data")
add_test(NAME macqc_fci_tests COMMAND macqc_fci_tests)
set_tests_properties(macqc_fci_tests PROPERTIES
  ENVIRONMENT "PSILAB_BASIS_DIR=${CMAKE_CURRENT_SOURCE_DIR}/basis" TIMEOUT 7200)
# ΨLab 3.1: coupled cluster (brute-force references, PySCF, Full CI).
add_executable(macqc_cc_tests backend/tests/test_cc.cpp)
target_link_libraries(macqc_cc_tests PRIVATE macqc)
target_include_directories(macqc_cc_tests PRIVATE third_party)
target_compile_definitions(macqc_cc_tests PRIVATE
  MACQC_TEST_DATA_DIR="${CMAKE_CURRENT_SOURCE_DIR}/backend/tests/data")
add_test(NAME macqc_cc_tests COMMAND macqc_cc_tests)
set_tests_properties(macqc_cc_tests PROPERTIES
  ENVIRONMENT "PSILAB_BASIS_DIR=${CMAKE_CURRENT_SOURCE_DIR}/basis" TIMEOUT 3600)
add_test(NAME cli_smoke COMMAND macqc_cli --print 1)
add_test(NAME cli_json COMMAND macqc_cli --json --n 2)
add_test(NAME cli_invalid COMMAND macqc_cli --R -1)
set_tests_properties(cli_invalid PROPERTIES WILL_FAIL TRUE)

add_executable(psilab_tests backend/tests/test_psilab.cpp)
target_link_libraries(psilab_tests PRIVATE psilab_text)
add_test(NAME psilab_tests COMMAND psilab_tests)
add_test(NAME psilab_example
  COMMAND psilab ${CMAKE_CURRENT_SOURCE_DIR}/examples/heh.inp ${CMAKE_CURRENT_BINARY_DIR}/heh.out)
add_test(NAME psilab_example_h2o
  COMMAND psilab ${CMAKE_CURRENT_SOURCE_DIR}/examples/h2o.inp ${CMAKE_CURRENT_BINARY_DIR}/h2o.out)
add_test(NAME psilab_missing_input COMMAND psilab does-not-exist.inp)
set_tests_properties(psilab_missing_input PROPERTIES WILL_FAIL TRUE)
