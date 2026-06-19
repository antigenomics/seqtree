cmake_minimum_required(VERSION 3.20)
project(seqtree LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()
if(NOT MSVC)
  set(CMAKE_CXX_FLAGS_RELEASE "-O3")
endif()  # MSVC Release already uses /O2; -O3 is not a valid MSVC flag

option(SEQTREE_TESTS "Build C++ tests" OFF)
option(SEQTREE_BENCH "Build C++ benchmarks" OFF)
option(SEQTREE_PYTHON "Build the pybind11 module" OFF)

find_package(Threads REQUIRED)

add_library(seqtree_core STATIC
  src/codec.cpp
  src/substitution_matrix.cpp
  src/positional_matrix.cpp
  src/trie.cpp
  src/index.cpp
  src/engine_seqtm.cpp
  src/engine_seqtrie.cpp
  src/searcher.cpp
)
target_include_directories(seqtree_core PUBLIC include PRIVATE src)
target_link_libraries(seqtree_core PUBLIC Threads::Threads)

if(SEQTREE_PYTHON)
  find_package(pybind11 CONFIG REQUIRED)
  pybind11_add_module(_core src/_bindings.cpp)
  target_link_libraries(_core PRIVATE seqtree_core)
  install(TARGETS _core DESTINATION seqtree)
endif()

if(SEQTREE_TESTS)
  enable_testing()
  add_executable(seqtree_tests
    tests/cpp/test_codec.cpp
    tests/cpp/test_matrix.cpp
    tests/cpp/test_trie.cpp
    tests/cpp/test_engines.cpp
    tests/cpp/test_edge.cpp
    tests/cpp/test_serialize.cpp
    tests/cpp/test_positional.cpp
  )
  target_include_directories(seqtree_tests PRIVATE tests/cpp src)
  target_link_libraries(seqtree_tests PRIVATE seqtree_core)
  add_test(NAME seqtree_tests COMMAND seqtree_tests)
endif()

if(SEQTREE_BENCH)
  add_executable(seqtree_bench bench/bench_seqtree.cpp)
  target_link_libraries(seqtree_bench PRIVATE seqtree_core)
endif()
