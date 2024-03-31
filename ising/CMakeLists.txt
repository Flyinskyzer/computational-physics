# test framework repository
include(FetchContent)
FetchContent_Declare(
    catch2
  # GIT_REPOSITORY https://github.com/catchorg/Catch2.git
    SOURCE_DIR /Catch2
  GIT_TAG v3.1.0
)
FetchContent_MakeAvailable(Catch2)

# main (irrelevant now; temporary settings here)cd
cmake_minimum_required(VERSION 3.14)
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_COMPILER "D:/mingw64/bin/g++.exe")
project("Test of Catch2" CXX)
add_executable(a.exe C:\Users\dell\Desktop\计算物理\ising\computational-physics\ising\main.cpp)
target_compile_options(a.exe PRIVATE -Wall -Wextra)

# unit tests
enable_testing()
add_executable(UnitTests_catch2 C:\Users\dell\Desktop\Project_IsingModel\Project_IsingModel\ising\UnitTests.cpp)
target_compile_options(UnitTests_catch2 PRIVATE -Wall -Wextra)
target_link_libraries(UnitTests_catch2 Catch2::Catch2WithMain)
add_test(UnitTests UnitTests_catch2)
