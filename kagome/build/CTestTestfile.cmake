# CMake generated Testfile for 
# Source directory: /home/enrui/CP2024_code/computational-physics/kagome
# Build directory: /home/enrui/CP2024_code/computational-physics/kagome/build
# 
# This file includes the relevant testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
add_test(UnitTests "UnitTests_catch2")
set_tests_properties(UnitTests PROPERTIES  _BACKTRACE_TRIPLES "/home/enrui/CP2024_code/computational-physics/kagome/CMakeLists.txt;22;add_test;/home/enrui/CP2024_code/computational-physics/kagome/CMakeLists.txt;0;")
subdirs("_deps/catch2-build")