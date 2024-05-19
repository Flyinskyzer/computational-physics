# CMake generated Testfile for 
# Source directory: /home/enrui/CP2024_code/computational-physics/all_lattice
# Build directory: /home/enrui/CP2024_code/computational-physics/all_lattice/build
# 
# This file includes the relevant testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
add_test(isingsquaretest "isingsquaretest")
set_tests_properties(isingsquaretest PROPERTIES  _BACKTRACE_TRIPLES "/home/enrui/CP2024_code/computational-physics/all_lattice/CMakeLists.txt;25;add_test;/home/enrui/CP2024_code/computational-physics/all_lattice/CMakeLists.txt;0;")
add_test(kagometest "kagometest")
set_tests_properties(kagometest PROPERTIES  _BACKTRACE_TRIPLES "/home/enrui/CP2024_code/computational-physics/all_lattice/CMakeLists.txt;30;add_test;/home/enrui/CP2024_code/computational-physics/all_lattice/CMakeLists.txt;0;")
add_test(cubictest "cubictest")
set_tests_properties(cubictest PROPERTIES  _BACKTRACE_TRIPLES "/home/enrui/CP2024_code/computational-physics/all_lattice/CMakeLists.txt;35;add_test;/home/enrui/CP2024_code/computational-physics/all_lattice/CMakeLists.txt;0;")
add_test(honeytest "honeytest")
set_tests_properties(honeytest PROPERTIES  _BACKTRACE_TRIPLES "/home/enrui/CP2024_code/computational-physics/all_lattice/CMakeLists.txt;40;add_test;/home/enrui/CP2024_code/computational-physics/all_lattice/CMakeLists.txt;0;")
add_test(square_MC_test "square_MC_test")
set_tests_properties(square_MC_test PROPERTIES  _BACKTRACE_TRIPLES "/home/enrui/CP2024_code/computational-physics/all_lattice/CMakeLists.txt;45;add_test;/home/enrui/CP2024_code/computational-physics/all_lattice/CMakeLists.txt;0;")
subdirs("_deps/catch2-build")
