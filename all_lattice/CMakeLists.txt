include(FetchContent)
FetchContent_Declare(catch2
  GIT_REPOSITORY https://github.com/catchorg/Catch2.git
  GIT_TAG v3.1.0
)
FetchContent_MakeAvailable(Catch2)
enable_testing()

cmake_minimum_required(VERSION 3.14)
set(CMAKE_CXX_STANDARD 14)
project(CP2024 CXX)
add_executable(kagome kagome.cpp)
target_compile_options(kagome PRIVATE -Wall -Wextra -O3 -DNDEBUG)
add_executable(square square.cpp)
target_compile_options(square PRIVATE -Wall -Wextra -O3 -DNDEBUG)
add_executable(cubic cubic.cpp)
target_compile_options(cubic PRIVATE -Wall -Wextra -O3 -DNDEBUG)
add_executable(honey honey.cpp)
target_compile_options(honey PRIVATE -Wall -Wextra -O3 -DNDEBUG)


add_executable(isingsquaretest utils/isingsquaretest.cpp)
target_compile_options(isingsquaretest PRIVATE -Wall -Wextra)
target_link_libraries(isingsquaretest Catch2::Catch2WithMain)
add_test(isingsquaretest isingsquaretest)

add_executable(kagometest utils/kagometest.cpp)
target_compile_options(kagometest PRIVATE -Wall -Wextra)
target_link_libraries(kagometest Catch2::Catch2WithMain)
add_test(kagometest kagometest)

add_executable(cubictest utils/cubictest.cpp)
target_compile_options(cubictest PRIVATE -Wall -Wextra)
target_link_libraries(cubictest Catch2::Catch2WithMain)
add_test(cubictest cubictest)

add_executable(honeytest utils/honeytest.cpp)
target_compile_options(honeytest PRIVATE -Wall -Wextra)
target_link_libraries(honeytest Catch2::Catch2WithMain)
add_test(honeytest honeytest)
