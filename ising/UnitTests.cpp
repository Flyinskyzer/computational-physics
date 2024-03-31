#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include "Ising.h"

// 对 Isingspin 类进行单元测试
TEST_CASE("Isingspin class tests") {
    Isingspin spin;

    SECTION("Initial size should be -1") {
        REQUIRE(spin._sz() == -1);
    }

    SECTION("Setting up spin should change size to 1") {
        spin.set_up();
        REQUIRE(spin._sz() == 1);
    }

    SECTION("Flipping spin should change its size") {
        spin.set_up();
        int old_sz = spin._sz();
        spin.flip();
        REQUIRE(spin._sz() == -old_sz);
    }

    SECTION("Flipping spin multiple times should revert to original size") {
        spin.set_up();
        int old_sz = spin._sz();
        spin.flip();
        spin.flip();
        REQUIRE(spin._sz() == old_sz);
    }

    // 其他 Isingspin 类函数的单元测试
    // ...
}


TEST_CASE("Isingsystem class tests") {
    Isingsystem system(4);  // 假设有4个自旋

    SECTION("Check the number of spins") {
        REQUIRE(system._n_spins() == 4);
    }

    SECTION("Setting up spins should initialize all spins") {
        system.set_up_spins();
        bool all_up = true;
        for (int i = 0; i < system._n_spins(); ++i) {
            if (system._sz(i) != 1) {
                all_up = false;
                break;
            }
        }
        REQUIRE(all_up == true);
    }

    SECTION("Flipping a spin should change its size") {
        system.set_up_spins();
        int old_sz = system._sz(0);
        system.flip_spin(0);
        REQUIRE(system._sz(0) == -old_sz);
    }

    SECTION("Flipping a spin multiple times should revert to original size") {
        system.set_up_spins();
        int old_sz = system._sz(0);
        system.flip_spin(0);
        system.flip_spin(0);
        REQUIRE(system._sz(0) == old_sz);
    }

    // 其他 Isingsystem 类函数的单元测试
    // ...
}
