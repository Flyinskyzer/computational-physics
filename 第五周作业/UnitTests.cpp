#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include "ising.hpp"
#include<vector>

// 对 Isingspin 类进行单元测试
TEST_CASE("Isingspin class tests","[single spin]") {
    Isingspin spin;

    SECTION("Initial size should be -1") {
        REQUIRE(spin._sz() == -1);
    }

    SECTION("Setting up spin should change the spin state to 1") {
        spin.set_up();
        REQUIRE(spin._sz() == 1);
    }

    SECTION("Setting down spin should change the spin to -1"){
        spin.set_dw();
        REQUIRE(spin._sz() == -1);
    }

    SECTION("Flipping spin should change its state") {
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

    SECTION("set the spin state as up"){
        spin.set_sz(1);
        REQUIRE(spin._sz()==1);
    }
    SECTION("set the spin state as down"){
        spin.set_sz(-1);
        REQUIRE(spin._sz()==-1);
    }
}


TEST_CASE("Isingsystem class tests") {
    Isingsystem system(4);  // 假设有4个自旋

    SECTION("Check the number of spins") {
        REQUIRE(system._n_spins() == 4);
    }
    SECTION("Setting up spins should initialize all spins") {
        system.set_up_spin(2);
        REQUIRE(system._spin()[2]._sz()==1);
    }   


    SECTION("return the J"){
        REQUIRE(system._J()==-1.0);
    }
    SECTION("set the single spin"){
        system.set_spin(2,1);
        REQUIRE(system._spin()[2]._sz()==1);
    }
    SECTION("set the spin state by a number")
    {
        system.set_state_by_code(3);
       REQUIRE(system._spin()[0]._sz()==1);
       REQUIRE(system._spin()[1]._sz()==1);
       REQUIRE(system._spin()[2]._sz()==-1);
       REQUIRE(system._spin()[3]._sz()==-1);
    }
    SECTION("calculate the energy"){
      system.set_state_by_code(3);
      REQUIRE(system.eval_energy_1D()==0);
    }
    SECTION("calculat the agular momentum"){
      system.set_state_by_code(3);
      REQUIRE(system.eval_mz_1D()==0);  
    }
}