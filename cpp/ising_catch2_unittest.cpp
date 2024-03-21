#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <iostream>
#include "/home/enrui/CP2024_code/computational-physics/CPP/ising.hpp"

TEST_CASE("IsingSpin", "[single spin]") {
IsingSpin spin;

SECTION("spin state (initial)") {
    REQUIRE(spin._sz() == 1);
}
SECTION("set spin state as up (1)") {
    spin.set_up();
    REQUIRE(spin._sz() == 1);
}
    SECTION("set spin state as down (1)") {
    spin.set_dw();
REQUIRE(spin._sz() == -1);
}
    SECTION("set spin state as up (2)") {
    spin.set_sz(1);
REQUIRE(spin._sz() == 1);
}
SECTION("set spin state as down (2)") {
    spin.set_sz(-1);
    REQUIRE(spin._sz() == -1);
}
SECTION("spin flip once") {
    spin.flip();
    REQUIRE(spin._sz() == -1);
}
SECTION("spin flip twice") {
    spin.flip();
    spin.flip();
    REQUIRE(spin._sz() == 1);
} 
};

TEST_CASE("class_isingsystem_test","[set status and caculate energy and mag]")
{
    IsingSystem oneDspin(10);

    SECTION("const J") {
    REQUIRE (oneDspin._J() == -1.0);
    }
    SECTION("number of spins") {
    REQUIRE (oneDspin._n_spins_() == 10);
    }
    SECTION("max state") {
    REQUIRE (oneDspin.long_maxrep_state() == (std::pow(2,10)-1));
    }
    SECTION("initial") {
        for (int i=0;i<=9;++i){
    REQUIRE (oneDspin._sz(i) == 1);}
    }
    SECTION("set spin state as down (1)") {
    oneDspin.set_dw_spin(3);
    REQUIRE (oneDspin._sz(3)==-1);
    }
    SECTION("set spin state as down (2)") {
    oneDspin.set_dw_spin(5);
    REQUIRE (oneDspin._sz(5)==-1);
    }
    SECTION("set spin state as down (3)") {
    oneDspin.set_dw_spin(8);
    REQUIRE (oneDspin._sz(8)==-1);
    }
    SECTION("set spin state as up (1)") {
    oneDspin.set_up_spin(8);
    REQUIRE (oneDspin._sz(8) == 1);
    }
    SECTION("set spin state ") {
    oneDspin.set_spin(9,-1);
    REQUIRE (oneDspin._sz(9) == -1);
    }
    SECTION("spin flip once") {
    oneDspin.flip_spin(9);
    REQUIRE(oneDspin._sz(9) == -1);
    }
    SECTION("spin flip twice") {
    oneDspin.flip_spin(9);
    oneDspin.flip_spin(9);
    REQUIRE(oneDspin._sz(9) == 1);
    }
    SECTION("set state") {
    oneDspin.set_state_by_code(7);
    for (int i=0;i<=2;++i){
    REQUIRE (oneDspin._sz(i) == 1);
    };
    for (int i=3;i<=9;++i){
    REQUIRE (oneDspin._sz(i) == -1);
    };
    }
    SECTION("E") {
        oneDspin.set_state_by_code(7);
    REQUIRE(oneDspin.eval_energy_1D() == -6);
    }
    SECTION("mz") {
        oneDspin.set_state_by_code(7);
    REQUIRE(oneDspin.eval_mz() == -4);
    }
}