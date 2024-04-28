#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include "kagome.h"

TEST_CASE("isingspinkagome","")
{
    Isingspinkagome kagomeunit;
    SECTION("setunit_by_0")
    {
        kagomeunit.setunit_by_code(0);
        REQUIRE(kagomeunit._unitstate()==0);
        REQUIRE(kagomeunit._unit_spinstate(0)==-1);
        REQUIRE(kagomeunit._unit_spinstate(1)==-1);
        REQUIRE(kagomeunit._unit_spinstate(2)==-1);
    }
    SECTION("setunit_by_1")
    {
        kagomeunit.setunit_by_code(1);
        REQUIRE(kagomeunit._unitstate()==1);
        REQUIRE(kagomeunit._unit_spinstate(0)==1);
        REQUIRE(kagomeunit._unit_spinstate(1)==-1);
        REQUIRE(kagomeunit._unit_spinstate(2)==-1);
    }
    SECTION("setunit_by_2")
    {
        kagomeunit.setunit_by_code(2);
        REQUIRE(kagomeunit._unitstate()==2);
        REQUIRE(kagomeunit._unit_spinstate(0)==-1);
        REQUIRE(kagomeunit._unit_spinstate(1)==1);
        REQUIRE(kagomeunit._unit_spinstate(2)==-1);
    }
    SECTION("setunit_by_3")
    {
        kagomeunit.setunit_by_code(3);
        REQUIRE(kagomeunit._unitstate()==3);
        REQUIRE(kagomeunit._unit_spinstate(0)==1);
        REQUIRE(kagomeunit._unit_spinstate(1)==1);
        REQUIRE(kagomeunit._unit_spinstate(2)==-1);
    }
}
TEST_CASE("Kagome system","")
{
    Isingsystem system({2,2});
    SECTION("set unit int units")
    {
        system.set_unitstate_in_total(3625);
        REQUIRE(system._state(0)==1);
        REQUIRE(system._state(1)==5);
        REQUIRE(system._state(2)==0);
        REQUIRE(system._state(3)==7);
        REQUIRE(system._spin_state(0,0)==1);
        REQUIRE(system._spin_state(0,1)==-1);
        REQUIRE(system._spin_state(0,2)==-1);
        REQUIRE(system._spin_state(1,0)==1);
        REQUIRE(system._spin_state(1,1)==-1);
        REQUIRE(system._spin_state(1,2)==1);
        REQUIRE(system._spin_state(2,0)==-1);
        REQUIRE(system._spin_state(2,1)==-1);
        REQUIRE(system._spin_state(2,2)==-1);
        REQUIRE(system._spin_state(3,0)==1);
        REQUIRE(system._spin_state(3,1)==1);
        REQUIRE(system._spin_state(3,2)==1);
    }
    SECTION("coordinate<->site")
    {
        REQUIRE(system.site_index({1,1})==3);
        vector<int> const a={0,1};
        REQUIRE(system.lattice_coordinate(2)==a);
    }
    
}
TEST_CASE("shift","")
{
    Isingsystem model2({6,6});
    SECTION("shift")
    {
        const vector<int> l_c2= {3,3};
        const vector<int> l_c2_n_x= {2,3};
        const vector<int> l_c2_p_x= {4,3};
        const vector<int> l_c2_n_y= {3,2};
        const vector<int> l_c2_p_y= {3,4};
        const vector<int> l_c2_n_x_p_y={2,4};
        const vector<int> l_c2_p_x_n_y={4,2};
        REQUIRE(model2.shift_neg_x(l_c2)==l_c2_n_x);
        REQUIRE(model2.shift_pos_x(l_c2)==l_c2_p_x);
        REQUIRE(model2.shift_neg_y(l_c2)==l_c2_n_y);
        REQUIRE(model2.shift_pos_y(l_c2)==l_c2_p_y);
        REQUIRE(model2.shift_neg_x_pos_y(l_c2)==l_c2_n_x_p_y);
        REQUIRE(model2.shift_pos_x_neg_y(l_c2)==l_c2_p_x_n_y);
    };
}

TEST_CASE("network")
{
    Isingsystem model({2,2});
    REQUIRE(model.NN(0,0)==1);
    REQUIRE(model.NN(0,1)==2);
    REQUIRE(model.NN(0,2)==1);
    REQUIRE(model.NN(0,3)==2);
    REQUIRE(model.NN(0,4)==3);
    REQUIRE(model.NN(0,5)==3);
}

TEST_CASE("calculate","")
{
    Isingsystem system({2,2});
    system.set_unitstate_in_total(3625);
    REQUIRE(system.momentum_eval()==0);
    REQUIRE(system.energy_eval()==0);
}