#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include "isingsquare.h"

TEST_CASE("isingsquare","coordinate->site_idx")
{
    const vector<int>sys_size={6,6};
    Isingsystem_square model(sys_size);

    SECTION("coordinate->site_idx")
    {
        const vector<int> l_c= {3,4};
        REQUIRE(model.site_index(l_c)==27);
        REQUIRE(model.latice_coordinate(27)== l_c);
    };
    SECTION("NN")
    {
        const int i = 21;
        REQUIRE(model.NN(i,0)==22);
        REQUIRE(model.NN(i,1)==27);
        REQUIRE(model.NN(i,2)==20);
        REQUIRE(model.NN(i,3)==15);
    };
}
TEST_CASE("shift","")
{
     const vector <int> sys_size2={6,6};
     Isingsystem_square model2(sys_size2);
      SECTION("coordinate->site_idx") 
    {
        const vector<int> l_c2= {3,3};
        const vector<int> l_c21= {0,3};
        const vector<int> l_c2_n_x= {2,3};
        const vector<int> l_c2_p_x= {4,3};
        const vector<int> l_c2_n_y= {3,2};
        const vector<int> l_c2_p_y= {3,4};
        const vector<int> l_c21_n_x= {5,3};
        REQUIRE(model2.shift_neg_x(l_c2)==l_c2_n_x);
        REQUIRE(model2.shift_pos_x(l_c2)==l_c2_p_x);
        REQUIRE(model2.shift_neg_y(l_c2)==l_c2_n_y);
        REQUIRE(model2.shift_pos_y(l_c2)==l_c2_p_y);
        REQUIRE(model2.shift_neg_x(l_c21)==l_c21_n_x);
    };
}
TEST_CASE("energy and momentum","")
{
    const vector <int> sys_size={6,6};
    Isingsystem_square model(sys_size);
    model.set_state_by_code(pow(2,36)-1);
    SECTION("E")
    {
        REQUIRE(model.energy_eval()==-72);
    };
    SECTION("M")
    {
        REQUIRE(model.momentum_eval()==360);
    };
}