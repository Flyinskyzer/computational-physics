#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <iostream>
#include <vector>
#include "honey.h"

TEST_CASE("class isingspin_honey_unit", "") 
{
    isingspin_honey_unit model;

    SECTION("set by num") 
    {
        model.unit_set_sz_by_num(0);
        vector <int> a1={-1,-1},a2={1,-1},a3={-1,1},a4={1,1};
        REQUIRE(model.unit_sz() == a1);
        REQUIRE(model.unit_mz() == -2);
        REQUIRE(model.unit_e() == -1);
        model.unit_set_sz_by_num(1);
        REQUIRE(model.unit_sz() == a2);
        REQUIRE(model.unit_mz() == 0);
        REQUIRE(model.unit_e() == 1);
        model.unit_set_sz_by_num(2);
        REQUIRE(model.unit_sz() == a3);
        REQUIRE(model.unit_mz() == 0);
        REQUIRE(model.unit_e() == 1);
        model.unit_set_sz_by_num(3);
        REQUIRE(model.unit_sz() == a4);
        REQUIRE(model.unit_mz() == 2);
        REQUIRE(model.unit_e() == -1);
    };
};

TEST_CASE("Ising_honey", "coordinate->site_idx") 
{
    const vector <int> sys_size={6,6};
    ising_system_honey model(sys_size);

    SECTION("coordinate->site_idx") 
    {
        const vector<int> l_c= {3,4};
        const vector<int> l_cx= {3,3};
        REQUIRE(model.site_index(l_c)==27);
        REQUIRE(model.lattice_coordinate(27)== l_c);
    };
    SECTION("NN")
    {
        constexpr int i = 21;
        REQUIRE(model.NN(i,0)==22);
        REQUIRE(model.NN(i,1)==27);
        REQUIRE(model.NN(i,2)==20);
        REQUIRE(model.NN(i,3)==15);
    }

};


TEST_CASE("shift", "") 
{
    const vector <int> sys_size2={6,6};
    ising_system_honey model2(sys_size2);

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
};

TEST_CASE("SET STATE AND CACULATE sz,mz,h,c", "") 
{
    SECTION("SIZE(3,3),NUM2134")
    {vector<int> sys_size = {3, 3};
    ising_system_honey model( sys_size);
    double N = 2*model._n_unitspins_();
    double H0 = model.eval_energy();
    long long max_state = model._maxrep_state();
    model.set_state_by_num(2134);
    vector <int> a=model._sz(5);
    REQUIRE(model.eval_energy()==-9);
    REQUIRE(model.eval_Mz());
    REQUIRE(N==18);
    REQUIRE(a[0]==-1);
    REQUIRE(a[1]==1);}
}