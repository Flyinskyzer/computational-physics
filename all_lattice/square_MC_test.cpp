#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include"MonteCarlosquare.h"
#include <iostream>
#include <random>

int generateRandomNumber(int min, int max) {
    std::mt19937 generator(std::random_device{}());
    std::uniform_int_distribution<int> distribution(min, max);
    return distribution(generator);
}

TEST_CASE("square")
{
    int n,themal,total;
    double T_spec,Tmax;
    int bin_start,bin_size,bin_gap,n_bins;
    n=4;T_spec=2;Tmax=8;themal=0;total=1000000;bin_start=20000;bin_size=800;bin_gap=200;n_bins=80;
    double h=-1.75538;
    double m2=0.868922;
    int h_count=0;
    int m2_count=0;
    int num=100;

    for(int i=0;i<num;i++){
        observe_sta result(n,T_spec,Tmax,themal,total,bin_start,bin_size,bin_gap,n_bins,1234567890+i/*generateRandomNumber(1,10000)*/);
        result.run_MC_given_T();
        result.statistic();
        result.calculate();
        double h_ave=result._h_average();
        double h_sig=result._h_sigama();
        double m2_ave=result._m2_average();
        double m2_sig=result._m2_sigama();
        if(h_ave-1.96*h_sig/8.9<=h&&h<=h_ave+1.96*h_sig/8.9){ h_count++;};
        if(m2_ave-1.96*m2_sig/8.9<=m2&&m2<=m2_ave+1.96*m2_sig/8.9){ m2_count++;};
    }

    SECTION("test of H"){
        REQUIRE(h_count<=99);
        REQUIRE(h_count>=91);
    }
    SECTION("test of M2"){
        REQUIRE(m2_count<=99);
        REQUIRE(m2_count>=91);
    }
}