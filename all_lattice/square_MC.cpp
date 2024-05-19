//#include"utils/isingsquare.h"
#include"MonteCarlosquare.h"
int main()
{
    int n,themal,total;
    double T_spec,Tmax;
    int bin_start,bin_size,bin_gap,n_bins;
    cin>>n>>T_spec>>Tmax>>themal>>total>>bin_start>>bin_size>>bin_gap>>n_bins;
    //MC result(n,T_spec,Tmax,themal,total);
    //result.run_MC_given_T();
    //result.run_MC();
    observe_sta result(n,T_spec,Tmax,themal,total,bin_start,bin_size,bin_gap,n_bins,1234567890);
    result.run_MC_given_T();
    result.statistic();
    result.calculate();
}
//边长 温度 最大温度 取样起点 取样终点 bin起点 bin大小 bin间隔 bin个数
//3 2 8 0 1000000 20000 800 200 80
    
    