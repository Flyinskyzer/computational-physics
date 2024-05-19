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
    observe_sta result(n,T_spec,Tmax,themal,total,bin_start,bin_size,bin_gap,n_bins);
    result.run_MC_given_T();
    result.statistic();
}
    
    
    