#ifndef MONTECARLO
#define MONTECARLO
#include<iostream>
#include<vector>
#include<random>
#include<string>
#include<cmath>
#include<sstream>
#include"utils/isingsquare.h"
#include<fstream>

using namespace std;
class RandomNumberGenerator{
private:
     const int seed_init;
     mt19937_64 RndEngine;
     uniform_real_distribution<double> DistUnif01;
     uniform_int_distribution<>DistUnifSite;
     uniform_int_distribution<>DistUnifState;

public:
   RandomNumberGenerator(const int seed_spec,const int n_spins_spec=0):seed_init(seed_spec),
   RndEngine(seed_init),DistUnif01(0.0,1.0),
   DistUnifSite(0,n_spins_spec-1),
   DistUnifState(0,static_cast<long long>(pow(2,n_spins_spec)-1)){};
   ~RandomNumberGenerator(){};

   double gen_rand01(){return DistUnif01(RndEngine);};//用于判定要不要翻转
   int gen_rand_site(){return DistUnifSite(RndEngine);};//用于选取翻转的位置
   long long gen_rand_state(){return DistUnifState(RndEngine);};//用于选取翻转的初态

   string save() const{
      stringstream save_data;
      save_data<<RndEngine;
      return save_data.str();
   }

   void load(string save_data){
      stringstream restore(save_data);
      restore>>RndEngine;
   }
};

class MC{
public:
   vector<int> mcs_samples;
   vector<double> hs;
   vector<double> ms; 
   vector<double> m2s;
   Isingsystem_square lattice;
   MC(int n,double T_spec,double Tmax_spec,int the_spec,int total_spec):
   lattice({n,n}),T(T_spec),Tmax(Tmax_spec),mcs_thermalization(the_spec),mcs_total(total_spec),randnum(123456890,n*n),
   mcs_samples(),hs(),ms(),m2s()
   {
     beta=1/T;
     //randomize_spins();
     lattice.set_state_by_code(1);
   }
     
   virtual ~MC(){};

   void run_MC(bool print_output=true){//输入print_out选择是否输出
   T=0;beta=1/T;
    while (T<=Tmax)
    {
        h_aver=0;
        m_aver=0;
        run_MC_given_T();
        if(print_output)
        {cout<<T<<" "<<h_aver<<" "<<m_aver<<" "<<m2_aver<<endl;};
        T=T+dT;
        beta=1/T;
    } 
   };

   void run_MC_given_T()
   {
      randomize_spins();//可以不要
      mcs_current=1;
      thermalize_given_T();
      mcs_sample=0;
      sampling_given_T();
   };

private:
   double T,Tmax,dT=0.05,beta,h_aver=0,m_aver=0,m2_aver=0;
   int mcs_current;
   int mcs_thermalization;
   int mcs_total;
   const int dt=10;//采样间隔
   int mcs_sample=0;
   RandomNumberGenerator randnum;

   
//随机初始状态
   void randomize_spins()
   {
    lattice.set_state_by_code(randnum.gen_rand_state());
   };

//随机翻转坐标
   void update_spins_Metropolis_randoms()
   {
      int flip_site;
      flip_site=randnum.gen_rand_site();
      double sz_i=lattice._sz(flip_site);
      //double Ei=lattice.energy_eval();
      lattice.flip_spin(flip_site);
      double sz_f=lattice._sz(flip_site);
      //double Ef=lattice.energy_eval();
      double dE=0;
      double p;
      double dE2;
      for(int i=0;i<=3;i++)
      {
         dE=dE+(sz_f-sz_i)*lattice._sz(lattice._spin_NN(flip_site,i));
      }
      dE=dE*lattice._J();
      //dE2=Ef-Ei;
      if(exp(-beta*dE)<1){
         p=randnum.gen_rand01();
         if(p>exp(-beta*dE)){lattice.flip_spin(flip_site);};
      }
      //cout<<p<<" "<<dE<<" "<<dE2<<endl;
   }

   void update_estimators()
   {
      
      if(mcs_current%dt==0) 
      {
         double h=lattice.energy_eval()/lattice._n_spins();
         double m=lattice.momentum_eval()/lattice._n_spins();
         double m2=m*m;
         mcs_sample++;

         /*ofstream outfile("sweep.txt", std::ios::app);
         outfile << mcs_sample << " " << h << " " << m <<" "<<m2<<endl;
         outfile.close();*/
         //cout<<mcs_current<<" "<<h<<" "<<m<<endl;
         mcs_samples.push_back(mcs_sample);
         ms.push_back(m);
         hs.push_back(h);
         m2s.push_back(m2);

         h_aver=(h_aver*(mcs_sample-1)+h)/(mcs_sample);
         m_aver=(m_aver*(mcs_sample-1)+m)/(mcs_sample);   
         m2_aver=(m2_aver*(mcs_sample-1)+m2)/(mcs_sample);
         
      }
      
   }

   void thermalize_given_T()
   {
      while(mcs_current<=mcs_thermalization){
         update_spins_Metropolis_randoms();
         mcs_current++;
      }
   };

   void sampling_given_T()
   {
      while(mcs_current<=mcs_total){
         update_spins_Metropolis_randoms();
         update_estimators();
         mcs_current++;
      }
   };
};

class observe_sta:public MC{
private:
    int bin_size;//bin大小
    int bin_gap;//bin间隔
    int bin_start;
    int n_bins;//bin个数
    vector<double> aver_bin_h;
    double aver_h;
    double error_h;

    vector<double> aver_bin_m;
    double aver_m;
    double error_m;

    vector<double> aver_bin_m2;
    double aver_m2;
    double error_m2;
public:
    observe_sta(int n,double T_spec,double Tmax_spec,int the_spec,int total_spec,
    int bin_start_spec,int bin_size_spec,int bin_gap_spec,int n_bins_spec)
    :MC(n,T_spec,Tmax_spec,the_spec,total_spec),
    bin_start(bin_start_spec),bin_size(bin_size_spec),bin_gap(bin_gap_spec),n_bins(n_bins_spec),
    aver_bin_h(),aver_bin_m(),aver_bin_m2()
    {aver_h=0;error_h=0;aver_m=0;error_m=0;aver_m2=0;error_m2=0;};
    void statistic(){
        for(int i=0;i<n_bins;i++){
            double h_bin=0,m_bin=0,m2_bin=0;
            for(int j=bin_start+i*(bin_size+bin_gap);j<bin_start+i*(bin_size+bin_gap)+bin_size;j++)
            {
                h_bin=h_bin+hs[j];
                m_bin=m_bin+ms[j];
                m2_bin=m2_bin+m2s[j];
            }
            h_bin=h_bin/bin_size;
            m_bin=m_bin/bin_size;
            m2_bin=m2_bin/bin_size;
            cout<<i<<" "<<h_bin<<" "<<m_bin<<" "<<m2_bin<<endl;
            aver_bin_h.push_back(h_bin);
            aver_bin_m.push_back(m_bin);
            aver_bin_m2.push_back(m2_bin);
        }
    }
};
//3 2 8 0 1000000 20000 800 200 80
#endif

