#ifndef Ising
#define Ising
#include<iostream>
#include<vector>
#include<cmath>
using namespace std;

class Isingspin{
private:
   int sz;

public:
   Isingspin():sz(-1){};
   ~Isingspin(){};

   int _sz() const {return sz;};
   void set_up() {sz=1;};
   void set_dw() {sz=-1;};
   void set_sz(int sz_spec){
    if(sz_spec==1||sz_spec==-1)
     { sz=sz_spec;}
    else
    cout<<"please enter the correct sz";}
   void flip(){sz=sz*(-1);};
};

class Isingsystem{

private:
   const double J;
   int n_spins;
   const long long maxrep_state;
   vector<Isingspin> spin;

public:
   Isingsystem(int n_spin_spec) : J(-1.0),n_spins(n_spin_spec),maxrep_state(static_cast<long long>(pow(2,n_spins)-1))
   {
      spin.resize(n_spins);
   };
   virtual ~Isingsystem(){};

   double _J() const { return J;};
   int _n_spins() const { return n_spins;};
   long long _maxrep_state() const { return maxrep_state;};

   int _sz(const int site_idx) const { return spin[site_idx]._sz();};
   void set_up_spin(const int site_idx) { spin[site_idx].set_up();};
   void set_dw_spin(const int site_idx) { spin[site_idx].set_dw();};
   void set_spin(const int site_idx,int s_spec) { spin[site_idx].set_sz(s_spec);};
   void flip_spin(const int site_idx){spin[site_idx].flip();};
   vector<Isingspin> _spin(){return spin;};

   void set_state_by_code(long long rep_state) {
      int i=0;
      while(rep_state>0&& i < n_spins){
         if(rep_state%2==1)
         {
            set_up_spin(i);
         }
         else{
            set_dw_spin(i);
         }
         i++;
         rep_state=rep_state/2;
      }
   };

   void display() const{
      for(int i=0;i<n_spins;i++){
         cout<<spin[i]._sz()<<'\t';
      }
   }

   double eval_mz_1D() const {
      double M=0;
      for(int i=0;i<n_spins;i++)
         M=M+spin[i]._sz();
      M=M/n_spins;
      return M;
   };

   double eval_energy_1D() const {
      double E=0.0;
      for(int i=0;i<n_spins-1;i++)
         E=E+0.5*J*spin[i]._sz()*spin[i+1]._sz();
      E=E+0.5*J*spin[n_spins-1]._sz()*spin[0]._sz();
      return E;
   };
};
#endif
