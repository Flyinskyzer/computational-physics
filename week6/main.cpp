#include "ising.hpp"

int main()
{
   int n;
   long long rep_state;
   cout<<"please enter the size of the one dimention ising model"<<endl;
   cin>>n;
   Isingsystem oneDspin(n);
   cout<<"please enter the number that decide the spin state"<<endl;
   cin>>rep_state;
   oneDspin.set_state_by_code(rep_state);
   oneDspin.display();
   cout<<"the energy is "<<oneDspin.eval_energy_1D()<<endl;
   cout<<"the magnetic moment is "<<oneDspin.eval_mz_1D()<<endl;
   return 0;
}
