#include"isingsquare.h"
int main(){
    Isingsystem_square spin2D({6,6});
    
    spin2D.set_state_by_code(pow(2,36)-1);
    spin2D.display();
    cout<<spin2D.energy_eval()<<endl;
    cout<<spin2D.momentum_eval()<<endl;
    
   /*cout<<spin2D.latice_coordinate(3)[0]<<'\t'<<spin2D.latice_coordinate(3)[1]<<endl;
   cout<<spin2D.NN(32,0)<<'\t'<<spin2D.NN(32,1)<<'\t'<<spin2D.NN(32,2)<<'\t'<<spin2D.NN(32,3)<<endl;
    return 0;*/
}


