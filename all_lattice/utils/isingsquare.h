#ifndef Isingsquqre
#define Isingsquare
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

class IsingSpinOnLattice :public Isingspin{
    private:
    vector<int> position;
    vector<int> NN={0,0,0,0};

    public:
    IsingSpinOnLattice()
    {
        set_dim(1);
        NN={-1};
    }
    ~IsingSpinOnLattice(){};
    void set_dim(int dim){position.assign(dim,0);};
    vector<int> _position ()const { return position;}
    vector<int> _NN() const {return NN;}
    int _NN(const int bond_idx) const {return NN[bond_idx];}
    void set_NN(const int bond_idx,const int site_idx){NN[bond_idx]=site_idx;};
};

class Isingsystem{

protected:
   const double J;
   int n_spins;
   const long long maxrep_state;
   vector<IsingSpinOnLattice> spin;

public:
   Isingsystem(int n_spin_spec) : J(-1.0),n_spins(n_spin_spec),maxrep_state(static_cast<long long>(pow(2,n_spins)-1))
   {
      spin.resize(n_spins);
      int dim=2;
      for(auto& each:spin)each.set_dim(dim);
   };
   virtual ~Isingsystem(){};

   //void set_dim(int dim=2){for(auto& each:spin)each.set_dim(dim);};//设置维度
   vector<int> _spin_position(const int site_idx) const {return spin[site_idx]._position();};
   vector<int> _spin_NN(const int site_idx) const {return spin[site_idx]._NN();};
   int _spin_NN(const int site_idx,const int bond_idx) const {return spin[site_idx]._NN(bond_idx);}

   double _J() const { return J;};
   int _n_spins() const { return n_spins;};
   long long _maxrep_state() const { return maxrep_state;};

   int _sz(const int site_idx) const { return spin[site_idx]._sz();};
   void set_up_spin(const int site_idx) { spin[site_idx].set_up();};
   void set_dw_spin(const int site_idx) { spin[site_idx].set_dw();};
   void set_spin(const int site_idx,int s_spec) { spin[site_idx].set_sz(s_spec);};
   void flip_spin(const int site_idx){spin[site_idx].flip();};

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
      cout<<endl;
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

class Isingsystem_square:public Isingsystem{
private:
   vector<int> system_size;
   void setup_NN(){
      for(int site_idx=0;site_idx<n_spins;site_idx++){
         vector<int>r=latice_coordinate(site_idx);
         spin[site_idx].set_NN(0,site_index(shift_pos_x(r)));
         spin[site_idx].set_NN(1,site_index(shift_pos_y(r)));
         spin[site_idx].set_NN(2,site_index(shift_neg_x(r)));
         spin[site_idx].set_NN(3,site_index(shift_neg_y(r)));
      }
   }
   //储存每个状态的能量和角动量
   vector<double>h;
   vector<double>m;
public:
   Isingsystem_square(const vector<int>system_size_spec):
     Isingsystem(system_size_spec[0]*system_size_spec[1]),system_size(system_size_spec){setup_NN();};
   ~Isingsystem_square(){};
   int site_index(const vector<int>latice_coordinate) const {
      int index;
      index=(latice_coordinate[1])*system_size[0]+latice_coordinate[0];
      return index;
   };
   vector<int>latice_coordinate(int site_index) const {
      vector<int>coordinate={0,0};
      coordinate[1]=site_index/system_size[0];
      coordinate[0]=site_index%system_size[0];
      
      return coordinate;
   };
   vector<int>shift_pos_x(const vector<int>r_spec) const{
      vector<int>r(r_spec);
      r[0]=(r[0]+1)%system_size[0];
      return r;
   };
   vector<int>shift_pos_y(const vector<int>r_spec) const{
      vector<int>r(r_spec);
      r[1]=(r[1]+1)%system_size[1];
      return r;
   };
   vector<int>shift_neg_x(const vector<int>r_spec) const{
      vector<int>r(r_spec);
      r[0]=(r[0]-1+system_size[0])%system_size[0];
      return r;
   };
    vector<int>shift_neg_y(const vector<int>r_spec) const{
      vector<int>r(r_spec);
      r[1]=(r[1]-1+system_size[1])%system_size[1];
      return r;
   };
   /*
   int NN(const int site_idx,const int bond_idx) const{
      vector<int>r=latice_coordinate(site_idx);
      switch(bond_idx){
         case 0:
         return site_index(shift_pos_x(r));
         break;

         case 1:
         return site_index(shift_pos_y(r));
         break;

         case 2:
         return site_index(shift_neg_x(r));
         break;

         case 3:
         return site_index(shift_neg_y(r));
         break;
      }

   }
   */
   int NN(const int site_idx,const int bond_idx) const{
      return spin[site_idx]._NN(bond_idx);
   }

   vector<IsingSpinOnLattice> _spin(){return spin;};

   double energy_eval(){
      double H=0;
      for(int i=0;i<n_spins;i++){
         for(int j=0;j<4;j++){
             H=H+spin[i]._sz()*spin[NN(i,j)]._sz();
         }
      }
      H=H*J*0.5;
      return H;
   }

   double momentum_eval(){
      double M=0;
      for(int i=0;i<n_spins;i++){
            M=M+spin[i]._sz();
      }
      return M;
   }
   void calculate(){
    for(long long i=0;i<=maxrep_state;i++)
    {
        set_state_by_code(i);
        h.push_back(energy_eval());
        m.push_back(momentum_eval());
    }
   }

   vector<double> _h(){return h;};
   vector<double> _m(){return m;};

}; 
#endif