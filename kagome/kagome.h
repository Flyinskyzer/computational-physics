#ifndef Isingkagome
#define Isingkagome
#include<iostream>
#include<vector>
#include<cmath>

using namespace std;
//一个自旋
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
//三个自旋组成一个单元unit，每个单元的状态用0--7表示
class Isingspinkagome :public Isingspin{
protected:
    vector<Isingspin> unit;
    int unitstate=0;//0--7
public:
    Isingspinkagome(){
        unit.resize(3);
        for(int i=0;i<=2;i++)
        {
            unit[i].set_sz(-1);
        }
        setunit_by_code(0);
    };
    ~Isingspinkagome(){};
    void setunit_by_code(int unitnum)//设置单元状态
    {
        unitstate=unitnum;
        for(int i=0;i<=2;i++){
            unit[i].set_sz(-1);
        };
        int i = 0;
        while (unitnum > 0 && i <=2)
        {
            if (unitnum % 2 == 1)
            {
                unit[i].set_sz(1);
            }
            else
            {
                unit[i].set_sz(-1);
            }
            i++;
            unitnum=unitnum/2;
        }
    }

    void unit_set_up(int i)
    {
        unit[i].set_up();
    }

    void unit_set_dw(int i)
    {
        unit[i].set_dw();
    }

    void unit_flip(int i)
    {
        unit[i].flip();
    }
    int _unitstate() const {return unitstate;};
    //返回unit中的某个自旋状态
    int _unit_spinstate(int i){
        return unit[i]._sz();
    }
};
//给出每个单元的坐标position，以及1*6的向量NN用于存储与该单元有关的单元的一维坐标
class IsingSpinOnLattice :public Isingspinkagome{
    private:
    vector<int> position;
    vector<int> NN={0,0,0,0,0,0};

    public:
    IsingSpinOnLattice()
    {
        set_dim(1);
    }
    ~IsingSpinOnLattice(){};
    void set_dim(int dim){position.assign(dim,0);};
    vector<int> _position ()const { return position;}
    vector<int> _NN() const {return NN;}
    int _NN(const int bond_idx) const {return NN[bond_idx];}
    void set_NN(const int bond_idx,const int site_idx){NN[bond_idx]=site_idx;};
    
};
//整个系统，给出向量units用于储存单元unit;
class Isingsystem{
protected:
   const double J;
   int n_units;//单元个数
   const long long maxrep_state;//系统的最大状态数
   vector<IsingSpinOnLattice> units;
   vector<int> system_size={0,0};

public:
   Isingsystem(const vector<int>system_size_spec): 
   J(-1.0),
   system_size(system_size_spec),
   n_units(system_size_spec[0]*system_size_spec[1]),
   maxrep_state(static_cast<long long>(pow(8,n_units)-1))
   {
      Isingspinkagome();
      units.resize(n_units);
      int dim=2;
      for(auto& each:units)each.set_dim(dim);
      setup_NN();
   };
   virtual ~Isingsystem(){};

   vector<int> _units_position(const int site_idx) const{return units[site_idx]._position();};
   vector<int> _units_NN(const int site_idx) const{return units[site_idx]._NN();};
   int _units_NN(const int site_idx,const int bond_idx) const{return units[site_idx]._NN(bond_idx);};

   double _J() const{return J;};
   int _n_units() const{return n_units;};
   long long _maxrep_state()const{return maxrep_state;};
   //设置units中某个unit的状态
   void set_units_by_code(int site_idx,int unitnum){
    units[site_idx].setunit_by_code(unitnum);
   };
   //返回units中某个unit的状态
   int _state(const int site_idx) const {return units[site_idx]._unitstate();};
   //返回units中的某个unit中的某个三品的状态
   int _spin_state(const int site_idx,const int j) {return units[site_idx]._unit_spinstate(j);};
   //设置整个units的状态
   void set_unitstate_in_total(long long rep_state){
    int i=0;
    for(int j=0;j<=n_units-1;j++){
        set_units_by_code(j,0);
    }
    while(rep_state>=0&&i<n_units){
        if(rep_state%8==0)
        {
            set_units_by_code(i,0);
        }
        if(rep_state%8==1)
        {
            set_units_by_code(i,1);
        }
        if(rep_state%8==2)
        {
            set_units_by_code(i,2);
        }
        if(rep_state%8==3)
        {
            set_units_by_code(i,3);
        }
        if(rep_state%8==4)
        {
            set_units_by_code(i,4);
        }
        if(rep_state%8==5)
        {
            set_units_by_code(i,5);
        }
        if(rep_state%8==6)
        {
            set_units_by_code(i,6);
        }
        if(rep_state%8==7)
        {
            set_units_by_code(i,7);
        }
        i++;
    rep_state=rep_state/8;
    }
   }

   void display() {
    for(int i=0;i<n_units;i++){
        cout<<units[i]._unitstate()<<" ";
    }
    cout<<endl;
    cout<<_n_units()<<endl;
    for(int i=0;i<n_units;i++){
        for(int j=0;j<=2;j++)
        {
        cout<<units[i]._unit_spinstate(j)<<" ";
        }
    }
   }

   //unit的二维坐标返回一维位置
   int site_index(const vector<int>lattice_coordinate) const{
    int index;
    index=(lattice_coordinate[1])*system_size[0]+lattice_coordinate[0];
    return index;
   }
   //unit的一维位置返回二维坐标
   vector<int>lattice_coordinate(int site_index) const {
      vector<int>coordinate={0,0};
      coordinate[1]=site_index/system_size[0];
      coordinate[0]=site_index%system_size[0];
      
      return coordinate;
   };

   //寻找与某晶格有关联的六个晶格

  
   //1
   //右
  vector<int>shift_pos_x(const vector<int>r_spec) const{
      vector<int>r(r_spec);
      r[0]=(r[0]+1)%system_size[0];
      return r;
   };
   //上
   vector<int>shift_pos_y(const vector<int>r_spec) const{
      vector<int>r(r_spec);
      r[1]=(r[1]+1)%system_size[1];
      return r;
   };

   //2
   //左
   vector<int>shift_neg_x(const vector<int>r_spec) const{
      vector<int>r(r_spec);
      r[0]=(r[0]-1+system_size[0])%system_size[0];
      return r;
   };
   
   //左上
   vector<int>shift_neg_x_pos_y(const vector<int>r_spec) const{
      vector<int>r(r_spec);
      r[1]=(r[1]+1)%system_size[1];
      r[0]=(r[0]-1+system_size[0])%system_size[0];
      return r;
   };

   //0
   //下
    vector<int>shift_neg_y(const vector<int>r_spec) const{
      vector<int>r(r_spec);
      r[1]=(r[1]-1+system_size[1])%system_size[1];
      return r;
   };
   //右下
   vector<int>shift_pos_x_neg_y(const vector<int>r_spec)const{
      vector<int>r(r_spec);
      r[1]=(r[1]-1+system_size[1])%system_size[1];
      r[0]=(r[0]+1)%system_size[0];
      return r;
   };

   //建立网络
   void setup_NN(){
    for(int site_idx=0;site_idx<n_units;site_idx++){
         vector<int>r=lattice_coordinate(site_idx);
         units[site_idx].set_NN(0,site_index(shift_pos_x(r)));
         units[site_idx].set_NN(1,site_index(shift_pos_y(r)));
         units[site_idx].set_NN(2,site_index(shift_neg_x(r)));
         units[site_idx].set_NN(3,site_index(shift_neg_y(r)));
         units[site_idx].set_NN(4,site_index(shift_neg_x_pos_y(r)));
         units[site_idx].set_NN(5,site_index(shift_pos_x_neg_y(r)));
      }
   };

   int NN(const int site_idx,const int bond_idx) const{
    return units[site_idx]._NN(bond_idx);
   }

   //计算某个状态下的能量和角动量

   double momentum_eval(){
    double M=0;
    for(int i=0;i<n_units;i++){
        for(int j=0;j<=2;j++){
            M=M+units[i]._unit_spinstate(j);
        }
    }
    return M;
   };

   void display_NN(){
    for(int i=0;i<n_units;i++)
    {
        for(int j=0;j<6;j++)
        {
            cout<<_units_NN(i,j)<<" ";
        }
        cout<<endl;
    }
   };

   double energy_eval(){
    double H=0;
    for(int i=0;i<n_units;i++){
        H=H+units[i]._unit_spinstate(2)*units[NN(i,2)]._unit_spinstate(1)
           +units[i]._unit_spinstate(2)*units[NN(i,4)]._unit_spinstate(0)
           +units[i]._unit_spinstate(1)*units[NN(i,1)]._unit_spinstate(0)
           +units[i]._unit_spinstate(1)*units[NN(i,0)]._unit_spinstate(2)
           +units[i]._unit_spinstate(0)*units[NN(i,3)]._unit_spinstate(1)
           +units[i]._unit_spinstate(0)*units[NN(i,5)]._unit_spinstate(2)
           +units[i]._unit_spinstate(0)*units[i]._unit_spinstate(1)
           +units[i]._unit_spinstate(0)*units[i]._unit_spinstate(2)
           +units[i]._unit_spinstate(1)*units[i]._unit_spinstate(2)
           +units[i]._unit_spinstate(1)*units[i]._unit_spinstate(0)
           +units[i]._unit_spinstate(2)*units[i]._unit_spinstate(1)
           +units[i]._unit_spinstate(2)*units[i]._unit_spinstate(0)
           ;
    }
    H=H*J*0.5;
    return H;
   }
};
#endif