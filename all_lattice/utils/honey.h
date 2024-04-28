#ifndef HONEY
#define HONEY

#include <iostream>
#include <cmath>
#include <cassert>
#include <stdio.h>
#include <vector>

using namespace std;

class IsingSpin
{
    private :
    int sz;

    public:
    IsingSpin() {sz =1;};
    ~IsingSpin() {};

    int _sz()const {return sz;};
    void set_up() {sz=1;};
    void set_dw() {sz=-1;};
    void set_sz(int sz_spec) {
        assert(sz_spec ==1||sz_spec ==-1);
        sz=sz_spec;
    };
    void flip() {sz *=-1;};
};

class IsingSystem 
{
    private :
        const double J;
        const int n_spins;
        const long long maxrep_state ;
        std::vector <IsingSpin> spin; //生成以IsingSpin为元素的一个数组
    
    public :
        IsingSystem(const int n_spins_spec) : J(-1.0), n_spins(n_spins_spec),maxrep_state(static_cast<long long>(std::pow(2, n_spins)) - 1)
        /*构造函数，：后面是初始化列表，相当于初始化J是常量，然后将初始化输入的值 n_spins_spec 赋值给类中量n_spins ，
        static_cast是用于类型转换，把括号内的数转换为long long型，这里就是用十进制的整数表示他的状态
        这个整数转换为二进制后有一串01,0表示自选向下，1表示自选向上，所以最大也就算2^n-1*/
        {
            spin.resize(n_spins);//vector的重新调整数组大小的函数，就是定义了有几个点
        };
        virtual ~IsingSystem() {};//西沟函数？这个沙没有也可以，不知道了


    double _J() const { return J; };//const函数定义参数，表明函数只能进行读取操作，不能进行写操作，就是相当于通过函数访问类中的私有成员？
    int _n_spins_() const { return n_spins;};
    long long _maxrep_state() const { return maxrep_state; };

    int _sz(const int site_idx) const { return spin[site_idx]._sz(); }//函数参数输入第几个量，然后就是数组中的第i个IsingSpin类，这个类里面的_sz()函数就是返回该点的自选方向
    void set_up_spin(const int site_idx) { spin[site_idx].set_up(); };//设定第i个自选向上 1
    void set_dw_spin(const int site_idx) { spin[site_idx].set_dw(); };//设定第i个自选向下 -1
    void set_spin(const int site_idx, int s_spec)
        {spin[site_idx].set_sz(s_spec); };//类中的set_sz函数，就没明白过，好像就是手动设定自选方向，输入1或者-1,sz就算输入量
    void flip_spin(const int site_idx) { spin[site_idx].flip(); };//第i个点自旋方向翻转

    void set_state_by_code(long long rep_state) //输入0~2^n-1的整数表示该一维体系所处状态，接下来进行二进制转换
    {
        int i=0;
        while (rep_state> 0) 
        {
           int remainder = rep_state % 2;  // 计算余数
           set_spin(i,2*remainder-1);   // 转换出来的二进制数第n-i位，然后设定第i为状态
           rep_state/= 2;   // 更新十进制数
           i++;               
        }
    
    // 添加必要的前导零以达到期望的位数
        while (i < n_spins) 
        {
            set_spin(i,-1);
            i++;
        }
    };

    double eval_mz() const 
    {
        double mz=0;
        for(int i=0;i<n_spins;i++)
        {
            mz=mz+_sz(i);
        }
        return mz;
    };

    double eval_energy_1D() const 
    {
        double e=0;
        for(int i=1;i<n_spins;i++)
        {
            e += _sz(i)*_sz(i-1);
        }
        e += _sz(n_spins-1)*_sz(0);
        e *= J;
        return e;
    };
};
class isingspin_honey_unit//晶格单位，两个粒子，新定义类
{
    private :
    vector<IsingSpin> unitspin;//以spin为元素，定义单位晶格大状态

    public:
    isingspin_honey_unit() {unitspin.resize(2);};//初始化大小为2,因为蜂窝的单位里面就两个
    ~isingspin_honey_unit() {};

    int unit_sz_i(const int unit_idx)const {return unitspin[unit_idx]._sz();};
    int unit_mz() {return unitspin[0]._sz()+unitspin[1]._sz();}//单位里面的总磁场;
    int unit_e() {return -unitspin[0]._sz()*unitspin[1]._sz();}//只考虑该单元两个的相互作用的能量
    vector<int> unit_sz()const {return {unit_sz_i(0),unit_sz_i(1)};};//返回一个表示状态的二维数对
    void unit_set_spin(const int unit_idx,int s_spec) {unitspin[unit_idx].set_sz(s_spec);};
    void unit_set_sz_by_num(int unit_sz_num) {//通过数字给单元赋值，0,1,2,3四个数字四个状态，懒得给输入加限定了
        int i=0;
        while (unit_sz_num> 0) 
        {
           int remainder = unit_sz_num % 2;  // 计算余数
           unitspin[i].set_sz(2*remainder-1);   // 转换出来的二进制数第n-i位，然后设定第i为状态
           unit_sz_num/= 2;   // 更新十进制数
           i++;               
        }
    // 添加必要的前导零以达到期望的位数
        while (i < 2) 
        {
            unitspin[i].set_sz(-1);
            i++;
        }
    };
};


class isingspin_honey_onlattice : public isingspin_honey_unit//以单元为节点构造晶格，类比就是把这两个粒子看作之前二维情况的一个点
{
    private:
    vector<int> position;//位置还是类似的而为有序数对
    vector<int> NN;//NN也一样，但是注意，周边单元的描述数字改成用 0,1,2,3描述不对，也可以直接用连接的粒子状态0,1

    public:
    isingspin_honey_onlattice()
    {
        set_dim(1);
        NN={0};
    }
    ~isingspin_honey_onlattice(){};
    void set_dim(int dim){position.assign(dim,0);};
    vector<int> _position ()const { return position;}
    vector<int> _NN() const {return NN;}
    int _NN(const int bond_idx) const {return NN[bond_idx];}
    void set_NN(const int bond_idx,const int site_idx){NN[bond_idx]=site_idx;};//注意现在知道的是周边的代表状态的数字0,1,2,3,不能直接得道周边的单元的第0,1的自旋状态
};

class Isingsystem_honey
{
    protected:
    const double J;
    const int n_unitspins;
    const long long maxrep_state;
    vector<isingspin_honey_onlattice> spin;


    public:
    Isingsystem_honey(const int n_unitspins_spec):J(1.0),n_unitspins(n_unitspins_spec),
    maxrep_state(static_cast<long long>(pow(4, n_unitspins)) - 1)
    {spin.resize(n_unitspins);};
    virtual ~Isingsystem_honey() {};
    void set_dim(int dim){for (auto &each:spin) each.set_dim(dim);};//给每个单元都setdim告诉他是几个维度的
    vector<int> _spin_position(const int site_idx) const {return spin[site_idx]._position();};//设定单元位置
    int _spin_NN(const int site_idx,const int bond_idx) const {return spin[site_idx]._NN(bond_idx);};//返回iNN的状态数字

    int _n_unitspins_() const { return n_unitspins;};
    long long _maxrep_state() const { return maxrep_state; };

    vector<int> _sz(const int site_idx) const { return spin[site_idx].unit_sz(); }//函数参数输入第几个量，然后就是数组中的第i个IsingSpin类，这个类里面的_sz()函数就是返回该点的自选方向
    void set_spin(const int site_idx, const int unit_idx,int s_spec)
        {spin[site_idx].unit_set_spin(unit_idx,s_spec); };//类中的set_sz函数，就没明白过，好像就是手动设定自选方向，输入1或者-1,sz就算输入量
    void honey_unit_set_spin(const int site_idx, int state_num) {spin[site_idx].unit_set_sz_by_num(state_num);}//通过数字0,1,2,3给第i给单元确定状态

    void set_state_by_num(long long rep_state) //输入0~4^n-1的整数表示该体系所处状态，接下来进行进制转换
    {
        int i=0;
        while (rep_state> 0) 
        {
           int remainder = rep_state % 4;  // 计算余数
           honey_unit_set_spin(i,remainder);   // 转换出来的二进制数第n-i位，然后设定第i为状态
           rep_state/= 4;   // 更新十进制数
           i++;               
        }

          // 添加必要的前导零以达到期望的位数
        while (i < n_unitspins) 
        {
            honey_unit_set_spin(i,0);
            i++;
        }
    }
};

class ising_system_honey: public Isingsystem_honey
{
    private:
    vector<int> system_size;
    void setup_NN ()  {
        for (int site_idx=0;site_idx<n_unitspins; site_idx++)
            {
                vector <int> r=lattice_coordinate (site_idx);
                spin[site_idx].set_NN(0,site_index(shift_pos_x(r)));
                spin[site_idx].set_NN(1,site_index(shift_pos_y(r)));
                spin[site_idx].set_NN(2,site_index(shift_neg_x(r)));
                spin[site_idx].set_NN(3,site_index(shift_neg_y(r)));
            }
    }
    vector<double>h;
    vector<double>m;

    public:
    ising_system_honey(const vector<int> sys_size):Isingsystem_honey(calc_n_spins(sys_size)),system_size(sys_size)
    { 
        setup_NN();
    };

    int calc_n_spins(const vector<int>& s_size)//计算矩阵大小，就是每个位置的值乘起来
    {
        int n_s = 1;
        for (auto &each : s_size) n_s *= each;
        return n_s;
    }
    ~ising_system_honey(){};
 
    int site_index(vector <int> lattice_coordinate) const
    {
        int _site_idx=lattice_coordinate[0]; 
        int w=1;
        for(int i=1;i<static_cast<int>(lattice_coordinate.size());i++)
        {
            for(int j=0;j<i;j++)
            {w*=system_size[j];}
            _site_idx=_site_idx+w*lattice_coordinate[i];
            w=1;
        };
        return _site_idx;
    };

    vector <int> lattice_coordinate(int site_idx) const
    {
        vector <int> l_c(system_size.size(),0);
        int w=1;
        for(int i=static_cast<int>(system_size.size()-1);i>=0;i--)
        {
            for(int j=0;j<i;j++)
            {w*=system_size[j];}
            l_c[i]=site_idx/w;
            site_idx=site_idx % w;
            w=1;
        };
        return l_c;
    };//前面的函数可以给出n维的情况，后面不想想了，就直接按老师的来，因为按老师这样，好像平移函数为给不出n维的，不想写了

    vector<int> shift_pos_x(const vector<int> r_spec)const {
        vector<int> r(r_spec);
        r[0]=(r[0]+1) % system_size[0];
        return r;}
    vector<int> shift_neg_x(const vector<int> r_spec)const {
        vector<int> r(r_spec);
        r[0]=((r[0]-1) % system_size[0]<0) ? ((r[0]-1) % system_size[0]+system_size[0]):((r[0]-1) % system_size[0]);
        return r;}
    vector<int> shift_pos_y(const vector<int> r_spec)const {
        vector<int> r(r_spec);
        r[1]=(r[1]+1) % system_size[1];
        return r;}
    vector<int> shift_neg_y(const vector<int> r_spec)const {
        vector<int> r(r_spec);
        r[1]=((r[1]-1) % system_size[1]<0) ? ((r[1]-1) % system_size[1]+system_size[1]):((r[1]-1) % system_size[1]);
        return r;}

        int NN (const int site_idx,const int bond_idx) const{
            return spin[site_idx]._NN(bond_idx);
        }

        void set_state(vector <bool> state_num)//这个是方便弄老师那个pi生成的测试
        {
            for (int site_idx=0;site_idx<n_unitspins; site_idx++)
            {honey_unit_set_spin(site_idx,2*state_num[site_idx]-1);}
        }


        double eval_energy()
        {
            double energy=0;
            for(int i=0;i<n_unitspins;i++)
            {
                energy += -spin[i].unit_sz_i(1)*spin[spin[i]._NN(0)].unit_sz_i(0);
                energy += -spin[i].unit_sz_i(1)*spin[spin[i]._NN(1)].unit_sz_i(0);
                energy += -spin[i].unit_sz_i(0)*spin[spin[i]._NN(2)].unit_sz_i(1);
                energy += -spin[i].unit_sz_i(0)*spin[spin[i]._NN(3)].unit_sz_i(1);
                energy += 2*spin[i].unit_e();
            }
            energy = energy/2;
            energy *= J;
            return energy;
        }

        double eval_Mz()
        {
            double mz=0;
            for(int i=0;i<n_unitspins;i++)
            {
                mz=mz+spin[i].unit_mz();
            }
            return mz;
        }
    void calculate(){
    for(long long i=0;i<=maxrep_state;i++)
    {
        set_state_by_num(i);
        h.push_back(eval_energy());
        m.push_back(eval_Mz());
    }
   };
   vector<double> _h(){return h;};
   vector<double> _m(){return m;};
};
#endif