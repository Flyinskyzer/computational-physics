#ifndef IsingSystem_hpp
#define IsingSystem_hpp

#include <iostream> //基本类
#include <cassert> //这里是断言函数
#include <cmath>
#include <vector>

//下面我们定义IsingSpin类，表示粒子的状态
class IsingSpin {
    private:
        int sz; /* +/-1 */ 

    public:
        IsingSpin() { sz = 1; };
        ~IsingSpin() {};

    int _sz() const { return sz; };//定义粒子状态函数
    void set_up() { sz = 1; };
    void set_dw() { sz = -1; };
    void set_sz(int sz_spec) {
    assert(sz_spec == 1 || sz_spec == -1);//保证sz_spec一定是+-1 //类似自检程序
    sz = sz_spec;
    };
    void flip() { sz *= -1; };//定义函数flip（改变正负）
};

//下面是关于磁矩和能量的计算
class IsingSystem
{
private:
    const double J; //耦合
    const int n_spins; //粒子数
    const long long maxrep_state;
    std::vector<IsingSpin> spin;//动态数组spin，表示n个粒子的状态

public:
    IsingSystem(const int n_spin_spec) : J(-1.0) /*J的初始值为-1*/, n_spins(n_spin_spec),
    maxrep_state(static_cast<long long> (std::pow(2, n_spins)) -1) 
    {spin.resize(n_spins);};//二进制转换：旋转配置数。此处将spin的大小调成n

    virtual ~IsingSystem() {};//虚函数


    //这一段应该是保留全局变量。让J，n，二进制不变
    double _J() const { return J; };//固定J的值，变量保留实现
    int _n_spins() const {return n_spins;}; //同上
    long long _maxrep_state() const {return maxrep_state;}; //这个和n相关，通上


    int _sz(const int site_idx/*第几个粒子*/ ) const {return spin[site_idx]._sz();};
    //函数_sz，用于获取特定位置site_idx处自旋的sz值。
    //参数site_idx表示自旋在spin向量中的索引位置。函数内部通过调用spin[site_idx]._sz()来获取该位置自旋的sz值

    void set_up_spin(const int site_idx) {spin[site_idx].set_up();};
    void set_dw_spin(const int site_idx) {spin[site_idx].set_dw();};

    void set_spin(const int site_idx/*第几个*/, int s_spec/*怎么转 +-1, 否则报错终止*/) //粒子状态设定
        { spin[site_idx].set_sz(s_spec);};

    void flip_spin(const int site_idx ) {spin[site_idx].flip();};//翻转 //！没找到在那里使用

    void set_state_by_code(long long rep_state) {
        int i = 0;

        //下面二进制转十进制
        while (rep_state > 0)
        {
            int remainder = rep_state % 2;
            set_spin(i,2*remainder - 1);
            rep_state /=2;
            i++;
        }

        while (i < n_spins)
        {
            set_spin(i,-1);//设定粒子状态
            i++;
        }
    };

    double eval_mz() const {
        int mz = 0;
        for (int i = 0; i < n_spins; i++)
        {
            mz = mz + _sz(i);
        }
        return mz;
    };

    double eval_energy_1D() const {
        double e=0;
        for(int i=1;i<n_spins;i++)
        {
            e += _sz(i)*_sz(i-1);
        }
        e += _sz(n_spins-1)*_sz(0);//环状
        e *= (-0.5*J);
        return e;
    };
};
#endif 