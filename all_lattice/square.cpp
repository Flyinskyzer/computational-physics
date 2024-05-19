#include"utils/observe.h"
#include"utils/isingsquare.h"

int main()
{
    int n,m,T;
    vector<double> EdivN;
    vector<double> EdivN2;
    vector<double> EdivN4;
    vector<double> CdivN;

    vector<double> MdivN;
    vector<double> MdivN2;
    vector<double> MdivN4;
    vector<double> U;

    vector<double>_T;

    cin>>n;//晶格边长
    cin>>m;//每个晶格有多少个自旋
    cin>>T;//最高温度

    Isingsystem_square square({n,n});
    square.calculate();
    vector<double> energy_eve_state=square._h();
    vector<double> monmentum_eve_state=square._m();

    observe H(energy_eve_state,energy_eve_state,T);
    H.aver_everytemper();
    H.aver_everylattice(n*n*m);
    EdivN=H._q_aver();
    EdivN2=H._q2_aver();
    EdivN4=H._q4_aver();
    CdivN=H._c_aver();
    _T=H._T();


    observe M(monmentum_eve_state,energy_eve_state,T);
    M.aver_everytemper();
    M.aver_everylattice(n*n*m);
    MdivN=M._q_aver();
    MdivN2=M._q2_aver();
    MdivN4=M._q4_aver();
    U=M._u_aver();

    for(int i=0;i<EdivN.size();i++){
        cout<<_T[i]<<" "<<EdivN[i]<<" "<<MdivN[i]<<" "<< MdivN2[i]<<endl;
    }
}