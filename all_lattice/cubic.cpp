#include"utils/observe.h"
#include"utils/cubic.h"

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

    cin>>n;//晶格边长
    cin>>T;//最高温度

    ising_system_3D cubic({n,n,n});
    cubic.calculate();
    vector<double> energy_eve_state=cubic._h();
    vector<double> monmentum_eve_state=cubic._m();

    observe H(energy_eve_state,energy_eve_state,T);
    H.aver_everytemper();
    H.aver_everylattice(n*n*n);
    EdivN=H._q_aver();
    EdivN2=H._q2_aver();
    EdivN4=H._q4_aver();
    CdivN=H._c_aver();

    observe M(monmentum_eve_state,energy_eve_state,T);
    M.aver_everytemper();
    M.aver_everylattice(n*n*n);
    MdivN=M._q_aver();
    MdivN2=M._q2_aver();
    MdivN4=M._q4_aver();
    U=M._u_aver();

    
    for(int i=0;i<EdivN.size();i++){
        cout<<EdivN[i]<<endl;
    }
}