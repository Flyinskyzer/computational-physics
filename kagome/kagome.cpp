#include"kagome.h"
#include <fstream>

int main()
{
    int n;
    cin>>n;
    const double k=1;
    Isingsystem kagome({n,n});
    vector<double>T;
    for (double i = 0.05; i <=4.05; i += 0.05)
    {
        T.push_back(i);
    }
    vector<double> H;
    H.resize(80);
    vector<double> H2;
    H2.resize(80);
    vector<double> C;
    C.resize(80);
    vector<double> M2;
    M2.resize(80);
    vector<double> M4;
    M4.resize(80);
    vector<double> U;
    U.resize(80);

    vector<double> h(pow(8,n*n));
    vector<double> m(pow(8,n*n));

    for(long long j=0;j<=pow(8,n*n)-1;j++){
        kagome.set_unitstate_in_total(j);
        h[j]=kagome.energy_eval();
        m[j]=kagome.momentum_eval();
    }

    double h0=h[0];

    for(int i=0;i<80;i++)
    {
        H[i] = 0;
        H2[i] = 0;
        C[i] = 0;
        M2[i] = 0;
        M4[i]=0;
        U[i]=0;
        double z = 0;

        for(long long j=0;j<= pow(8, n*n) - 1; j++)
        {
            z = z + exp(-(h[j]-h0) / (k * T[i]));
            H[i] = H[i] + h[j] * exp(-(h[j]-h0) / (k * T[i]));
            H2[i] = H2[i] + h[j] * h[j] * exp(-(h[j]-h0) / (k * T[i]));
            M2[i] = M2[i] + m[j] * m[j] * exp(-(h[j]-h0) / (k * T[i]));
            M4[i]=M4[i]+pow(m[j],4)* exp(-(h[j]-h0) / (k * T[i]));
        }
        H[i] = H[i] / z;
        M2[i] = M2[i] / z;
        M4[i]=M4[i]/z;
        U[i]=M4[i]/(M2[i]*M2[i]);
        M2[i] = M2[i] /pow(3*pow(n, 2),2);
        H2[i] = H2[i] / z;
        C[i] = (H2[i] - H[i] * H[i]) / (T[i] * T[i]);
        C[i] = C[i] / (3*pow(n, 2));
        H[i] = H[i] / (3*pow(n, 2)) ;
        cout<<H[i]<<" "<<C[i]<<" "<<M2[i]<< " "<<U[i]<<" "<<endl;
    }
return 0;
}