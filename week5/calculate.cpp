#include "isingsquare.h"
#include <fstream>

int main()
{
    int n = 5;
    const double k = 1;
    Isingsystem_square spin2D({n, n});
    vector<double> T;
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
    
    vector<double> h(pow(2,n*n));
    vector<double> m(pow(2,n*n));

    for(long long j = 0; j <= pow(2, n*n) - 1; j++){
        spin2D.set_state_by_code(j);
        h[j] = spin2D.energy_eval();
        m[j]= spin2D.momentum_eval();
    }
    
    spin2D.set_state_by_code(0);
    double h0=spin2D.energy_eval();

    for (int i = 0; i < 80; i++)
    {
        H[i] = 0;
        H2[i] = 0;
        C[i] = 0;
        M2[i] = 0;
        double z = 0;

        for (long long j = 0; j <= pow(2, n*n) - 1; j++)
        {
            /*spin2D.set_state_by_code(j);
            double h = spin2D.energy_eval();
            double m = spin2D.momentum_eval();
            z = z + exp(-(h-h0) / (k * T[i]));
            H[i] = H[i] + h * exp(-(h-h0) / (k * T[i]));
            H2[i] = H2[i] + h * h * exp(-(h-h0) / (k * T[i]));
            M2[i] = M2[i] + m * m * exp(-(h-h0) / (k * T[i]));*/
            z = z + exp(-(h[j]-h0) / (k * T[i]));
            H[i] = H[i] + h[j] * exp(-(h[j]-h0) / (k * T[i]));
            H2[i] = H2[i] + h[j] * h[j] * exp(-(h[j]-h0) / (k * T[i]));
            M2[i] = M2[i] + m[j] * m[j] * exp(-(h[j]-h0) / (k * T[i]));
        }
        H[i] = H[i] / z;
        M2[i] = M2[i] / z;
        M2[i] = M2[i] / pow(n, 4);
        H2[i] = H2[i] / z;
        C[i] = (H2[i] - H[i] * H[i]) / (T[i] * T[i]);
        C[i] = C[i] / pow(n, 2);
        H[i] = H[i] / pow(n, 2);
        /*std::cout<<H[i]<<''<<M2[i]<<''<< C[i]<<''<<T[i]<<''<<H2[i]<<endl;*/
    }
    std::ofstream outfile("ising/output5.txt");
    outfile << "Temperature" << std::endl;
    /*for (size_t i = 0; i < T.size(); ++i)
    {
        outfile << T[i] << '\t';
    }*/
    outfile << "M2" << std::endl;
    for (size_t i = 0; i < M2.size(); ++i)
    {
        outfile << M2[i] << ',';
    }
    cout<<endl;
    outfile << "C" << std::endl;
    for (size_t i = 0; i < C.size(); ++i)
    {
        outfile << C[i] << ',';
    }
    cout<<endl;
    /*outfile << "H2" << std::endl;
    for (size_t i = 0; i < H2.size(); ++i)
    {
        outfile << H2[i] << '\t';
    }*/
    outfile << "H" << std::endl;
    for (size_t i = 0; i < H.size(); ++i)
    {
        outfile << H[i] << ',';
    }
    outfile.close();
    return 0;
}