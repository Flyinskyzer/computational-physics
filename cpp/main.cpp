#include <iostream>
#include <vector>
#include "ising.hpp"

using namespace std;
int main(int argc, const char *argv[])
{
    int n;
    long long spin_config;
    cout << "The number of particles:";
    cin >> n;
    cout << "The number of spin configuration:";
    cin >> spin_config;
    IsingSystem oneDspin(n);
    oneDspin.set_state_by_code(spin_config);

    cout << "The magnetization:" << oneDspin.eval_mz() << endl;
    cout << "The energy:" << oneDspin.eval_energy_1D() << endl;
    cin.get();
    return 0;
}