#ifndef OBSERVE
#define OBSERVE
#include<vector>
#include<cmath>
using namespace std;
class observe
{
private:
    vector<double> q;
    double Tup;//最高温度
    vector<double> T;//不同温度
    int n;//温度数
    int num;//状态数
    int k=1;
    vector<double> q_aver,q2_aver,q4_aver,c_aver,u_aver;//不同温度
public:
    observe(){};
    observe(vector<double>q_spec,double tup):q(q_spec),Tup(tup){
        for(double i=0.05;i<=Tup;i++){
            T.push_back(i);
        }
        n=T.size();
        num=static_cast<long long>(q.size());
    }

    void aver_everytemper()
    {
        for(int i=0;i<n;i++)
        {
            double Q=0;
            double Q2=0;
            double Q4=0;
            double C=0;
            double z=0;
            double U=0;

            for(long long j=0;j<num;j++)
            {
                z=z+exp(-(q[j]-q[0])/(k*T[i]));
                Q=Q+q[j]*exp(-(q[j]-q[0])/(k*T[i]));
                Q2=Q2+q[j]*q[j]*exp(-(q[j]-q[0])/(k*T[i]));
                Q4=Q4+q[j]*q[j]*q[j]*q[j]*exp(-(q[j]-q[0])/(k*T[i]));
            }
            Q=Q/z;
            Q2=Q2/z;
            Q4=Q4/z;
            C=(Q2-Q*Q)/(T[i]*T[i]);
            U=Q4/(Q2*Q2);


            q_aver.push_back(Q);
            q2_aver.push_back(Q2);
            q4_aver.push_back(Q4);
            c_aver.push_back(C);
            u_aver.push_back(U);
        }
    }

    void aver_everylattice(int n)//输入晶格数
    {
        for(int i=0;i<n;i++){
            q_aver[i]=q_aver[i]/n;
            c_aver[i]=c_aver[i]/n;
            q2_aver[i]=q2_aver[i]/(n*n);
            q4_aver[i]=q4_aver[i]/pow(n,4);
        }
    }

    vector<double> _q_aver(){return q_aver;};
    vector<double> _q2_aver(){return q2_aver;};
    vector<double> _q4_aver(){return q4_aver;};
    vector<double> _c_aver(){return c_aver;};
    vector<double> _u_aver(){return u_aver;};
    vector<double> _T(){return T;};
};
#endif