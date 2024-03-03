#include<iostream>
#include<vector>
using namespace std;
int main()
{
   int n;
   int a,c;
   int b=0;
   std::vector<int> v;
   v.resize(1000);
   v[0]=0;
   v[1]=1;
   cin>>n;
   for(int i=2;i<=n;i++)
   {
      v[i]=v[i-1]+v[i-2];
   }
   cout<<v[n]<<endl;
}