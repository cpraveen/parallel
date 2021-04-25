#include <iostream>
#include <omp.h>

using namespace std;

int main()
{
   int n = 100;
   double a[n], b[n], c[n];

   for(int i=0; i<n; ++i)
   {
      a[i] = i;
      b[i] = 2*i;
   }

   #pragma omp parallel
   for(int i=0; i<n; ++i)
   {
      c[i] = a[i] + b[i];
   }

   for(int i=0; i<n; ++i)
      cout << a[i] << " " << b[i] << " " << c[i] << endl;

   return 0;
}
