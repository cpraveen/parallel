// Solve
//    -Laplace(u) = 1  in  (0,1)x(0,1)
//              u = 0  on  boundary
#include <iostream>
#include <cmath>
#include <omp.h>

#define uold(i,j) uold_[j + (i)*N]
#define unew(i,j) unew_[j + (i)*N]
#define rhs(i,j)  rhs_[j + (i)*N]

#define swap(a, b) { auto tmp_ptr = a; a = b; b = tmp_ptr;}

using namespace std;

int main()
{
   const int N = 500;
   const double eps = 1.0e-10; // convergence threshold
   const double h = 1.0 / (N - 1);
   double maxdelta = 2.0 * eps;
   double *rhs_, *uold_, *unew_;

   uold_ = new double[N*N];
   unew_ = new double[N*N];
   rhs_ = new double[N*N];

   #pragma omp parallel for collapse(2)
   for(int i=0; i<N; ++i)
      for(int j=0; j<N; ++j)
      {
         uold(i,j) = 0.0;
         unew(i,j) = 0.0;
         rhs(i,j) = 1.0;
      }

   int iter = 0;
   while(maxdelta > eps)
   {
      maxdelta = 0.0;
      #pragma omp parallel for collapse(2) reduction(max:maxdelta)
      for(int i=1; i<N-1; ++i)
         for(int j=1; j<N-1; ++j)
         {
            unew(i,j) = (  uold(i+1,j) + uold(i-1,j)
                          + uold(i,j+1) + uold(i,j-1)
                          + h*h*rhs(i,j) ) * 0.25;
            maxdelta = fmax(maxdelta, fabs(unew(i,j)-uold(i,j)));
         }
      swap(uold_, unew_);
      ++iter;
      cout << iter << "  " << maxdelta << endl;
   }

   return 0;
}
