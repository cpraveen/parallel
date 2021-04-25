#include <iostream>
#include <omp.h>

using namespace std;

void say_hello(int my_thread_num, int total_threads)
{
   cout << "Hello I am thread " << my_thread_num << " among " << total_threads
        << " threads\n";
}

int main()
{
   cout << "I am the master, and I am alone\n";

   #pragma omp parallel
   {
      say_hello(omp_get_thread_num(), omp_get_num_threads());
   }
   return 0;
}
