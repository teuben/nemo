//    g++ sum_tbb.cc -o sum_tbb -ltbb

#include <oneapi/tbb.h>   

int main(){
  int n = 1001;
  int sum = oneapi::tbb::parallel_reduce(
					 oneapi::tbb::blocked_range<int>(1,n), 0,
					 [](oneapi::tbb::blocked_range<int> const& r, int init) -> int {
					   for (int v = r.begin(); v != r.end(); v++) {
					     init += v;
					   }
					   return init;
					 },
					 [](int lhs, int rhs) -> int {
					   return lhs + rhs;
					 }
					 );

    printf("N=%d Sum: %d\n", n, sum);
    return 0;
}
