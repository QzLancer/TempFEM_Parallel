#include <stdio.h>
#include <omp.h>

int main() {
  int j;

  //#pragma omp parallel for
  //#pragma omp taskwait
  //#pragma omp barrier
  //#pragma omp for ordered
  //#pragma omp parallel for num_threads(40)

  omp_set_num_threads(40);
  #pragma omp parallel for

  for(j=1;j<41;j++){
     printf("j=%d.\n",j);
     if( j % 2 == 0)
         printf("j=%d.\n",j);
  }
  printf("OK!\n");
  return 0;
}
