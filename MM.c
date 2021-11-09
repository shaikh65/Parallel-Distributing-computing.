#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <omp.h>
#define NUM_THREADS 3

double generate_random(int *seed);
void matrix_multiplication(int l, int m, int n);
int main(int argc, char *argv[]);



int main(int argc, char *argv[])
{
 printf("\n\n-----------------------------------------PDC Project-------------------------\n");
 printf("\n----------------------Matrix multiplication in-parallel using OpenMP---------\n");
    int id,l,m,n;
    printf("\n");
    printf(" Number of processors available = %d\n",omp_get_num_procs());
    printf(" Max number of threads = %d\n",omp_get_max_threads());
    //l=m=n= 4;
    matrix_multiplication(100, 100, 100);
    return 0;
}


	

void matrix_multiplication(int l, int m, int n){
    double *A,*B,*C;
    int i,j,k,ops,seed;
    double rate;
    double time_begin;
    double funcTime;
    double time_elapsed;
    double time_stop;

    A = (double *)malloc(l * n * sizeof(double));
    B = (double *)malloc(l * m * sizeof(double));
    C = (double *)malloc(m * n * sizeof(double));
    

    seed = 123456789;
    for (k = 0; k < l * m; k++){
    
        B[k] = generate_random(&seed);
    }
    for (k = 0; k < m * n; k++){
    
        C[k] = generate_random(&seed);
    }


    funcTime = omp_get_wtime();
    time_begin = omp_get_wtime();
    omp_set_num_threads(2);
    int nthreads, id;
    
#pragma omp parallel shared(A, B, C, l, m, n) private(i, j, k)
{
#pragma omp for
    for (j = 0; j < n; j++)
    {
        for (i = 0; i < l; i++)
        {
            A[i + j * l] = 0.0;
            for (k = 0; k < m; k++)
            {
                A[i + j * l] = A[i + j * l] + B[i + k * l] * C[k + j * m];
            }
        }
    }
    
    
     // for (j = 0; j < n; j++)
   // {
   // printf("%f",A[j]);
  //  }

}
    time_stop = omp_get_wtime();

 
    ops = l * n * (2 * m);
    time_elapsed = time_stop - time_begin -(time_begin - funcTime);
    rate = (double)(ops) / time_elapsed / 1000000.0;
    
    printf("\n");
    printf("Matrix multiplication timing:\n");
    printf(" A = B * C\n");
    
    printf(" L = %d\n", l);
    printf(" M = %d\n", m);
    printf(" N = %d\n", n);
    
    printf(" floating point operations =: %d\n", ops);  
 /*Floating point OPS is a measure of computer performance, useful in fields of scientific computations that require floating-point calculations. 
 For such cases it is a more accurate measure than measuring instructions per second.*/
    printf(" Time of execution = %f\n", time_elapsed);
    printf(" Rate = FLOPS / dT = %f\n", rate);
    free(A);
    free(B);
    free(C);
    return;
}

double generate_random(int *seed){
    int k;
    double r;
    k = *seed / 127773;
    *seed = 16807 * (*seed - k * 127773) - k * 2836;
    if (*seed < 0){
        *seed = *seed + 2147483647;
    }
    r = (double)(*seed) * 4.656612875E-10;
    return r;
}


