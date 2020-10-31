#include "../include/for_you_to_do.h"

int get_block_size(){
    //return the block size you'd like to use 
    /*add your code here */
    return 129;
  
}

/**
 * 
 * this function computes LU factorization
 * for a square matrix
 * 
 * syntax 
 *  
 *  input : 
 *      A     n by n , square matrix
 *      ipiv  1 by n , vector
 *      n            , length of vector / size of matrix
 *  
 *  output :
 *      return -1 : if the matrix A is singular (max pivot == 0)
 *      return  0 : return normally 
 * 
 **/


int mydgetrf(double *A, int *ipiv, int n) 
{
    /* add your code here */
    // int i, maxind;
    // double max;
    // for(i = 0; i < n; i++){
    //     maxind = i;
    //     max = fabs(A[i*n + i]);
    //     int t;
    //     for(t = i+1; t < n; t++){
    //         if(fabs(A[t*n + i]) > max){
    //             maxind = t;
    //             max = fabs(A[t*n + i]);
    //         }
    //     }

    //     if(max == 0){
    //         printf("LU factoration failed: coefficient matrix is singular");
    //         return -1;
    //     }
    //     else{
    //         if(maxind != i){
    //             // save pivoting information
    //             int temps= ipiv[i];
    //             ipiv[i] = ipiv[maxind];
    //             ipiv[maxind] = temps;

    //             //swap row for matrix method 1
    //             // int j;
    //             // for(j = 0; j < n; j++){
    //             //     double k;
    //             //     k = A[i * n + j];
    //             //     A[i * n + j] = A[maxind * n + j];
    //             //     A[maxind * n + j] = k;
    //             // }

    //             //swap row method 2 -- need to test which one is faster
    //             double trow[n];
    //             memcpy(trow, A + i * n, n*sizeof(double));
    //             memcpy(A + i * n, A + maxind * n, n*sizeof(double));
    //             memcpy(A + maxind * n, trow, n*sizeof(double));
    //         }

    //     }

    //     int j;
    //     for(j = i + 1; j <n;j++){
    //         A[j*n + i] = A[j*n + i] / A[i*n + i];
    //         int k;
    //         for(k =  i + 1; k < n; k++){
    //             A[j*n + k] = A[j*n + k] - A[j*n + i] * A[i*n + k];
    //         }
    //     }
    // }

    // return 0;
}

/**
 * 
 * this function computes triangular matrix - vector solver
 * for a square matrix . according to lecture slides, this
 * function computes forward AND backward subtitution in the
 * same function.
 * 
 * syntax 
 *  
 *  input :
 *      UPLO  'L' or 'U' , denotes whether input matrix is upper
 *                         lower triangular . ( forward / backward
 *                         substitution )
 * 
 *      A     n by n     , square matrix
 * 
 *      B     1 by n     , vector
 * 
 *      ipiv  1 by n     , vector , denotes interchanged index due
 *                                  to pivoting by mydgetrf()
 * 
 *      n                , length of vector / size of matrix
 *  
 *  output :
 *      none
 * 
 **/
void mydtrsv(char UPLO, double *A, double *B, int n, int *ipiv)
{
    int i,j;
    double sum;
    /* add your code here */
    if(UPLO == 'L'){
        double y[n];
        
        y[0] = B[ipiv[0]];
        for(i = 1; i < n; i++ ){
           sum = 0;
           for(j = 0; j < i; j++ ){
               sum += y[j] * A[i*n + j];
           } 
           y[i] = B[ipiv[i]] - sum;
        }

        for(i = 0; i < n; i++){
            B[i] = y[i];
        }
    }

    else if(UPLO == 'U'){
        double x[n];
        int i,j;
        double sum;

        x[n-1] = B[n-1] / A[(n-1) * n + n - 1];

        for(i = n - 2; i >= 0; i-- ){
            sum = 0;
            for(j = i + 1; j < n; j++){
                sum += x[j] * A[i*n + j];
            }
            x[i] = (B[i] - sum) / A[i*n + i];
        }
        for(i = 0; i < n; i++){
            B[i] = x[i];
        }
    }
    return;
}

/**
 * 
 * Same function as what you used in lab1, cache_part4.c : optimal( ... ).
 * 
 **/
void mydgemm(double *A, double *B, double *C, int n, int i, int j, int k, int b)
{
    int ic,jc,kc,i1,j1,k1;
    register int m;
    int block_size = 3;

    for(i1 = i; i1 < n; i1 += b){
        if(i1 >= n){
            printf("i1 error\n");
            return;
        }
        // printf("i = %i", i1);
        for(j1 = j;j1 < n; j1 += b){
            // printf("j = %i", j1);
            if(j1 >= n){
                printf("j1 error\n");
                return;
            }
            for(k1 = k; k1 < k1 + b; k1 += b){
                if(k1 >= n){
                    printf("k1 error\n");
                    return;
                }

                for (ic=i1; ic<(ic + b); ic+=block_size){
                    if(ic >= n){
                        printf("ic error\n");
                        return;
                    }
                    for (jc=j1; jc<(jc + b); jc+=block_size) {

                        if(jc >= n){
                            printf("jc error\n");
                            return;
                        }

                        //9 register used for matrix C
                        register double c00 = C[ic * n + jc];
                        register double c01 = C[ic * n + (jc + 1)];
                        register double c02 = C[ic * n + (jc + 2)];

                        register double c10 = C[(ic + 1) * n + jc];
                        register double c11 = C[(ic + 1) * n + (jc + 1)];
                        register double c12 = C[(ic + 1) * n + (jc + 2)];

                        register double c20 = C[(ic + 2) * n + jc];
                        register double c21 = C[(ic + 2) * n + (jc + 1)];
                        register double c22 = C[(ic + 2) * n + (jc + 2)];

                        // double test = c00 + c01 + c02 + c10 + c11 + c12 + c20 + c21 + c22;

                        // printf("test = %f", test);


                        //6 registers INIT for A and B matrix
                        register double a00;
                        register double a10;
                        register double a20;
                        // register double a30;

                        register double b00;
                        register double b01;
                        register double b02;
                        // register double b03;

                        for (kc=k1; kc<(kc + b); kc+=block_size){
                            if(kc >= n){
                                printf("kc error\n");
                                printf("%i\n",kc); 
                                printf("%i\n",k1);
                                return;
                            }
                            for(m = 0; m < block_size; m++){
                                if(m >= n){
                                    printf("m error\n");
                                    return;
                                }
                                a00 = A[ic * n + kc + m];
                                a10 = A[(ic + 1)*n + kc + m];
                                a20 = A[(ic + 2)*n + kc + m];

                                b00 = B[(kc + m) * n + (jc)];
                                b01 = B[(kc + m) * n + (jc + 1)];
                                b02 = B[(kc + m) * n + (jc + 2)];

                                //Start doing the computing process
                                c00 -= a00 * b00;
                                c01 -= a00 * b01;
                                c02 -= a00 * b02;
                                c10 -= a10 * b00;
                                c11 -= a10 * b01;
                                c12 -= a10 * b02;
                                c20 -= a20 * b00;
                                c21 -= a20 * b01;
                                c22 -= a20 * b02;
                            }
                            

                        }
                        //Write back the value to matrix C
                        C[ic * n + jc] = c00;
                        C[ic * n + (jc + 1)] = c01;
                        C[ic * n + (jc + 2)] = c02;

                        C[(ic + 1) * n + jc] = c10;
                        C[(ic + 1) * n + (jc + 1)] = c11;
                        C[(ic + 1) * n + (jc + 2)] = c12;

                        C[(ic + 2) * n + jc] = c20;
                        C[(ic + 2) * n + (jc + 1)] = c21;
                        C[(ic + 2) * n + (jc + 2)] = c22;

                    }
                }

            }
        }
    }
    return;
}

/**
 * 
 * this function computes triangular matrix - vector solver
 * for a square matrix using block gepp introduced in course
 * lecture .
 * 
 * just implement the block algorithm you learned in class.
 * 
 * syntax 
 *  
 *  input :
 *      
 * 
 *      A     n by n     , square matrix
 * 
 *      B     1 by n     , vector
 * 
 *      ipiv  1 by n     , vector , denotes interchanged index due
 *                                  to pivoting by mydgetrf()
 * 
 *      n                , length of vector / size of matrix
 *  
 *  output :
 *      return -1 : if the matrix A is singular (max pivot == 0)
 *      return  0 : return normally 
 * 
 **/
int mydgetrf_block(double *A, int *ipiv, int n, int b) 
{

    /* add your code here */

    int i,j,k,ic,t, maxind;
    double max;

    for(ic = 0; ic <n;ic +=b){
        for(i = ic; i < ic+b ; i++){
            maxind = i;
            max = fabs(A[i*n + i]);
            for(t = i+1; t < n; t++){
                if(fabs(A[t*n + i]) > max){
                    maxind = t;
                    max = fabs(A[t*n + i]);
                }
            }

            if(max == 0){
                printf("LU factoration failed: coefficient matrix is singular");
                return -1;
            }
            else{
                if(maxind != i){
                    // save pivoting information
                    int temps= ipiv[i];
                    ipiv[i] = ipiv[maxind];
                    ipiv[maxind] = temps;

                    //swap row for matrix method 1
                    // int j;
                    // for(j = 0; j < n; j++){
                    //     double k;
                    //     k = A[i * n + j];
                    //     A[i * n + j] = A[maxind * n + j];
                    //     A[maxind * n + j] = k;
                    // }

                    //swap row method 2 -- need to test which one is faster
                    double trow[n];
                    memcpy(trow, A + i * n, n*sizeof(double));
                    memcpy(A + i * n, A + maxind * n, n*sizeof(double));
                    memcpy(A + maxind * n, trow, n*sizeof(double));
                }

            }

            for(j = i + 1; j <n;j++){
                A[j*n + i] = A[j*n + i] / A[i*n + i];

                //block version
                for(k = i + 1; k < i + b; k++){
                    A[j*n + k] = A[j*n + k] - A[j*n + i] * A[i*n + k];
                }

                //naive version - to test the top part of the code work
                // for(k = i + 1; k < n; k++){
                //     A[j*n + k] = A[j*n + k] - A[j*n + i] * A[i*n + k];
                // }
            }
        }

        //update A(ib:end, end+1:n), basically same method as before, use the value store in A(ib:n, ib:end)
        register double total;
        //end = ic + b
        for(i = ic; i < ic + b; i++){
            for(j= ic;j < n;j++){
                total = 0;
                for(k = ic; k < i; k++){
                    // A[i*n - j] -= A[i*n + k] * A[k*n + j];
                    total += A[i*n + k] * A[k*n + j];
                }
                A[i*n + j] -= total;
            }
        }

        // update A(end + 1: n , end + 1 : n)
        // end = ic + b
        mydgemm(A, A, A,n, ic + b, ic + b, ic, b);
    }

    

    return 0;
}

