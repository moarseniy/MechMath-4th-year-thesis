#include "init.h"

#include <stdio.h>

#include <cuda_runtime.h>
#include <cusparse_v2.h>



#include <iostream>
#include <cuda.h>
#include <cusolverSp.h>

using namespace std;

void TestCudaSolve(int *h_csrRowPtrA, int *h_csrColIndA, double *h_csrValA, int n, int nnz, double *h_b, double *h_x) {

    cusolverSpHandle_t handle;
    cusolverStatus_t status;
    cusparseStatus_t status2;

    status = cusolverSpCreate(&handle);


    cusparseMatDescr_t descr;
    status2 = cusparseCreateMatDescr(&descr);



    double* d_csrValA, *d_b, *d_x;
    int* d_csrRowPtrA, *d_csrColIndA;
    cudaMalloc((void**)&d_csrValA, nnz * sizeof(double));
    cudaMalloc((void**)&d_b, n * sizeof(double));
    cudaMalloc((void**)&d_x, n * sizeof(double));
    cudaMalloc((void**)&d_csrRowPtrA, (n + 1) * sizeof(int));
    cudaMalloc((void**)&d_csrColIndA, nnz * sizeof(int));


    cudaMemcpy(d_csrValA, h_csrValA, nnz * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(d_csrRowPtrA, h_csrRowPtrA, (n + 1) * sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(d_csrColIndA, h_csrColIndA, nnz * sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(d_b, h_b, n * sizeof(double), cudaMemcpyHostToDevice);



    cout<<"start solving...\n";
    double tol = 1e-16;
    int reorder = 1;
    int singularity = 0;
    status = cusolverSpDcsrlsvqr(handle, n, nnz, descr, d_csrValA, d_csrRowPtrA, d_csrColIndA, d_b, tol,
                     reorder, d_x, &singularity);

    cout<<"end solving...\n";
    cudaMemcpy(h_x, d_x, n * sizeof(double), cudaMemcpyDeviceToHost);
    //cout<<"singularity = "<<singularity<<"\n";


    cudaFree(d_csrValA);
    cudaFree(d_csrRowPtrA);
    cudaFree(d_csrColIndA);
    cudaFree(d_b);
    cudaFree(d_x);
    cusolverSpDestroy(handle);


}

void SortCOO(int n, int nnz) {

}



void setUpDescriptor(cusparseMatDescr_t& descrA, cusparseMatrixType_t matrixType, cusparseIndexBase_t indexBase) {
    cusparseCreateMatDescr(&descrA);
    cusparseSetMatType(descrA, matrixType);
    cusparseSetMatIndexBase(descrA, indexBase);
}


void setUpDescriptorLU(cusparseMatDescr_t& descrLU, cusparseMatrixType_t matrixType, cusparseIndexBase_t indexBase, cusparseFillMode_t fillMode, cusparseDiagType_t diagType) {
    cusparseCreateMatDescr(&descrLU);
    cusparseSetMatType(descrLU, matrixType);
    cusparseSetMatIndexBase(descrLU, indexBase);
    cusparseSetMatFillMode(descrLU, fillMode);
    cusparseSetMatDiagType(descrLU, diagType);
}


void memoryQueryLU(csrilu02Info_t& info_A, csrsv2Info_t& info_L, csrsv2Info_t& info_U, cusparseHandle_t handle, const int N, const int nnz, cusparseMatDescr_t descrA, cusparseMatDescr_t descr_L,
    cusparseMatDescr_t descr_U, double* d_A, int* d_A_RowIndices, int* d_A_ColIndices, cusparseOperation_t matrixOperation, void** pBuffer) {

    cusparseCreateCsrilu02Info(&info_A);
    cusparseCreateCsrsv2Info(&info_L);
    cusparseCreateCsrsv2Info(&info_U);

    int pBufferSize_M, pBufferSize_L, pBufferSize_U;
    cusparseDcsrilu02_bufferSize(handle, N, nnz, descrA, d_A, d_A_RowIndices, d_A_ColIndices, info_A, &pBufferSize_M);
    cusparseDcsrsv2_bufferSize(handle, matrixOperation, N, nnz, descr_L, d_A, d_A_RowIndices, d_A_ColIndices, info_L, &pBufferSize_L);
    cusparseDcsrsv2_bufferSize(handle, matrixOperation, N, nnz, descr_U, d_A, d_A_RowIndices, d_A_ColIndices, info_U, &pBufferSize_U);

    int pBufferSize = max(pBufferSize_M, max(pBufferSize_L, pBufferSize_U));
    cudaMalloc((void**)pBuffer, pBufferSize);

}


// ANALYSIS FUNCTION FOR LU DECOMPOSITION
void analysisLUDecomposition(csrilu02Info_t& info_A, csrsv2Info_t& info_L, csrsv2Info_t& info_U, cusparseHandle_t handle, const int N, const int nnz, cusparseMatDescr_t descrA, cusparseMatDescr_t descr_L,
    cusparseMatDescr_t descr_U, double* d_A, int* d_A_RowIndices, int* d_A_ColIndices, cusparseOperation_t matrixOperation, cusparseSolvePolicy_t solvePolicy1, cusparseSolvePolicy_t solvePolicy2, void* pBuffer) {

    int structural_zero;

    cusparseDcsrilu02_analysis(handle, N, nnz, descrA, d_A, d_A_RowIndices, d_A_ColIndices, info_A, solvePolicy1, pBuffer);
    cusparseStatus_t status = cusparseXcsrilu02_zeroPivot(handle, info_A, &structural_zero);
    if (CUSPARSE_STATUS_ZERO_PIVOT == status) { printf("A(%d,%d) is missing\n", structural_zero, structural_zero); }

    cusparseDcsrsv2_analysis(handle, matrixOperation, N, nnz, descr_L, d_A, d_A_RowIndices, d_A_ColIndices, info_L, solvePolicy1, pBuffer);
    cusparseDcsrsv2_analysis(handle, matrixOperation, N, nnz, descr_U, d_A, d_A_RowIndices, d_A_ColIndices, info_U, solvePolicy2, pBuffer);

}


// COMPUTE LU DECOMPOSITION FOR SPARSE MATRICES

void computeSparseLU(csrilu02Info_t& info_A, cusparseHandle_t handle, const int N, const int nnz, cusparseMatDescr_t descrA, double* d_A, int* d_A_RowIndices,
    int* d_A_ColIndices, cusparseSolvePolicy_t solutionPolicy, void* pBuffer) {

    int numerical_zero;

    cusparseDcsrilu02(handle, N, nnz, descrA, d_A, d_A_RowIndices, d_A_ColIndices, info_A, solutionPolicy, pBuffer);
    cusparseStatus_t status = cusparseXcsrilu02_zeroPivot(handle, info_A, &numerical_zero);
    if (CUSPARSE_STATUS_ZERO_PIVOT == status) { printf("U(%d,%d) is zero\n", numerical_zero, numerical_zero); }

}


void Prepare_CSR(int *rows, int *ptr, int nnz, int n) {
    cusparseHandle_t    handle;
    cusparseCreate(&handle);

    int *d_rows;
    cudaMalloc(&d_rows, nnz * sizeof(int));
    cudaMemcpy(d_rows, rows, nnz * sizeof(int), cudaMemcpyHostToDevice);

    int *d_ptr;
    cudaMalloc(&d_ptr, (n + 1) * sizeof(int));
    cudaMemcpy(d_ptr, ptr, (n + 1) * sizeof(int), cudaMemcpyHostToDevice);

    cusparseXcoo2csr(handle, d_rows, nnz, n, d_ptr, CUSPARSE_INDEX_BASE_ZERO);

    cudaMemcpy(ptr, d_ptr, (n + 1) * sizeof(int), cudaMemcpyDeviceToHost);
}

void LU_GPU_SOLVE(int *h_A_RowIndices, int *h_A_ColIndices, double *h_A, int n, int nnz, double *h_x, double *result)
{

    cusparseHandle_t    handle;

    cusparseMatDescr_t  descrA = 0;
    cusparseMatDescr_t  descr_L = 0;
    cusparseMatDescr_t  descr_U = 0;

    csrilu02Info_t      info_A = 0;
    csrsv2Info_t        info_L = 0;
    csrsv2Info_t        info_U = 0;

    void* pBuffer = 0;

    cudaEvent_t start, stop;
    cusparseCreate(&handle);
    cudaEventCreate(&start);
    cudaEventCreate(&stop);
    float gpuTime = 0.0;

    cudaEventRecord(start, 0);


    const int Nrows = n;
    const int Ncols = n;
    const int N = Nrows;


    double* d_x;
    cudaMalloc(&d_x, Nrows * sizeof(double));
    cudaMemcpy(d_x, h_x, Nrows * sizeof(double), cudaMemcpyHostToDevice);



    setUpDescriptor(descrA, CUSPARSE_MATRIX_TYPE_GENERAL, CUSPARSE_INDEX_BASE_ONE);

//    for (int i = 0; i < nnz; i++) {
//        printf("%f ", h_A[i]);
//    }


    for (int i = 0; i < nnz; i++) {
        //h_A_ColIndices[i]++;
        //printf("%d ", h_A_ColIndices[i]);
    }


    double* d_A;
    cudaMalloc(&d_A, nnz * sizeof(*d_A));

    int* d_A_RowIndices;
    cudaMalloc(&d_A_RowIndices, (Nrows + 1) * sizeof(*d_A_RowIndices));

    int* d_A_ColIndices;
    cudaMalloc(&d_A_ColIndices, nnz * sizeof(*d_A_ColIndices));

    cudaMemcpy(d_A, h_A, nnz * sizeof(*h_A), cudaMemcpyHostToDevice);
    cudaMemcpy(d_A_RowIndices, h_A_RowIndices, (Nrows + 1) * sizeof(*h_A_RowIndices), cudaMemcpyHostToDevice);
    cudaMemcpy(d_A_ColIndices, h_A_ColIndices, nnz * sizeof(*h_A_ColIndices), cudaMemcpyHostToDevice);


    cudaMemcpy(h_A, d_A, nnz * sizeof(*h_A), cudaMemcpyDeviceToHost);

    setUpDescriptorLU(descr_L, CUSPARSE_MATRIX_TYPE_GENERAL, CUSPARSE_INDEX_BASE_ONE, CUSPARSE_FILL_MODE_LOWER, CUSPARSE_DIAG_TYPE_UNIT);
    setUpDescriptorLU(descr_U, CUSPARSE_MATRIX_TYPE_GENERAL, CUSPARSE_INDEX_BASE_ONE, CUSPARSE_FILL_MODE_UPPER, CUSPARSE_DIAG_TYPE_NON_UNIT);


    memoryQueryLU(info_A, info_L, info_U, handle, N, nnz, descrA, descr_L, descr_U, d_A, d_A_RowIndices, d_A_ColIndices, CUSPARSE_OPERATION_NON_TRANSPOSE, &pBuffer);


    analysisLUDecomposition(info_A, info_L, info_U, handle, N, nnz, descrA, descr_L, descr_U, d_A, d_A_RowIndices, d_A_ColIndices, CUSPARSE_OPERATION_NON_TRANSPOSE, CUSPARSE_SOLVE_POLICY_NO_LEVEL,CUSPARSE_SOLVE_POLICY_USE_LEVEL, pBuffer);


    computeSparseLU(info_A, handle, N, nnz, descrA, d_A, d_A_RowIndices, d_A_ColIndices, CUSPARSE_SOLVE_POLICY_NO_LEVEL, pBuffer);


    double* d_z;
    cudaMalloc(&d_z, N * sizeof(double));

    const double alpha = 1.;
    cusparseDcsrsv2_solve(handle, CUSPARSE_OPERATION_NON_TRANSPOSE, N, nnz, &alpha, descr_L, d_A, d_A_RowIndices, d_A_ColIndices, info_L, d_x, d_z, CUSPARSE_SOLVE_POLICY_NO_LEVEL, pBuffer);


    double* d_y;
    cudaMalloc(&d_y, N * sizeof(double));

    cusparseDcsrsv2_solve(handle, CUSPARSE_OPERATION_NON_TRANSPOSE, N, nnz, &alpha, descr_U, d_A, d_A_RowIndices, d_A_ColIndices, info_U, d_z, d_y, CUSPARSE_SOLVE_POLICY_USE_LEVEL, pBuffer);


    double* h_y = (double*)malloc(Ncols * sizeof(double));
    cudaMemcpy(h_x, d_y, N * sizeof(double), cudaMemcpyDeviceToHost);

    cudaEventRecord(stop, 0);
    cudaEventSynchronize(stop);
    cudaEventElapsedTime(&gpuTime, start, stop);
    printf("GPU time = %.4f \n", gpuTime);

    printf("\n\nFinal result\n");
    for (int i = 0; i < N; i++) {
        printf("%f ", h_x[i]);
    }
}




///////////////////////////////////////////////////////
__global__
void cgm_gpu(double *z_k, double *r_k, double *Az,
             int* x, int *y, double *data, double *b,
             double *x_k, const int n, int sparse_size, double *partialSum) {

    double mf = 0.0, alpha, beta, eps = 0.00001, Spz, Spr, Spr1;

    int tx = threadIdx.x;
    int i = tx + blockIdx.x * blockDim.x;

    if (i < n) {
        partialSum[tx] = b[i] * b[i];
    }

    int stride;
    for (stride = blockDim.x/2; stride > 0;  stride >>= 1) {
        __syncthreads();
        if (tx < stride) {
            partialSum[tx] += partialSum[tx + stride];
        }
    }

    if (tx == 0) {
        //b[blockIdx.x] = partialSum[tx];
        mf = partialSum[tx];
    }



    x_k[i] = 0.2;
    Az[i] = 0.0;

    if (i < sparse_size) {
        Az[x[i]] += data[i] * x_k[y[i]];
    }
    r_k[i] = b[i] - Az[i];
    z_k[i] = r_k[i];

    //do{
        Spz=0.0;
        Spr=0.0;
        Az[i] = 0.0;
        if (i < sparse_size) {
            Az[x[i]] += data[i] * z_k[y[i]];
        }
        //Spz
        if (i < n) {
            partialSum[tx] = Az[i] * z_k[i];
        }

        for (stride = blockDim.x/2; stride > 0;  stride >>= 1) {
            __syncthreads();
            if (tx < stride) {
                partialSum[tx] += partialSum[tx + stride];
            }
        }
        if (tx == 0) {
            //b[blockIdx.x] = partialSum[tx];
            Spz = partialSum[tx];
            printf("Spz=%f\n", Spz);
        }
        //Spr
        if (i < n) {
            partialSum[tx] = r_k[i] * r_k[i];
        }

        for (stride = blockDim.x/2; stride > 0;  stride >>= 1) {
            __syncthreads();
            if (tx < stride) {
                partialSum[tx] += partialSum[tx + stride];
            }
        }

        if (tx == 0) {
            //b[blockIdx.x] = partialSum[tx];
            Spr = partialSum[tx];
            printf("Spr=%f\n", Spr);
        }
        ////
        alpha = Spr / Spz;
        Spr1 = 0.0;
        x_k[i] += alpha * z_k[i];
        r_k[i] -= alpha * Az[i];
        //Spr1
        if (i < n) {
            partialSum[tx] = r_k[i] * r_k[i];
            printf("%f ", r_k[i]);
        }

        for (stride = blockDim.x/2; stride > 0;  stride >>= 1) {
            __syncthreads();
            if (tx < stride) {
                partialSum[tx] += partialSum[tx + stride];
            }
        }

        if (tx == 0) {
            //b[blockIdx.x] = partialSum[tx];
            Spr1 = partialSum[tx];
            printf("Spr1=%f\n", Spr1);
        }
        ////
        beta = Spr1 / Spr;
        z_k[i] = r_k[i] + beta * z_k[i];


   // } while (Spr1 / mf > eps * eps);

    if (i == 0)
        printf("GPU CGM SUCCESS\n");
}

void callCGM_GPU(int *x, int *y, double *data, double *b, double *x_k, int n, int sparse_size) {
    double *z_k, *r_k, *Az;
    double *d_z_k, *d_r_k, *d_Az, *d_data, *d_b, *d_x_k, *partialSum;
    int *d_x, *d_y;
    z_k = (double*)malloc(n * sizeof(double));
    r_k = (double*)malloc(n * sizeof(double));
    Az = (double*)malloc(n * sizeof(double));


    cudaMalloc(&d_z_k, n * sizeof(double));
    cudaMalloc(&d_r_k, n * sizeof(double));
    cudaMalloc(&d_Az, n * sizeof(double));
    cudaMalloc(&d_x, sparse_size * sizeof(int));
    cudaMalloc(&d_y, sparse_size * sizeof(int));
    cudaMalloc(&d_data, sparse_size * sizeof(double));
    cudaMalloc(&d_b, n * sizeof(double));
    cudaMalloc(&d_x_k, n * sizeof(double));
    cudaMalloc(&partialSum, n * sizeof(double));

    cudaMemcpy(d_z_k, z_k, n * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(d_r_k, r_k, n * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(d_Az, Az, n * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(d_x, x, sparse_size * sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(d_y, y, sparse_size * sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(d_data, data, sparse_size * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(d_b, b, n * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(d_x_k, x_k, n * sizeof(double), cudaMemcpyHostToDevice);

    cgm_gpu<<<1, n>>>(d_z_k, d_r_k, d_Az, d_x, d_y, d_data, d_b, d_x_k, n, sparse_size, partialSum);

    cudaMemcpy(x_k, d_x_k, n * sizeof(double), cudaMemcpyDeviceToHost);


//    cudaFree(d_x);
//    cudaFree(d_y);
//    free(x);
//    free(y);
}



__global__
void saxpy(int n, float a, float *x, float *y)
{
  int i = blockIdx.x*blockDim.x + threadIdx.x;
  if (i < n) y[i] = a*x[i] + y[i];
}

void callCudaKernel()
{
  int N = 1<<8;
  float *x, *y, *d_x, *d_y;
  x = (float*)malloc(N*sizeof(float));
  y = (float*)malloc(N*sizeof(float));

  cudaMalloc(&d_x, N*sizeof(float)); 
  cudaMalloc(&d_y, N*sizeof(float));

  for (int i = 0; i < N; i++) {
    x[i] = 1.0f;
    y[i] = 2.0f;
  }

  cudaMemcpy(d_x, x, N*sizeof(float), cudaMemcpyHostToDevice);
  cudaMemcpy(d_y, y, N*sizeof(float), cudaMemcpyHostToDevice);

  // Perform SAXPY on 1M elements
  saxpy<<<(N+255)/256, 256>>>(N, 2.0f, d_x, d_y);

  cudaMemcpy(y, d_y, N*sizeof(float), cudaMemcpyDeviceToHost);

  float maxError = 0.0f;
  for (int i = 0; i < N; i++) {
    maxError = max(maxError, abs(y[i]-4.0f));
    printf("%f ", y[i]);
  }
  printf("Max error: %f\n", maxError);

  cudaFree(d_x);
  cudaFree(d_y);
  free(x);
  free(y);
}

