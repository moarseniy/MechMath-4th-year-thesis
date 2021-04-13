#ifndef CUDA_LIB_
#define CUDA_LIB_



void TestCudaSolve(int *h_csrRowPtrA, int *h_csrColIndA, double *h_csrValA, int n, int nnz, double *h_b, double *h_x);
void Prepare_CSR(int *rows, int *ptr, int nnz, int n);
void LU_GPU_SOLVE(int *h_A_RowIndices, int *h_A_ColIndices, double *h_A, int n, int nnz, double *h_x, double *result);

void callCudaKernel();
void callCGM_GPU(int *x, int *y, double *data, double *b, double *x_k, int n, int sparse_size);



#endif
