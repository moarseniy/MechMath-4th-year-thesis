#ifndef CUDA_LIB_
#define CUDA_LIB_

#include <iostream>
//#include "femfunc.h"

void CudaSolve(int *h_csrRowPtrA, int *h_csrColIndA, float *h_csrValA, int n, int nnz, float *h_b, float *h_x);
void Prepare_CSR(int *rows, int *ptr, int nnz, int n);
void SortCOO(int *h_cooRows, int *h_cooCols, float *h_cooVals, int n, int nnz);

//void FiniteElementMethodCUDA(float *D, int *h_elements, int elementsCount, int *h_nodesX, int *h_nodesY, int *h_nodesZ, int nodesCount);
void FiniteElementMethodCUDA(float *h_D, int *h_elements,
                             int *h_elements0, int *h_elements1, int *h_elements2, int *h_elements3, int elementsCount,
                             float *h_nodesX, float *h_nodesY, float *h_nodesZ, int nodesCount,
                             int *h_colors, int colorsCount, int *h_constraints, int constraintsCount, float *h_b);

void FiniteElementMethodCUDA2(int *h_x_s, int *h_y_s, int *h_rowSizes, float *h_D, int *h_elements,
                              int *h_elements0, int *h_elements1, int *h_elements2, int *h_elements3, int elementsCount,
                              float *h_nodesX, float *h_nodesY, float *h_nodesZ, int nodesCount,
                              int *h_colors, int colorsCount, int *h_constraints, int constraintsCount, float *h_b, int nnz_size, float *h_x);

void LU_GPU_SOLVE(int *h_A_RowIndices, int *h_A_ColIndices, float *h_A, int n, int nnz, float *h_x, float *result);

void callCGM_GPU_NEW(float *h_A, float *h_b, int SIZE);

void callCudaKernel();
void callCGM_GPU(int *x, int *y, float *data, float *b, float *x_k, int n, int sparse_size);



#endif
