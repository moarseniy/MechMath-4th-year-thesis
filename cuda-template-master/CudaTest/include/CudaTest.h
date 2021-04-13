#ifndef CUDA_LIB_
#define CUDA_LIB_
void callCudaKernel();
void callCGM_GPU(int *x, int *y, double *data, double *b, double *x_k, int n, int sparse_size);



#endif
