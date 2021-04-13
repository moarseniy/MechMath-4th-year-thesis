#include "init.h"
#include <stdio.h>


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

