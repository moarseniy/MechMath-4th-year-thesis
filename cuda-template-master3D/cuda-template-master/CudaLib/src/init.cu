
#include <stdio.h>

#include <cuda_runtime.h>
#include <cusparse_v2.h>



#include <iostream>
#include <cuda.h>
#include <cusolverSp.h>

#include <math.h>
#include <vector>
#include "init.h"
//#include "femfunc.h"

#include "Linal2.h"

using namespace std;

__device__
float det3(float a0, float a1, float a2, float a3,
           float a4, float a5, float a6, float a7, float a8) {
    return a0*a4*a8 +
                    a1*a6*a5 +
                    a2*a3*a7 -
                    a6*a4*a2 -
                    a0*a5*a7 -
                    a1*a3*a8;
}

__device__
float det3x3(float *c) {
    return c[0]*c[4]*c[8] +
            c[1]*c[6]*c[5] +
            c[2]*c[3]*c[7] -
            c[6]*c[4]*c[2] -
            c[0]*c[5]*c[7] -
            c[1]*c[3]*c[8];
}

__device__
float det4x4(float *c) {
    float v1 = det3(c[5], c[6], c[7], c[9], c[10], c[11], c[13], c[14], c[15]);
    float v2 = det3(c[1], c[2], c[3], c[9], c[10], c[11], c[13], c[14], c[15]);
    float v3 = det3(c[1], c[2], c[3], c[5], c[6], c[7], c[13], c[14], c[15]);
    float v4 = det3(c[1], c[2], c[3], c[5], c[6], c[7], c[9], c[10], c[11]);
    return v1 - v2 + v3 - v4;
}

__device__
float det(float *c, int size) {
    if (size == 1) {
        return c[0];
    } else if (size == 2) {
        return c[0 + 0 * 2] * c[1 + 1 * 2] - c[0 + 1 * 2] * c[1 + 0 * 2];
    } else if (size == 3) {
        return c[0]*c[4]*c[8] +
                c[1]*c[6]*c[5] +
                c[2]*c[3]*c[7] -
                c[6]*c[4]*c[2] -
                c[0]*c[5]*c[7] -
                c[1]*c[3]*c[8];
    } else if (size == 4) {
        float v1 = det3(c[5], c[6], c[7], c[9], c[10], c[11], c[13], c[14], c[15]);
        float v2 = det3(c[1], c[2], c[3], c[9], c[10], c[11], c[13], c[14], c[15]);
        float v3 = det3(c[1], c[2], c[3], c[5], c[6], c[7], c[13], c[14], c[15]);
        float v4 = det3(c[1], c[2], c[3], c[5], c[6], c[7], c[9], c[10], c[11]);
        return v1 - v2 + v3 - v4;
    }
}

__device__
void Get_matrix(float *a, int n, float *c, int indRow, int indCol) {
    //float *a = (float*)malloc(3 * 3 * sizeof (float));
    int ki = 0;
    for (int i = 0; i < n; i++) {
        if (i != indRow) {
            for (int j = 0, kj = 0; j < n; j++) {
                if (j != indCol) {
                    a[kj + ki * 3] = c[j + i * n];
                    kj++;
                }
            }
            ki++;
        }
    }

    //return a;
}

__device__
void inverse(float *ic, float *b, int size) {
    //float *ic = (float*)malloc(4 * 4 * sizeof (float));
    //printf("%f", b[0]);
    float determinant = det4x4(b);

    if (determinant) {
        for (int i = 0; i < size; i++) {
            for (int j = 0; j < size; j++) {
                float *temp =(float*)malloc(3 * 3 * sizeof (float));
                //__shared__ float temp[3 * 3];
                Get_matrix(temp, size, b, i, j);
                ic[j + i * size] = ((i + j + 2) % 2 == 0 ? 1.0 : -1.0) * det3x3(temp) / determinant;
                free(temp);
            }
        }
    }

    float swap;
    for (int i = 0; i < size; i++) {
        for (int j = 0; j < size; j++) {
            if (i > j) {
               swap = ic[j + i * size];
               ic[j + i * size] = ic[i + j * size];
               ic[i + j * size] = swap;
            }
        }
    }

   //return ic;
}

__global__
void CalculateLocalSets(int elementsCount, const float *nodesX, const float *nodesY, const float *nodesZ, int nodesCount,
                        int *elements, int k, int sumColors, float *D, int *K_x, int *K_y, float *K_value, int *constraints, int constraintCount) {
    int index = blockIdx.x * blockDim.x + threadIdx.x;
    //printf("!%d ", index);
    //printf("%d ", elementsCount);

    if (index < elementsCount) {
        //__shared__ float C[4 * 4];
        //__shared__ float IC[4 * 4];

        float *C = (float*)malloc(4 * 4 * sizeof(float));
        float *IC = (float*)malloc(4 * 4 * sizeof(float));

        C[0 + 0 * 4] = C[0 + 1 * 4] = C[0 + 2 * 4] = C[0 + 3 * 4] = 1.0;
        C[1 + 0 * 4] = nodesX[elements[4 * (index + sumColors) + 0]]; C[1 + 1 * 4] = nodesX[elements[4 * (index + sumColors) + 1]]; C[1 + 2 * 4] = nodesX[elements[4 * (index + sumColors) + 2]]; C[1 + 3 * 4] = nodesX[elements[4 * (index + sumColors) + 3]];
        C[2 + 0 * 4] = nodesY[elements[4 * (index + sumColors) + 0]]; C[2 + 1 * 4] = nodesY[elements[4 * (index + sumColors) + 1]]; C[2 + 2 * 4] = nodesY[elements[4 * (index + sumColors) + 2]]; C[2 + 3 * 4] = nodesY[elements[4 * (index + sumColors) + 3]];
        C[3 + 0 * 4] = nodesZ[elements[4 * (index + sumColors) + 0]]; C[3 + 1 * 4] = nodesZ[elements[4 * (index + sumColors) + 1]]; C[3 + 2 * 4] = nodesZ[elements[4 * (index + sumColors) + 2]]; C[3 + 3 * 4] = nodesZ[elements[4 * (index + sumColors) + 3]];

        float determinant = det4x4(C);

//        printf("%d %d %d %d ", elements[4 * (index + sumColors) + 0], elements[4 * (index + sumColors) + 1],
//                elements[4 * (index + sumColors) + 2], elements[4 * (index + sumColors) + 3]);
//        printf("%d %d %d %d ", 4 * (index + sumColors) + 0, 4 * (index + sumColors) + 1,
//                4 * (index + sumColors) + 2, 4 * (index + sumColors) + 3);



        //printf("%f ", nodesX[elements[4 * (index + sumColors) + 3]]);

        //printf("%f ", C[1]);

        //printf("%d ", index);



        //__syncthreads();

        inverse(IC, C, 4);
        free(C);

        //inverse(C, IC, 4);
        //printf("%f ", IC[1]);

        //__syncthreads();

        //__shared__ float B[6 * 12];
        float *B = (float*)malloc(6 * 12 * sizeof(float));

        for (int i = 0; i < 4; i++) {
            B[(3 * i + 0) + 0 * 12] = IC[i + 1 * 4];
            B[(3 * i + 1) + 0 * 12] = 0.0;
            B[(3 * i + 2) + 0 * 12] = 0.0;

            B[(3 * i + 0) + 1 * 12] = 0.0;
            B[(3 * i + 1) + 1 * 12] = IC[i + 2 * 4];
            B[(3 * i + 2) + 1 * 12] = 0.0;

            B[(3 * i + 0) + 2 * 12] = 0.0;
            B[(3 * i + 1) + 2 * 12] = 0.0;
            B[(3 * i + 2) + 2 * 12] = IC[i + 3 * 4];

            B[(3 * i + 0) + 3 * 12] = IC[i + 2 * 4];
            B[(3 * i + 1) + 3 * 12] = IC[i + 1 * 4];
            B[(3 * i + 2) + 3 * 12] = 0.0;

            B[(3 * i + 0) + 4 * 12] = 0.0;
            B[(3 * i + 1) + 4 * 12] = IC[i + 3 * 4];
            B[(3 * i + 2) + 4 * 12] = IC[i + 2 * 4];

            B[(3 * i + 0) + 5 * 12] = IC[i + 3 * 4];
            B[(3 * i + 1) + 5 * 12] = 0.0;
            B[(3 * i + 2) + 5 * 12] = IC[i + 1 * 4];
        }

        free(IC);

        //__syncthreads();


        //printf("%d-%f ",i + sumColors, determinant);

//        for (int i = 0; i < 6; i++) {
//            for (int j = 0; j < 12; j++) {
//                printf("%f ", B[j + i * 12]);
//            }
//            //printf("\n");
//        }
        //printf("%f ", B[0]);


        //transpose
        //__shared__ float B_T[12 * 6];
        float *B_T = (float*)malloc(12 * 6 * sizeof(float));
        for (int i = 0; i < 12; i++) {
            for (int j = 0; j < 6; j++) {
                B_T[j + i * 6] = B[i + j * 12];
                //printf("%f ", B_T[j + i * 6]);
            }
            //printf("\n");
        }
        //printf("%f-%f ", B[2], B_T[2]);

        //__syncthreads();

        //product B_T * D
        float *temp = (float*)malloc(12 * 6 * sizeof(float));
        //__shared__ float temp[12 * 6];

        for (int i = 0; i < 12; i++) {
            for (int j = 0; j < 6; j++) {
                temp[j + i * 6] = 0.0;
                for (int k = 0; k < 6; k++) {
                    temp[j + i * 6] += B_T[k + i * 6] * D[j + k * 6];
                }
            }
        }

        free(B_T);

        //__syncthreads();

        //product (B_T * D) * B
        float *K = (float*)malloc(12 * 12 * sizeof(float));
        //__shared__ float K[12 * 12];

        for (int i = 0; i < 12; i++) {
            for (int j = 0; j < 12; j++) {
                K[j + i * 12] = 0.0;
                for (int k = 0; k < 6; k++) {
                    K[j + i * 12] += temp[k + i * 6] * B[j + k * 12];
                }
            }
        }

        free(temp);
        free(B);

        //__syncthreads();

        //scale K * |det(C)| / 6.0
        for (int i = 0; i < 12; i++) {
            for (int j = 0; j < 12; j++) {
                //K[j + i * 12] *= (determinant > 0 ? determinant : (-1 * determinant)) / 6.0;
                K[j + i * 12] *= (fabs(determinant)) / 6.0;
                //printf("%f ", K[j + i * 12]);
            }
        }



//        for (int i = 0; i < 12; i ++) {
//            for (int j = 0; j < 12; j++) {
//                printf("%f ", K[j + i * 12]);
//            }
//            printf("\n");
//        }

        //printf("%f ", K_value[0]);

        //__syncthreads();
        for (int i = 0; i < 4; i++) {
            for (int j = 0; j < 4; j++) {
                int idi = 3 * elements[4 * (index + sumColors) + i];
                int idj = 3 * elements[4 * (index + sumColors) + j];
                //printf("%d-%d ", idi, idj);
                K_x[(idj + 0) + (idi + 0) * 3 * nodesCount] = idi + 0;
                K_y[(idj + 0) + (idi + 0) * 3 * nodesCount] = idj + 0;
                K_value[(idj + 0) + (idi + 0) * 3 * nodesCount] += K[(3 * j + 0) + (3 * i + 0) * 12];

                K_x[(idj + 1) + (idi + 0) * 3 * nodesCount] = idi + 0;
                K_y[(idj + 1) + (idi + 0) * 3 * nodesCount] = idj + 1;
                K_value[(idj + 1) + (idi + 0) * 3 * nodesCount] += K[(3 * j + 1) + (3 * i + 0) * 12];

                K_x[(idj + 2) + (idi + 0) * 3 * nodesCount] = idi + 0;
                K_y[(idj + 2) + (idi + 0) * 3 * nodesCount] = idj + 2;
                K_value[(idj + 2) + (idi + 0) * 3 * nodesCount] += K[(3 * j + 2) + (3 * i + 0) * 12];

                K_x[(idj + 0) + (idi + 1) * 3 * nodesCount] = idi + 1;
                K_y[(idj + 0) + (idi + 1) * 3 * nodesCount] = idj + 0;
                K_value[(idj + 0) + (idi + 1) * 3 * nodesCount] += K[(3 * j + 0) + (3 * i + 1) * 12];

                K_x[(idj + 1) + (idi + 1) * 3 * nodesCount] = idi + 1;
                K_y[(idj + 1) + (idi + 1) * 3 * nodesCount] = idj + 1;
                K_value[(idj + 1) + (idi + 1) * 3 * nodesCount] += K[(3 * j + 1) + (3 * i + 1) * 12];

                K_x[(idj + 2) + (idi + 1) * 3 * nodesCount] = idi + 1;
                K_y[(idj + 2) + (idi + 1) * 3 * nodesCount] = idj + 2;
                K_value[(idj + 2) + (idi + 1) * 3 * nodesCount] += K[(3 * j + 2) + (3 * i + 1) * 12];

                K_x[(idj + 0) + (idi + 2) * 3 * nodesCount] = idi + 2;
                K_y[(idj + 0) + (idi + 2) * 3 * nodesCount] = idj + 0;
                K_value[(idj + 0) + (idi + 2) * 3 * nodesCount] += K[(3 * j + 0) + (3 * i + 2) * 12];

                K_x[(idj + 1) + (idi + 2) * 3 * nodesCount] = idi + 2;
                K_y[(idj + 1) + (idi + 2) * 3 * nodesCount] = idj + 1;
                K_value[(idj + 1) + (idi + 2) * 3 * nodesCount] += K[(3 * j + 1) + (3 * i + 2) * 12];

                K_x[(idj + 2) + (idi + 2) * 3 * nodesCount] = idi + 2;
                K_y[(idj + 2) + (idi + 2) * 3 * nodesCount] = idj + 2;
                K_value[(idj + 2) + (idi + 2) * 3 * nodesCount] += K[(3 * j + 2) + (3 * i + 2) * 12];

                for (int t = 0; t < constraintCount; t++) {
                    for (int i1 = 0; i1 < 3; i1++) {
                        for (int j1 = 0; j1 < 3; j1++) {
                            if (idi + i1 == constraints[t] || idj + j1 == constraints[t]) {
                                if (idi + i1 == idj + j1) {
                                    K_value[(idj + j1) + (idi + i1) * 3 * nodesCount] = 1.0;
                                } else {
                                    K_value[(idj + j1) + (idi + i1) * 3 * nodesCount] = 0.0;
                                }
                            }
                        }
                    }
                }
//                K_x[(idj + 0) + (idi + 0) * 3 * nodesCount] = idi + 0;
//                K_y[(idj + 0) + (idi + 0) * 3 * nodesCount] = idj + 0;
//                K_value[(idj + 0) + (idi + 0) * 3 * nodesCount] += K[(3 * j + 0) + (3 * i + 0) * 12];

//                K_x[(idj + 1) + (idi + 0) * 3 * nodesCount] = idi + 0;
//                K_y[(idj + 1) + (idi + 0) * 3 * nodesCount] = idj + 1;
//                K_value[(idj + 1) + (idi + 0) * 3 * nodesCount] += K[(3 * j + 1) + (3 * i + 0) * 12];

//                K_x[(idj + 2) + (idi + 0) * 3 * nodesCount] = idi + 0;
//                K_y[(idj + 2) + (idi + 0) * 3 * nodesCount] = idj + 2;
//                K_value[(idj + 2) + (idi + 0) * 3 * nodesCount] += K[(3 * j + 2) + (3 * i + 0) * 12];

//                K_x[(idj + 0) + (idi + 1) * 3 * nodesCount] = idi + 1;
//                K_y[(idj + 0) + (idi + 1) * 3 * nodesCount] = idj + 0;
//                K_value[(idj + 0) + (idi + 1) * 3 * nodesCount] += K[(3 * j + 0) + (3 * i + 1) * 12];

//                K_x[(idj + 1) + (idi + 1) * 3 * nodesCount] = idi + 1;
//                K_y[(idj + 1) + (idi + 1) * 3 * nodesCount] = idj + 1;
//                K_value[(idj + 1) + (idi + 1) * 3 * nodesCount] += K[(3 * j + 1) + (3 * i + 1) * 12];

//                K_x[(idj + 2) + (idi + 1) * 3 * nodesCount] = idi + 1;
//                K_y[(idj + 2) + (idi + 1) * 3 * nodesCount] = idj + 2;
//                K_value[(idj + 2) + (idi + 1) * 3 * nodesCount] += K[(3 * j + 2) + (3 * i + 1) * 12];

//                K_x[(idj + 0) + (idi + 2) * 3 * nodesCount] = idi + 2;
//                K_y[(idj + 0) + (idi + 2) * 3 * nodesCount] = idj + 0;
//                K_value[(idj + 0) + (idi + 2) * 3 * nodesCount] += K[(3 * j + 0) + (3 * i + 2) * 12];

//                K_x[(idj + 1) + (idi + 2) * 3 * nodesCount] = idi + 2;
//                K_y[(idj + 1) + (idi + 2) * 3 * nodesCount] = idj + 1;
//                K_value[(idj + 1) + (idi + 2) * 3 * nodesCount] += K[(3 * j + 1) + (3 * i + 2) * 12];

//                K_x[(idj + 2) + (idi + 2) * 3 * nodesCount] = idi + 2;
//                K_y[(idj + 2) + (idi + 2) * 3 * nodesCount] = idj + 2;
//                K_value[(idj + 2) + (idi + 2) * 3 * nodesCount] += K[(3 * j + 2) + (3 * i + 2) * 12];

            }
            //__syncthreads();
        }

        free(K);
        //__syncthreads();

    }
    //__syncthreads();

}

__global__
void ApplyConstraintsCuda(int *K_x, int *K_y, float *K_value, int *constraints, int constraintsCount, int i) {
    int index = blockIdx.x * blockDim.x + threadIdx.x;
    if (index < constraintsCount) {
        //printf("[%d]%d ",index, constraints[index]);
        if (K_x[i] == constraints[index] || K_y[i] == constraints[index]) {
            if (K_x[i] == K_y[i]) {
                K_value[i] = 1.0;
            } else {
                K_value[i] = 0.0;
            }
        }
    }
    __syncthreads();
}

__global__
void CountNonZeroValues(int *K_x, int *K_y, float *K_value, int SIZE, int *nnz) {
    nnz[0] = 0;
    float epsilon = 1e-50;
    for (int i = 0; i < SIZE; i++) {
        if (K_value[i] != 0.0) {
            //printf("%f ", K_value[i]);
            nnz[0]++;
        }
    }
    //printf("(func)%d %d\n", nnz[0], elementsCount);
    //__syncthreads();
}

__global__
void ConvertToCSR(int *K_x, int *K_y, float *K_value, int *ptr, int *ind, float *data, int SIZE, int nodesCount) {
    int k = 0;
    float epsilon = 1e-50;

    for (int i = 0; i < 3 * nodesCount + 1; i++) {
        ptr[i] = 0;
    }

    for (int i = 0; i < SIZE; i++) {
        if (K_value[i] != 0.0) {
            data[k] = K_value[i];
            ind[k] = K_y[i];
            ptr[K_x[i] + 1]++;
            k++;
        }
    }
    printf("k=%d\n", k);

    for (int i = 0; i < 3 * nodesCount; i++) {
        ptr[i + 1] += ptr[i];
    }
    printf("ConvertSuccess\n");
    __syncthreads();
}

__global__
void TestCudaFunc(int a) {
    int index = blockIdx.x * blockDim.x + threadIdx.x;
    if (index < a) {
        printf("%d ", index);
    }
}

void FiniteElementMethodCUDA(float *h_D, int *h_elements,
                             int *h_elements0, int *h_elements1, int *h_elements2, int *h_elements3, int elementsCount,
                             float *h_nodesX, float *h_nodesY, float *h_nodesZ, int nodesCount,
                             int *h_colors, int colorsCount, int *h_constraints, int constraintsCount, float *h_b) {
    cudaEvent_t start, stop;
    cudaEventCreate(&start);
    cudaEventCreate(&stop);
    float gpuTime = 0.0;
    cudaEventRecord(start, 0);

    int *d_elements, *d_elements0, *d_elements1, *d_elements2, *d_elements3, *d_colors, *d_K_x, *d_K_y, *d_constraints;
    float *d_D, *d_nodesX, *d_nodesY, *d_nodesZ,  *d_K_value, *d_b, *d_x;

    int SIZE = 3 * nodesCount * 3 * nodesCount;

    int *d_nnz;
    int *h_nnz = new int[1];
    cudaMalloc((void**)&d_nnz, 1 * sizeof(int));

    int *h_K_x = new int[SIZE];
    int *h_K_y = new int[SIZE];
    float *h_K_value = new float[SIZE];

    cudaMalloc((void**)&d_b, 3 * nodesCount * sizeof(float));
    cudaMemcpy(d_b, h_b, 3 * nodesCount * sizeof(float), cudaMemcpyHostToDevice);

    cudaMalloc((void**)&d_x, 3 * nodesCount * sizeof(float));

    //cout << "elementsCount = " << elementsCount << endl;
    cudaMalloc((void**)&d_elements, 4 * elementsCount * sizeof(int));
    cudaMemcpy(d_elements, h_elements, 4 * elementsCount * sizeof(int), cudaMemcpyHostToDevice);

    cudaMalloc((void**)&d_D, 6 * 6 * sizeof(float));
    cudaMemcpy(d_D, h_D, 6 * 6 * sizeof(float), cudaMemcpyHostToDevice);

    cudaMalloc((void**)&d_K_x, SIZE * sizeof(int));
    cudaMalloc((void**)&d_K_y, SIZE * sizeof(int));
    cudaMalloc((void**)&d_K_value, SIZE * sizeof(float));

//    cudaMalloc((void**)&d_elements0, elementsCount * sizeof(int));
//    cudaMalloc((void**)&d_elements1, elementsCount * sizeof(int));
//    cudaMalloc((void**)&d_elements2, elementsCount * sizeof(int));
//    cudaMalloc((void**)&d_elements3, elementsCount * sizeof(int));
//    cudaMemcpy(d_elements0, h_elements0, elementsCount * sizeof(int), cudaMemcpyHostToDevice);
//    cudaMemcpy(d_elements1, h_elements1, elementsCount * sizeof(int), cudaMemcpyHostToDevice);
//    cudaMemcpy(d_elements2, h_elements2, elementsCount * sizeof(int), cudaMemcpyHostToDevice);
//    cudaMemcpy(d_elements3, h_elements3, elementsCount * sizeof(int), cudaMemcpyHostToDevice);


    cudaMalloc((void**)&d_nodesX, nodesCount * sizeof(float));
    cudaMalloc((void**)&d_nodesY, nodesCount * sizeof(float));
    cudaMalloc((void**)&d_nodesZ, nodesCount * sizeof(float));
    cudaMemcpy(d_nodesX, h_nodesX, nodesCount * sizeof(float), cudaMemcpyHostToDevice);
    cudaMemcpy(d_nodesY, h_nodesY, nodesCount * sizeof(float), cudaMemcpyHostToDevice);
    cudaMemcpy(d_nodesZ, h_nodesZ, nodesCount * sizeof(float), cudaMemcpyHostToDevice);

    cudaMalloc((void**)&d_colors, colorsCount * sizeof(int));
    cudaMemcpy(d_colors, h_colors, colorsCount * sizeof(int), cudaMemcpyHostToDevice);

    cudaMalloc((void**)&d_constraints, constraintsCount * sizeof(int));
    cudaMemcpy(d_constraints, h_constraints, constraintsCount * sizeof(int), cudaMemcpyHostToDevice);

    cudaEventRecord(stop, 0);
    cudaEventSynchronize(stop);
    cudaEventElapsedTime(&gpuTime, start, stop);
    printf("GPU(Read data) = %.4f ms\n", gpuTime);

    ////////
//    int size = 1000, temp_sum = 5;
//    int *test1 = new int[size];

//    for (int k = 0; k < 10; k++) {
        //TestCudaFunc<<<1, 50>>>(temp_sum);
//        temp_sum += 2;

        //cudaDeviceSynchronize();
        //cout << endl;
//    }

    ///////

    int sumColors = 0;
    cout << "\ncolorsCount = " << colorsCount << endl;
    for (int k = 0; k < colorsCount; k++) {
        CalculateLocalSets<<<(255+h_colors[k])/256, 256>>> (h_colors[k], d_nodesX, d_nodesY, d_nodesZ, nodesCount,
                                                  d_elements, k, sumColors, d_D, d_K_x, d_K_y, d_K_value, d_constraints, constraintsCount);
        sumColors += h_colors[k];
        //cout << endl << "!!!" << k << " " << h_colors[k] << endl;
        //break;
        cudaDeviceSynchronize();
    }


    cudaEventRecord(stop, 0);
    cudaEventSynchronize(stop);
    cudaEventElapsedTime(&gpuTime, start, stop);
    printf("GPU(ColorStiffnessMatrix) = %.4f ms\n", gpuTime);

//    cudaMemcpy(h_K_x, d_K_x, SIZE * sizeof(int), cudaMemcpyDeviceToHost);
//    cudaMemcpy(h_K_y, d_K_y, SIZE * sizeof(int), cudaMemcpyDeviceToHost);
//    cudaMemcpy(h_K_value, d_K_value, SIZE * sizeof(float), cudaMemcpyDeviceToHost);

//    for (int i = 0; i < SIZE; i++) {
//        //if (h_K_value[i] != 0.0) {
//            cout << h_K_x[i] << " " << h_K_y[i] << " " << h_K_value[i] << endl;
//        //}
//    }

    cudaFree(d_nodesX);
    cudaFree(d_nodesY);
    cudaFree(d_nodesZ);


    CountNonZeroValues<<<1, 1>>> (d_K_x, d_K_y, d_K_value, SIZE, d_nnz);
    cudaMemcpy(h_nnz, d_nnz, 1 * sizeof(int), cudaMemcpyDeviceToHost);
    cout << "\nNONZERO=" << h_nnz[0] << endl;


//    for (int i = 0; i < SIZE; i++) {
//        ApplyConstraintsCuda<<<(255+constraintsCount)/256, 256>>> (d_K_x, d_K_y, d_K_value, d_constraints, constraintsCount, i);
//    }

//    cudaEventRecord(stop, 0);
//    cudaEventSynchronize(stop);
//    cudaEventElapsedTime(&gpuTime, start, stop);
//    printf("GPU(ApplyConstraints) = %.4f ms\n", gpuTime);


//    CountNonZeroValues<<<1, 1>>> (d_K_x, d_K_y, d_K_value, SIZE, d_nnz);
//    cudaMemcpy(h_nnz, d_nnz, 1 * sizeof(int), cudaMemcpyDeviceToHost);
//    cout << "NONZERO(After Constraints)=" << h_nnz[0] << endl;


    int *d_ptr, *d_ind, *d_row;
    float *d_data;

    int *h_ptr = new int[3 * nodesCount + 1];
    float *h_data = new float[h_nnz[0]];
    int *h_ind = new int[h_nnz[0]];

    cudaMalloc((void**)&d_ptr, (3 * nodesCount + 1) * sizeof(int));
    //cudaMalloc((void**)&d_row, h_nnz[0] * sizeof(int));
    cudaMalloc((void**)&d_ind, h_nnz[0] * sizeof(int));
    cudaMalloc((void**)&d_data, h_nnz[0] * sizeof(float));

    ConvertToCSR<<<1, 1>>> (d_K_x, d_K_y, d_K_value, d_ptr, d_ind, d_data, SIZE, nodesCount);

    cudaEventRecord(stop, 0);
    cudaEventSynchronize(stop);
    cudaEventElapsedTime(&gpuTime, start, stop);
    printf("GPU(ConvertToCSR) = %.4f ms\n", gpuTime);






    cudaFree(d_D);
    cudaFree(d_K_value);
    cudaFree(d_K_x);
    cudaFree(d_K_y);
    cudaFree(d_colors);
    cudaFree(d_constraints);
    cudaFree(d_elements);

//    cudaMemcpy(h_ptr, d_ptr, (3 * nodesCount + 1) * sizeof(int), cudaMemcpyDeviceToHost);
//    for (int i = 0; i < 3 * nodesCount + 1; i++) {
//        cout << h_ptr[i] << " ";
//    }
//    cudaMemcpy(h_ind, d_ind, h_nnz[0] * sizeof(int), cudaMemcpyDeviceToHost);
//    for (int i = 0; i < h_nnz[0]; i++) {
//        cout << h_ind[i] << " ";
//    }
//    cudaMemcpy(h_data, d_data, h_nnz[0] * sizeof(float), cudaMemcpyDeviceToHost);
//    for (int i = 0; i < h_nnz[0]; i++) {
//        cout << h_data[i] << " ";
//    }




    cusolverSpHandle_t handle;
    cusolverSpCreate(&handle);
    cusparseMatDescr_t descr;
    cusparseCreateMatDescr(&descr);

//    cudaMemcpy(d_csrValA, h_csrValA, nnz * sizeof(float), cudaMemcpyHostToDevice);
//    cudaMemcpy(d_csrRowPtrA, h_csrRowPtrA, (n + 1) * sizeof(int), cudaMemcpyHostToDevice);
//    cudaMemcpy(d_csrColIndA, h_csrColIndA, nnz * sizeof(int), cudaMemcpyHostToDevice);
//    cudaMemcpy(d_b, h_b, n * sizeof(float), cudaMemcpyHostToDevice);



    cout<<"start solving...\n";
    float tol = 1e-16;
    int reorder = 1;
    int singularity = 0;
    cusolverSpScsrlsvqr(handle, 3 * nodesCount, h_nnz[0], descr, d_data, d_ptr, d_ind, d_b, tol,
                     reorder, d_x, &singularity);

    cout << "end solving...\n";

    cudaEventRecord(stop, 0);
    cudaEventSynchronize(stop);
    cudaEventElapsedTime(&gpuTime, start, stop);
    printf("GPU(SOLVER) = %.4f ms\n", gpuTime);

    float *h_x = new float[3 * nodesCount];
    cudaMemcpy(h_x, d_x, 3 * nodesCount * sizeof(float), cudaMemcpyDeviceToHost);

    for (int i = 0; i < 3 * nodesCount; i++) {
        cout << h_x[i] << " ";
    }


    cudaFree(d_b);
    cudaFree(d_x);
    cudaFree(d_data);
    cudaFree(d_ptr);
    cudaFree(d_ind);
    cudaFree(d_nnz);
    delete [] h_nnz;

}



void CudaSolve(int *h_csrRowPtrA, int *h_csrColIndA, float *h_csrValA, int n, int nnz, float *h_b, float *h_x) {

    cusolverSpHandle_t handle;
    cusolverStatus_t status;
    cusparseStatus_t status2;

    status = cusolverSpCreate(&handle);


    cusparseMatDescr_t descr;
    status2 = cusparseCreateMatDescr(&descr);



    float* d_csrValA, *d_b, *d_x;
    int* d_csrRowPtrA, *d_csrColIndA;
    cudaMalloc((void**)&d_csrValA, nnz * sizeof(float));
    cudaMalloc((void**)&d_b, n * sizeof(float));
    cudaMalloc((void**)&d_x, n * sizeof(float));
    cudaMalloc((void**)&d_csrRowPtrA, (n + 1) * sizeof(int));
    cudaMalloc((void**)&d_csrColIndA, nnz * sizeof(int));


    cudaMemcpy(d_csrValA, h_csrValA, nnz * sizeof(float), cudaMemcpyHostToDevice);
    cudaMemcpy(d_csrRowPtrA, h_csrRowPtrA, (n + 1) * sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(d_csrColIndA, h_csrColIndA, nnz * sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(d_b, h_b, n * sizeof(float), cudaMemcpyHostToDevice);



    cout<<"start solving...\n";
    float tol = 1e-16;
    int reorder = 1;
    int singularity = 0;
    status = cusolverSpScsrlsvqr(handle, n, nnz, descr, d_csrValA, d_csrRowPtrA, d_csrColIndA, d_b, tol,
                     reorder, d_x, &singularity);

    cout<<"end solving...\n";
    cudaMemcpy(h_x, d_x, n * sizeof(float), cudaMemcpyDeviceToHost);
    //cout<<"singularity = "<<singularity<<"\n";


    cudaFree(d_csrValA);
    cudaFree(d_csrRowPtrA);
    cudaFree(d_csrColIndA);
    cudaFree(d_b);
    cudaFree(d_x);
    cusolverSpDestroy(handle);


}

void SortCOO(int *h_cooRows, int *h_cooCols, float *h_cooVals, int n, int nnz) {
    cusparseHandle_t handle = NULL;
    cudaStream_t stream = NULL;


    int *h_P = new int[nnz];

    int *d_cooRows = NULL;
    int *d_cooCols = NULL;
    int *d_P       = NULL;
    float *d_cooVals = NULL;
    float *d_cooVals_sorted = NULL;
    size_t pBufferSizeInBytes = 0;
    void *pBuffer = NULL;

    cudaStreamCreateWithFlags(&stream, cudaStreamNonBlocking);
    cusparseCreate(&handle);
    cusparseSetStream(handle, stream);

    cusparseXcoosort_bufferSizeExt(
            handle,
            n,
            n,
            nnz,
            d_cooRows,
            d_cooCols,
            &pBufferSizeInBytes
        );

     printf("pBufferSizeInBytes = %lld bytes \n", (long long)pBufferSizeInBytes);

     cudaMalloc( &d_cooRows, sizeof(int)*nnz);
     cudaMalloc( &d_cooCols, sizeof(int)*nnz);
        cudaMalloc( &d_P      , sizeof(int)*nnz);
        cudaMalloc( &d_cooVals, sizeof(float)*nnz);
        cudaMalloc( &d_cooVals_sorted, sizeof(float)*nnz);
        cudaMalloc( &pBuffer, sizeof(char)* pBufferSizeInBytes);

        cudaMemcpy(d_cooRows, h_cooRows, sizeof(int)*nnz   , cudaMemcpyHostToDevice);
        cudaMemcpy(d_cooCols, h_cooCols, sizeof(int)*nnz   , cudaMemcpyHostToDevice);
        cudaMemcpy(d_cooVals, h_cooVals, sizeof(float)*nnz, cudaMemcpyHostToDevice);
        cudaDeviceSynchronize();


        cusparseCreateIdentityPermutation(handle, nnz, d_P);


        cusparseXcoosortByRow(
            handle,
            n,
            n,
            nnz,
            d_cooRows,
            d_cooCols,
            d_P,
            pBuffer
        );

        cusparseSgthr(
            handle,
            nnz,
            d_cooVals,
            d_cooVals_sorted,
            d_P,
            CUSPARSE_INDEX_BASE_ZERO
        );

        cudaDeviceSynchronize();
        cudaMemcpy(h_cooRows, d_cooRows, sizeof(int)*nnz   , cudaMemcpyDeviceToHost);
        cudaMemcpy(h_cooCols, d_cooCols, sizeof(int)*nnz   , cudaMemcpyDeviceToHost);
        cudaMemcpy(h_P,       d_P      , sizeof(int)*nnz   , cudaMemcpyDeviceToHost);
        cudaMemcpy(h_cooVals, d_cooVals_sorted, sizeof(float)*nnz, cudaMemcpyDeviceToHost);
        cudaDeviceSynchronize();


//        printf("sorted coo: \n");
//        for(int j = 0 ; j < nnz; j++){
//            printf("(%d, %d, %f) \n", h_cooRows[j], h_cooCols[j], h_cooVals[j] );
//        }

//        for(int j = 0 ; j < nnz; j++){
//            printf("P[%d] = %d \n", j, h_P[j] );
//        }

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
    cusparseMatDescr_t descr_U, float* d_A, int* d_A_RowIndices, int* d_A_ColIndices, cusparseOperation_t matrixOperation, void** pBuffer) {

    cusparseCreateCsrilu02Info(&info_A);
    cusparseCreateCsrsv2Info(&info_L);
    cusparseCreateCsrsv2Info(&info_U);

    int pBufferSize_M, pBufferSize_L, pBufferSize_U;
//    cusparseDcsrilu02_bufferSize(handle, N, nnz, descrA, d_A, d_A_RowIndices, d_A_ColIndices, info_A, &pBufferSize_M);
//    cusparseDcsrsv2_bufferSize(handle, matrixOperation, N, nnz, descr_L, d_A, d_A_RowIndices, d_A_ColIndices, info_L, &pBufferSize_L);
//    cusparseDcsrsv2_bufferSize(handle, matrixOperation, N, nnz, descr_U, d_A, d_A_RowIndices, d_A_ColIndices, info_U, &pBufferSize_U);

    int pBufferSize = max(pBufferSize_M, max(pBufferSize_L, pBufferSize_U));
    cudaMalloc((void**)pBuffer, pBufferSize);

}


// ANALYSIS FUNCTION FOR LU DECOMPOSITION
void analysisLUDecomposition(csrilu02Info_t& info_A, csrsv2Info_t& info_L, csrsv2Info_t& info_U, cusparseHandle_t handle, const int N, const int nnz, cusparseMatDescr_t descrA, cusparseMatDescr_t descr_L,
    cusparseMatDescr_t descr_U, float* d_A, int* d_A_RowIndices, int* d_A_ColIndices, cusparseOperation_t matrixOperation, cusparseSolvePolicy_t solvePolicy1, cusparseSolvePolicy_t solvePolicy2, void* pBuffer) {

    int structural_zero;

//    cusparseDcsrilu02_analysis(handle, N, nnz, descrA, d_A, d_A_RowIndices, d_A_ColIndices, info_A, solvePolicy1, pBuffer);
    cusparseStatus_t status = cusparseXcsrilu02_zeroPivot(handle, info_A, &structural_zero);
    if (CUSPARSE_STATUS_ZERO_PIVOT == status) { printf("A(%d,%d) is missing\n", structural_zero, structural_zero); }

//    cusparseDcsrsv2_analysis(handle, matrixOperation, N, nnz, descr_L, d_A, d_A_RowIndices, d_A_ColIndices, info_L, solvePolicy1, pBuffer);
//    cusparseDcsrsv2_analysis(handle, matrixOperation, N, nnz, descr_U, d_A, d_A_RowIndices, d_A_ColIndices, info_U, solvePolicy2, pBuffer);

}


// COMPUTE LU DECOMPOSITION FOR SPARSE MATRICES

void computeSparseLU(csrilu02Info_t& info_A, cusparseHandle_t handle, const int N, const int nnz, cusparseMatDescr_t descrA, float* d_A, int* d_A_RowIndices,
    int* d_A_ColIndices, cusparseSolvePolicy_t solutionPolicy, void* pBuffer) {

    int numerical_zero;

//    cusparseDcsrilu02(handle, N, nnz, descrA, d_A, d_A_RowIndices, d_A_ColIndices, info_A, solutionPolicy, pBuffer);
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

void LU_GPU_SOLVE(int *h_A_RowIndices, int *h_A_ColIndices, float *h_A, int n, int nnz, float *h_x, float *result)
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


    float* d_x;
    cudaMalloc(&d_x, Nrows * sizeof(float));
    cudaMemcpy(d_x, h_x, Nrows * sizeof(float), cudaMemcpyHostToDevice);



    setUpDescriptor(descrA, CUSPARSE_MATRIX_TYPE_GENERAL, CUSPARSE_INDEX_BASE_ONE);

//    for (int i = 0; i < nnz; i++) {
//        printf("%f ", h_A[i]);
//    }


    for (int i = 0; i < nnz; i++) {
        //h_A_ColIndices[i]++;
        //printf("%d ", h_A_ColIndices[i]);
    }


    float* d_A;
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


    float* d_z;
    cudaMalloc(&d_z, N * sizeof(float));

    const float alpha = 1.;
//    cusparseDcsrsv2_solve(handle, CUSPARSE_OPERATION_NON_TRANSPOSE, N, nnz, &alpha, descr_L, d_A, d_A_RowIndices, d_A_ColIndices, info_L, d_x, d_z, CUSPARSE_SOLVE_POLICY_NO_LEVEL, pBuffer);


    float* d_y;
    cudaMalloc(&d_y, N * sizeof(float));

//    cusparseDcsrsv2_solve(handle, CUSPARSE_OPERATION_NON_TRANSPOSE, N, nnz, &alpha, descr_U, d_A, d_A_RowIndices, d_A_ColIndices, info_U, d_z, d_y, CUSPARSE_SOLVE_POLICY_USE_LEVEL, pBuffer);


    float* h_y = (float*)malloc(Ncols * sizeof(float));
    cudaMemcpy(h_x, d_y, N * sizeof(float), cudaMemcpyDeviceToHost);

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
void cgm_gpu(float *z_k, float *r_k, float *Az,
             int* x, int *y, float *data, float *b,
             float *x_k, const int n, int sparse_size, float *partialSum) {

    float mf = 0.0, alpha, beta, eps = 0.00001, Spz, Spr, Spr1;

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

void callCGM_GPU(int *x, int *y, float *data, float *b, float *x_k, int n, int sparse_size) {
    float *z_k, *r_k, *Az;
    float *d_z_k, *d_r_k, *d_Az, *d_data, *d_b, *d_x_k, *partialSum;
    int *d_x, *d_y;
    z_k = (float*)malloc(n * sizeof(float));
    r_k = (float*)malloc(n * sizeof(float));
    Az = (float*)malloc(n * sizeof(float));


    cudaMalloc(&d_z_k, n * sizeof(float));
    cudaMalloc(&d_r_k, n * sizeof(float));
    cudaMalloc(&d_Az, n * sizeof(float));
    cudaMalloc(&d_x, sparse_size * sizeof(int));
    cudaMalloc(&d_y, sparse_size * sizeof(int));
    cudaMalloc(&d_data, sparse_size * sizeof(float));
    cudaMalloc(&d_b, n * sizeof(float));
    cudaMalloc(&d_x_k, n * sizeof(float));
    cudaMalloc(&partialSum, n * sizeof(float));

    cudaMemcpy(d_z_k, z_k, n * sizeof(float), cudaMemcpyHostToDevice);
    cudaMemcpy(d_r_k, r_k, n * sizeof(float), cudaMemcpyHostToDevice);
    cudaMemcpy(d_Az, Az, n * sizeof(float), cudaMemcpyHostToDevice);
    cudaMemcpy(d_x, x, sparse_size * sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(d_y, y, sparse_size * sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(d_data, data, sparse_size * sizeof(float), cudaMemcpyHostToDevice);
    cudaMemcpy(d_b, b, n * sizeof(float), cudaMemcpyHostToDevice);
    cudaMemcpy(d_x_k, x_k, n * sizeof(float), cudaMemcpyHostToDevice);

    cgm_gpu<<<1, n>>>(d_z_k, d_r_k, d_Az, d_x, d_y, d_data, d_b, d_x_k, n, sparse_size, partialSum);

    cudaMemcpy(x_k, d_x_k, n * sizeof(float), cudaMemcpyDeviceToHost);


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

