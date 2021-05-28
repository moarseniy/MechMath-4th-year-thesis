#include "CudaTest.h"
#include <stdio.h>


cusparseHandle_t    handle;

cusparseMatDescr_t  descrA = 0;
cusparseMatDescr_t  descr_L = 0;
cusparseMatDescr_t  descr_U = 0;

csrilu02Info_t      info_A = 0;
csrsv2Info_t        info_L = 0;
csrsv2Info_t        info_U = 0;

void LU_GPU_SOLVE() {

}
