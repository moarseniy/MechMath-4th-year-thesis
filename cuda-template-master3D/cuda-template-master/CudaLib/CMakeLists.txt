cmake_minimum_required(VERSION 3.5)
find_package(CUDA QUIET REQUIRED)
#findCUDAToolkit()

MESSAGE(STATUS "${CUDA_NVCC_FLAGS}")

set(CUDA_64_BIT_DEVICE_CODE ON)
set(CUDA_NVCC_FLAGS "${CUDA_NVCC_FLAGS} -std=c++11 -arch=sm_52" )

MESSAGE(STATUS "${CUDA_NVCC_FLAGS}")

#INCLUDE_DIRECTORIES(
#    include
#)

CUDA_ADD_LIBRARY(
    simplecuda
    include/init.h
    src/init.cu
)

TARGET_INCLUDE_DIRECTORIES(simplecuda PUBLIC include)
#TARGET_INCLUDE_DIRECTORIES(simplecuda PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/src")
TARGET_LINK_LIBRARIES(simplecuda Linal2 ${CUDA_LIBRARIES} ${CUDA_CUBLAS_LIBRARIES} ${CUDA_cusparse_LIBRARY} ${CUDA_cusolver_LIBRARY})
