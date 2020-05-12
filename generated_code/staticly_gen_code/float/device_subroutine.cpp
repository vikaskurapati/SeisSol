#include "gemmgen_aux.h"
__global__ void
    __launch_bounds__(64)
kernel_sgemm_NT_NT_m35_56_n9_9_k9_pps_47bd4ea(const float ** A, int ExtraOffsetA, const float ** B, int ExtraOffsetB, float * C, int ExtraOffsetC, unsigned NumElements) {
  if ((threadIdx.y + blockDim.y * blockIdx.x) < NumElements) {
    const float* GlobMatA = &A[(threadIdx.y + blockDim.y * blockIdx.x)][0 + ExtraOffsetA];
    const float* GlobMatB = &B[(threadIdx.y + blockDim.y * blockIdx.x)][0 + ExtraOffsetB];
    float* GlobMatC = &C[(threadIdx.y + blockDim.y * blockIdx.x) * 360 + 0 + ExtraOffsetC];
    __shared__ float Scratch[81];
    float* ShrMatB = &Scratch[threadIdx.y * 81];

    // using ExtendedPatchLoader
    ShrMatB[threadIdx.x + 0] = GlobMatB[threadIdx.x + 0];
    if (threadIdx.x < 17) {
      ShrMatB[threadIdx.x + 64] = GlobMatB[threadIdx.x + 64];
    }
    __syncthreads();

    if (threadIdx.x < 35) {

      float Results[9] = {0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f};
      float Value;

      for (int k = 0; k < 9; ++k) {
        Value = GlobMatA[threadIdx.x + 56 * k];

#pragma unroll
        for (int n = 0; n < 9; ++n) {
          Results[n] += Value * ShrMatB[k + 9 * n];
        }
      }

#pragma unroll
      for (int n = 0; n < 9; ++n) {
        GlobMatC[threadIdx.x + 40 * n] = Results[n];
      }
    }
  }
}
void sgemm_NT_NT_m35_56_n9_9_k9_pps_47bd4ea(const float ** A, int ExtraOffsetA, const float ** B, int ExtraOffsetB, float * C, int ExtraOffsetC, unsigned NumElements) {
  dim3 Block(64, 1, 1);
  dim3 Grid((NumElements + 1 - 1) / 1, 1, 1);
  kernel_sgemm_NT_NT_m35_56_n9_9_k9_pps_47bd4ea<<<Grid,Block>>>(A, ExtraOffsetA, B, ExtraOffsetB, C, ExtraOffsetC, NumElements);
  CHECK_ERR;
}
__global__ void
    __launch_bounds__(64)
kernel_sgemm_NT_NT_m53_56_n9_40_k35_nsp_41d708c(const float * __restrict__ A, int ExtraOffsetA, const float * B, int ExtraOffsetB, float ** C, int ExtraOffsetC, unsigned NumElements) {
  if ((threadIdx.y + blockDim.y * blockIdx.x) < NumElements) {
    const float* GlobMatA = &A[1 + ExtraOffsetA];
    const float* GlobMatB = &B[(threadIdx.y + blockDim.y * blockIdx.x) * 360 + 0 + ExtraOffsetB];
    float* GlobMatC = &C[(threadIdx.y + blockDim.y * blockIdx.x)][1 + ExtraOffsetC];
    __shared__ float Scratch[355];
    float* ShrMatB = &Scratch[threadIdx.y * 355];

    // using ExtendedPatchLoader
#pragma unroll
    for (int i = 0; i < 5; ++i) {
      ShrMatB[threadIdx.x + i * 64] = GlobMatB[threadIdx.x + i * 64];
    }
    if (threadIdx.x < 35) {
      ShrMatB[threadIdx.x + 320] = GlobMatB[threadIdx.x + 320];
    }
    __syncthreads();

    if (threadIdx.x < 53) {

      float Results[9] = {0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f};
      float Value;

      for (int k = 0; k < 35; ++k) {
        Value = GlobMatA[threadIdx.x + 56 * k];

#pragma unroll
        for (int n = 0; n < 9; ++n) {
          Results[n] += Value * ShrMatB[k + 40 * n];
        }
      }

#pragma unroll
      for (int n = 0; n < 9; ++n) {
        GlobMatC[threadIdx.x + 56 * n] = Results[n] + GlobMatC[threadIdx.x + 56 * n];
      }
    }
  }
}
void sgemm_NT_NT_m53_56_n9_40_k35_nsp_41d708c(const float * __restrict__ A, int ExtraOffsetA, const float * B, int ExtraOffsetB, float ** C, int ExtraOffsetC, unsigned NumElements) {
  dim3 Block(64, 1, 1);
  dim3 Grid((NumElements + 1 - 1) / 1, 1, 1);
  kernel_sgemm_NT_NT_m53_56_n9_40_k35_nsp_41d708c<<<Grid,Block>>>(A, ExtraOffsetA, B, ExtraOffsetB, C, ExtraOffsetC, NumElements);
  CHECK_ERR;
}
__global__ void
    __launch_bounds__(64)
kernel_sgemm_NT_NT_m54_56_n9_40_k35_nsp_b4599dc(const float * __restrict__ A, int ExtraOffsetA, const float * B, int ExtraOffsetB, float ** C, int ExtraOffsetC, unsigned NumElements) {
  if ((threadIdx.y + blockDim.y * blockIdx.x) < NumElements) {
    const float* GlobMatA = &A[1 + ExtraOffsetA];
    const float* GlobMatB = &B[(threadIdx.y + blockDim.y * blockIdx.x) * 360 + 0 + ExtraOffsetB];
    float* GlobMatC = &C[(threadIdx.y + blockDim.y * blockIdx.x)][1 + ExtraOffsetC];
    __shared__ float Scratch[355];
    float* ShrMatB = &Scratch[threadIdx.y * 355];

    // using ExtendedPatchLoader
#pragma unroll
    for (int i = 0; i < 5; ++i) {
      ShrMatB[threadIdx.x + i * 64] = GlobMatB[threadIdx.x + i * 64];
    }
    if (threadIdx.x < 35) {
      ShrMatB[threadIdx.x + 320] = GlobMatB[threadIdx.x + 320];
    }
    __syncthreads();

    if (threadIdx.x < 54) {

      float Results[9] = {0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f};
      float Value;

      for (int k = 0; k < 35; ++k) {
        Value = GlobMatA[threadIdx.x + 56 * k];

#pragma unroll
        for (int n = 0; n < 9; ++n) {
          Results[n] += Value * ShrMatB[k + 40 * n];
        }
      }

#pragma unroll
      for (int n = 0; n < 9; ++n) {
        GlobMatC[threadIdx.x + 56 * n] = Results[n] + GlobMatC[threadIdx.x + 56 * n];
      }
    }
  }
}
void sgemm_NT_NT_m54_56_n9_40_k35_nsp_b4599dc(const float * __restrict__ A, int ExtraOffsetA, const float * B, int ExtraOffsetB, float ** C, int ExtraOffsetC, unsigned NumElements) {
  dim3 Block(64, 1, 1);
  dim3 Grid((NumElements + 1 - 1) / 1, 1, 1);
  kernel_sgemm_NT_NT_m54_56_n9_40_k35_nsp_b4599dc<<<Grid,Block>>>(A, ExtraOffsetA, B, ExtraOffsetB, C, ExtraOffsetC, NumElements);
  CHECK_ERR;
}
__global__ void
    __launch_bounds__(64)
kernel_sgemm_NT_NT_m55_56_n9_40_k35_nsp_e7d889b(const float * __restrict__ A, int ExtraOffsetA, const float * B, int ExtraOffsetB, float ** C, int ExtraOffsetC, unsigned NumElements) {
  if ((threadIdx.y + blockDim.y * blockIdx.x) < NumElements) {
    const float* GlobMatA = &A[1 + ExtraOffsetA];
    const float* GlobMatB = &B[(threadIdx.y + blockDim.y * blockIdx.x) * 360 + 0 + ExtraOffsetB];
    float* GlobMatC = &C[(threadIdx.y + blockDim.y * blockIdx.x)][1 + ExtraOffsetC];
    __shared__ float Scratch[355];
    float* ShrMatB = &Scratch[threadIdx.y * 355];

    // using ExtendedPatchLoader
#pragma unroll
    for (int i = 0; i < 5; ++i) {
      ShrMatB[threadIdx.x + i * 64] = GlobMatB[threadIdx.x + i * 64];
    }
    if (threadIdx.x < 35) {
      ShrMatB[threadIdx.x + 320] = GlobMatB[threadIdx.x + 320];
    }
    __syncthreads();

    if (threadIdx.x < 55) {

      float Results[9] = {0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f};
      float Value;

      for (int k = 0; k < 35; ++k) {
        Value = GlobMatA[threadIdx.x + 56 * k];

#pragma unroll
        for (int n = 0; n < 9; ++n) {
          Results[n] += Value * ShrMatB[k + 40 * n];
        }
      }

#pragma unroll
      for (int n = 0; n < 9; ++n) {
        GlobMatC[threadIdx.x + 56 * n] = Results[n] + GlobMatC[threadIdx.x + 56 * n];
      }
    }
  }
}
void sgemm_NT_NT_m55_56_n9_40_k35_nsp_e7d889b(const float * __restrict__ A, int ExtraOffsetA, const float * B, int ExtraOffsetB, float ** C, int ExtraOffsetC, unsigned NumElements) {
  dim3 Block(64, 1, 1);
  dim3 Grid((NumElements + 1 - 1) / 1, 1, 1);
  kernel_sgemm_NT_NT_m55_56_n9_40_k35_nsp_e7d889b<<<Grid,Block>>>(A, ExtraOffsetA, B, ExtraOffsetB, C, ExtraOffsetC, NumElements);
  CHECK_ERR;
}
__global__ void
    __launch_bounds__(96)
kernel_sgemm_T_NT_m9_9_n9_9_k9_ppp_db79476(const float ** A, int ExtraOffsetA, const float ** B, int ExtraOffsetB, float ** C, int ExtraOffsetC, unsigned NumElements) {
  if ((threadIdx.y + blockDim.y * blockIdx.x) < NumElements) {
    const float* GlobMatA = &A[(threadIdx.y + blockDim.y * blockIdx.x)][0 + ExtraOffsetA];
    const float* GlobMatB = &B[(threadIdx.y + blockDim.y * blockIdx.x)][0 + ExtraOffsetB];
    float* GlobMatC = &C[(threadIdx.y + blockDim.y * blockIdx.x)][0 + ExtraOffsetC];
    __shared__ float Scratch[540];
    float* ShrMatB = &Scratch[threadIdx.y * 180];
    float* ShrMatA = &Scratch[threadIdx.y * 180 + 81];

    // using ExtendedPatchLoader
    ShrMatB[threadIdx.x + 0] = GlobMatB[threadIdx.x + 0];
    ShrMatB[threadIdx.x + 32] = GlobMatB[threadIdx.x + 32];
    if (threadIdx.x < 17) {
      ShrMatB[threadIdx.x + 64] = GlobMatB[threadIdx.x + 64];
    }

    // using ExtendedTransposePatchLoader
    int Index;
    Index = threadIdx.x + 0;
    ShrMatA[(Index % 9) * 11 + Index / 9] = GlobMatA[threadIdx.x + 0];
    Index = threadIdx.x + 32;
    ShrMatA[(Index % 9) * 11 + Index / 9] = GlobMatA[threadIdx.x + 32];
    Index = threadIdx.x + 64;
    ShrMatA[(Index % 9) * 11 + Index / 9] = GlobMatA[threadIdx.x + 64];
    if (threadIdx.x < 3) {
      Index = threadIdx.x + 96;
      ShrMatA[(Index % 9) * 11 + Index / 9] = GlobMatA[threadIdx.x + 96];
    }
    __syncthreads();

    if (threadIdx.x < 9) {

      float Results[9] = {0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f};
      float Value;

      for (int k = 0; k < 9; ++k) {
        Value = ShrMatA[threadIdx.x + 11 * k];

#pragma unroll
        for (int n = 0; n < 9; ++n) {
          Results[n] += Value * ShrMatB[k + 9 * n];
        }
      }

#pragma unroll
      for (int n = 0; n < 9; ++n) {
        GlobMatC[threadIdx.x + 9 * n] = Results[n];
      }
    }
  }
}
void sgemm_T_NT_m9_9_n9_9_k9_ppp_db79476(const float ** A, int ExtraOffsetA, const float ** B, int ExtraOffsetB, float ** C, int ExtraOffsetC, unsigned NumElements) {
  dim3 Block(32, 3, 1);
  dim3 Grid((NumElements + 3 - 1) / 3, 1, 1);
  kernel_sgemm_T_NT_m9_9_n9_9_k9_ppp_db79476<<<Grid,Block>>>(A, ExtraOffsetA, B, ExtraOffsetB, C, ExtraOffsetC, NumElements);
  CHECK_ERR;
}
__global__ void
    __launch_bounds__(96)
kernel_sgemm_NT_T_m9_9_n9_9_k9_ppp_0c970c0(float fluxScale, const float ** A, int ExtraOffsetA, const float ** B, int ExtraOffsetB, float ** C, int ExtraOffsetC, unsigned NumElements) {
  if ((threadIdx.y + blockDim.y * blockIdx.x) < NumElements) {
    const float* GlobMatA = &A[(threadIdx.y + blockDim.y * blockIdx.x)][0 + ExtraOffsetA];
    const float* GlobMatB = &B[(threadIdx.y + blockDim.y * blockIdx.x)][0 + ExtraOffsetB];
    float* GlobMatC = &C[(threadIdx.y + blockDim.y * blockIdx.x)][0 + ExtraOffsetC];
    __shared__ float Scratch[243];
    float* ShrMatB = &Scratch[threadIdx.y * 81];

    // using ExtendedPatchLoader
    ShrMatB[threadIdx.x + 0] = GlobMatB[threadIdx.x + 0];
    ShrMatB[threadIdx.x + 32] = GlobMatB[threadIdx.x + 32];
    if (threadIdx.x < 17) {
      ShrMatB[threadIdx.x + 64] = GlobMatB[threadIdx.x + 64];
    }
    __syncthreads();

    if (threadIdx.x < 9) {

      float Results[9] = {0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f};
      float Value;

      for (int k = 0; k < 9; ++k) {
        Value = GlobMatA[threadIdx.x + 9 * k];

#pragma unroll
        for (int n = 0; n < 9; ++n) {
          Results[n] += Value * ShrMatB[n + 9 * k];
        }
      }

#pragma unroll
      for (int n = 0; n < 9; ++n) {
        GlobMatC[threadIdx.x + 9 * n] = fluxScale * Results[n];
      }
    }
  }
}
void sgemm_NT_T_m9_9_n9_9_k9_ppp_0c970c0(float fluxScale, const float ** A, int ExtraOffsetA, const float ** B, int ExtraOffsetB, float ** C, int ExtraOffsetC, unsigned NumElements) {
  dim3 Block(32, 3, 1);
  dim3 Grid((NumElements + 3 - 1) / 3, 1, 1);
  kernel_sgemm_NT_T_m9_9_n9_9_k9_ppp_0c970c0<<<Grid,Block>>>(fluxScale, A, ExtraOffsetA, B, ExtraOffsetB, C, ExtraOffsetC, NumElements);
  CHECK_ERR;
}
__global__ void
    __launch_bounds__(64)
kernel_sgemm_NT_NT_m56_56_n6_56_k56_npp_87cf7c3(const float * __restrict__ A, int ExtraOffsetA, const float ** B, int ExtraOffsetB, float ** C, int ExtraOffsetC, unsigned NumElements) {
  if ((threadIdx.y + blockDim.y * blockIdx.x) < NumElements) {
    const float* GlobMatA = &A[0 + ExtraOffsetA];
    const float* GlobMatB = &B[(threadIdx.y + blockDim.y * blockIdx.x)][0 + ExtraOffsetB];
    float* GlobMatC = &C[(threadIdx.y + blockDim.y * blockIdx.x)][0 + ExtraOffsetC];
    __shared__ float Scratch[336];
    float* ShrMatB = &Scratch[threadIdx.y * 336];

    // using ExtendedPatchLoader
#pragma unroll
    for (int i = 0; i < 5; ++i) {
      ShrMatB[threadIdx.x + i * 64] = GlobMatB[threadIdx.x + i * 64];
    }
    if (threadIdx.x < 16) {
      ShrMatB[threadIdx.x + 320] = GlobMatB[threadIdx.x + 320];
    }
    __syncthreads();

    if (threadIdx.x < 56) {

      float Results[6] = {0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f};
      float Value;

      for (int k = 0; k < 56; ++k) {
        Value = GlobMatA[threadIdx.x + 56 * k];

#pragma unroll
        for (int n = 0; n < 6; ++n) {
          Results[n] += Value * ShrMatB[k + 56 * n];
        }
      }

#pragma unroll
      for (int n = 0; n < 6; ++n) {
        GlobMatC[threadIdx.x + 56 * n] = Results[n];
      }
    }
  }
}
void sgemm_NT_NT_m56_56_n6_56_k56_npp_87cf7c3(const float * __restrict__ A, int ExtraOffsetA, const float ** B, int ExtraOffsetB, float ** C, int ExtraOffsetC, unsigned NumElements) {
  dim3 Block(64, 1, 1);
  dim3 Grid((NumElements + 1 - 1) / 1, 1, 1);
  kernel_sgemm_NT_NT_m56_56_n6_56_k56_npp_87cf7c3<<<Grid,Block>>>(A, ExtraOffsetA, B, ExtraOffsetB, C, ExtraOffsetC, NumElements);
  CHECK_ERR;
}
__global__ void
    __launch_bounds__(64)
kernel_sgemm_NT_NT_m56_56_n6_56_k56_npp_da40531(const float * __restrict__ A, int ExtraOffsetA, const float ** B, int ExtraOffsetB, float ** C, int ExtraOffsetC, unsigned NumElements) {
  if ((threadIdx.y + blockDim.y * blockIdx.x) < NumElements) {
    const float* GlobMatA = &A[0 + ExtraOffsetA];
    const float* GlobMatB = &B[(threadIdx.y + blockDim.y * blockIdx.x)][0 + ExtraOffsetB];
    float* GlobMatC = &C[(threadIdx.y + blockDim.y * blockIdx.x)][0 + ExtraOffsetC];
    __shared__ float Scratch[336];
    float* ShrMatB = &Scratch[threadIdx.y * 336];

    // using ExtendedPatchLoader
#pragma unroll
    for (int i = 0; i < 5; ++i) {
      ShrMatB[threadIdx.x + i * 64] = GlobMatB[threadIdx.x + i * 64];
    }
    if (threadIdx.x < 16) {
      ShrMatB[threadIdx.x + 320] = GlobMatB[threadIdx.x + 320];
    }
    __syncthreads();

    if (threadIdx.x < 56) {

      float Results[6] = {0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f};
      float Value;

      for (int k = 0; k < 56; ++k) {
        Value = GlobMatA[threadIdx.x + 56 * k];

#pragma unroll
        for (int n = 0; n < 6; ++n) {
          Results[n] += Value * ShrMatB[k + 56 * n];
        }
      }

#pragma unroll
      for (int n = 0; n < 6; ++n) {
        GlobMatC[threadIdx.x + 56 * n] = Results[n] + GlobMatC[threadIdx.x + 56 * n];
      }
    }
  }
}
void sgemm_NT_NT_m56_56_n6_56_k56_npp_da40531(const float * __restrict__ A, int ExtraOffsetA, const float ** B, int ExtraOffsetB, float ** C, int ExtraOffsetC, unsigned NumElements) {
  dim3 Block(64, 1, 1);
  dim3 Grid((NumElements + 1 - 1) / 1, 1, 1);
  kernel_sgemm_NT_NT_m56_56_n6_56_k56_npp_da40531<<<Grid,Block>>>(A, ExtraOffsetA, B, ExtraOffsetB, C, ExtraOffsetC, NumElements);
  CHECK_ERR;
}
__global__ void
    __launch_bounds__(32)
kernel_sgemm_NT_NT_m21_24_n9_56_k56_nps_dd17eae(const float * __restrict__ A, int ExtraOffsetA, const float ** B, int ExtraOffsetB, float * C, int ExtraOffsetC, unsigned NumElements) {
  if ((threadIdx.y + blockDim.y * blockIdx.x) < NumElements) {
    const float* GlobMatA = &A[0 + ExtraOffsetA];
    const float* GlobMatB = &B[(threadIdx.y + blockDim.y * blockIdx.x)][0 + ExtraOffsetB];
    float* GlobMatC = &C[(threadIdx.y + blockDim.y * blockIdx.x) * 216 + 0 + ExtraOffsetC];
    __shared__ float Scratch[504];
    float* ShrMatB = &Scratch[threadIdx.y * 504];

    // using ExtendedPatchLoader
#pragma unroll
    for (int i = 0; i < 15; ++i) {
      ShrMatB[threadIdx.x + i * 32] = GlobMatB[threadIdx.x + i * 32];
    }
    if (threadIdx.x < 24) {
      ShrMatB[threadIdx.x + 480] = GlobMatB[threadIdx.x + 480];
    }
    __syncthreads();

    if (threadIdx.x < 21) {

      float Results[9] = {0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f};
      float Value;

      for (int k = 0; k < 56; ++k) {
        Value = GlobMatA[threadIdx.x + 24 * k];

#pragma unroll
        for (int n = 0; n < 9; ++n) {
          Results[n] += Value * ShrMatB[k + 56 * n];
        }
      }

#pragma unroll
      for (int n = 0; n < 9; ++n) {
        GlobMatC[threadIdx.x + 24 * n] = Results[n];
      }
    }
  }
}
void sgemm_NT_NT_m21_24_n9_56_k56_nps_dd17eae(const float * __restrict__ A, int ExtraOffsetA, const float ** B, int ExtraOffsetB, float * C, int ExtraOffsetC, unsigned NumElements) {
  dim3 Block(32, 1, 1);
  dim3 Grid((NumElements + 1 - 1) / 1, 1, 1);
  kernel_sgemm_NT_NT_m21_24_n9_56_k56_nps_dd17eae<<<Grid,Block>>>(A, ExtraOffsetA, B, ExtraOffsetB, C, ExtraOffsetC, NumElements);
  CHECK_ERR;
}
__global__ void
    __launch_bounds__(64)
kernel_sgemm_NT_NT_m21_24_n9_9_k9_sps_899f278(const float * A, int ExtraOffsetA, const float ** B, int ExtraOffsetB, float * C, int ExtraOffsetC, unsigned NumElements) {
  if ((threadIdx.y + blockDim.y * blockIdx.x) < NumElements) {
    const float* GlobMatA = &A[(threadIdx.y + blockDim.y * blockIdx.x) * 216 + 0 + ExtraOffsetA];
    const float* GlobMatB = &B[(threadIdx.y + blockDim.y * blockIdx.x)][0 + ExtraOffsetB];
    float* GlobMatC = &C[(threadIdx.y + blockDim.y * blockIdx.x) * 216 + 0 + ExtraOffsetC];
    __shared__ float Scratch[162];
    float* ShrMatB = &Scratch[threadIdx.y * 81];

    // using ExtendedPatchLoader
    ShrMatB[threadIdx.x + 0] = GlobMatB[threadIdx.x + 0];
    ShrMatB[threadIdx.x + 32] = GlobMatB[threadIdx.x + 32];
    if (threadIdx.x < 17) {
      ShrMatB[threadIdx.x + 64] = GlobMatB[threadIdx.x + 64];
    }
    __syncthreads();

    if (threadIdx.x < 21) {

      float Results[9] = {0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f};
      float Value;

      for (int k = 0; k < 9; ++k) {
        Value = GlobMatA[threadIdx.x + 24 * k];

#pragma unroll
        for (int n = 0; n < 9; ++n) {
          Results[n] += Value * ShrMatB[k + 9 * n];
        }
      }

#pragma unroll
      for (int n = 0; n < 9; ++n) {
        GlobMatC[threadIdx.x + 24 * n] = Results[n];
      }
    }
  }
}
void sgemm_NT_NT_m21_24_n9_9_k9_sps_899f278(const float * A, int ExtraOffsetA, const float ** B, int ExtraOffsetB, float * C, int ExtraOffsetC, unsigned NumElements) {
  dim3 Block(32, 2, 1);
  dim3 Grid((NumElements + 2 - 1) / 2, 1, 1);
  kernel_sgemm_NT_NT_m21_24_n9_9_k9_sps_899f278<<<Grid,Block>>>(A, ExtraOffsetA, B, ExtraOffsetB, C, ExtraOffsetC, NumElements);
  CHECK_ERR;
}
__global__ void
    __launch_bounds__(64)
kernel_sgemm_NT_NT_m56_56_n9_24_k21_nsp_41f2e57(const float * __restrict__ A, int ExtraOffsetA, const float * B, int ExtraOffsetB, float ** C, int ExtraOffsetC, unsigned NumElements) {
  if ((threadIdx.y + blockDim.y * blockIdx.x) < NumElements) {
    const float* GlobMatA = &A[0 + ExtraOffsetA];
    const float* GlobMatB = &B[(threadIdx.y + blockDim.y * blockIdx.x) * 216 + 0 + ExtraOffsetB];
    float* GlobMatC = &C[(threadIdx.y + blockDim.y * blockIdx.x)][0 + ExtraOffsetC];
    __shared__ float Scratch[213];
    float* ShrMatB = &Scratch[threadIdx.y * 213];

    // using ExtendedPatchLoader
    ShrMatB[threadIdx.x + 0] = GlobMatB[threadIdx.x + 0];
    ShrMatB[threadIdx.x + 64] = GlobMatB[threadIdx.x + 64];
    ShrMatB[threadIdx.x + 128] = GlobMatB[threadIdx.x + 128];
    if (threadIdx.x < 21) {
      ShrMatB[threadIdx.x + 192] = GlobMatB[threadIdx.x + 192];
    }
    __syncthreads();

    if (threadIdx.x < 56) {

      float Results[9] = {0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f};
      float Value;

      for (int k = 0; k < 21; ++k) {
        Value = GlobMatA[threadIdx.x + 56 * k];

#pragma unroll
        for (int n = 0; n < 9; ++n) {
          Results[n] += Value * ShrMatB[k + 24 * n];
        }
      }

#pragma unroll
      for (int n = 0; n < 9; ++n) {
        GlobMatC[threadIdx.x + 56 * n] = Results[n] + GlobMatC[threadIdx.x + 56 * n];
      }
    }
  }
}
void sgemm_NT_NT_m56_56_n9_24_k21_nsp_41f2e57(const float * __restrict__ A, int ExtraOffsetA, const float * B, int ExtraOffsetB, float ** C, int ExtraOffsetC, unsigned NumElements) {
  dim3 Block(64, 1, 1);
  dim3 Grid((NumElements + 1 - 1) / 1, 1, 1);
  kernel_sgemm_NT_NT_m56_56_n9_24_k21_nsp_41f2e57<<<Grid,Block>>>(A, ExtraOffsetA, B, ExtraOffsetB, C, ExtraOffsetC, NumElements);
  CHECK_ERR;
}
__global__ void
    __launch_bounds__(64)
kernel_sgemm_NT_NT_m21_24_n9_24_k21_nss_5830167(const float * __restrict__ A, int ExtraOffsetA, const float * B, int ExtraOffsetB, float * C, int ExtraOffsetC, unsigned NumElements) {
  if ((threadIdx.y + blockDim.y * blockIdx.x) < NumElements) {
    const float* GlobMatA = &A[0 + ExtraOffsetA];
    const float* GlobMatB = &B[(threadIdx.y + blockDim.y * blockIdx.x) * 216 + 0 + ExtraOffsetB];
    float* GlobMatC = &C[(threadIdx.y + blockDim.y * blockIdx.x) * 216 + 0 + ExtraOffsetC];
    __shared__ float Scratch[426];
    float* ShrMatB = &Scratch[threadIdx.y * 213];

    // using ExtendedPatchLoader
#pragma unroll
    for (int i = 0; i < 6; ++i) {
      ShrMatB[threadIdx.x + i * 32] = GlobMatB[threadIdx.x + i * 32];
    }
    if (threadIdx.x < 21) {
      ShrMatB[threadIdx.x + 192] = GlobMatB[threadIdx.x + 192];
    }
    __syncthreads();

    if (threadIdx.x < 21) {

      float Results[9] = {0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f};
      float Value;

      for (int k = 0; k < 21; ++k) {
        Value = GlobMatA[threadIdx.x + 24 * k];

#pragma unroll
        for (int n = 0; n < 9; ++n) {
          Results[n] += Value * ShrMatB[k + 24 * n];
        }
      }

#pragma unroll
      for (int n = 0; n < 9; ++n) {
        GlobMatC[threadIdx.x + 24 * n] = Results[n];
      }
    }
  }
}
void sgemm_NT_NT_m21_24_n9_24_k21_nss_5830167(const float * __restrict__ A, int ExtraOffsetA, const float * B, int ExtraOffsetB, float * C, int ExtraOffsetC, unsigned NumElements) {
  dim3 Block(32, 2, 1);
  dim3 Grid((NumElements + 2 - 1) / 2, 1, 1);
  kernel_sgemm_NT_NT_m21_24_n9_24_k21_nss_5830167<<<Grid,Block>>>(A, ExtraOffsetA, B, ExtraOffsetB, C, ExtraOffsetC, NumElements);
  CHECK_ERR;
}
__global__ void
    __launch_bounds__(576)
kernel_scopyAddScale_m56_56_n9_56_pp_88237f4(float Scale, const float ** A, int ExtraOffsetA, float ** B, int ExtraOffsetB, unsigned NumElements) {
  if ((threadIdx.z + blockDim.z * blockIdx.x) < NumElements) {
    const float* GlobMatA = &A[(threadIdx.z + blockDim.z * blockIdx.x)][ExtraOffsetA + 0 + threadIdx.y * 56];
    float* GlobMatB = &B[(threadIdx.z + blockDim.z * blockIdx.x)][ExtraOffsetB + 0 + threadIdx.y * 56];
    if (threadIdx.x < 56) {
      GlobMatB[threadIdx.x] = Scale * GlobMatA[threadIdx.x];
    }
  }
}
__global__ void
    __launch_bounds__(576)
kernel_initialize_m56_56_p_8eacdf4(float ** A, int ExtraOffsetA, unsigned NumElements) {
  if ((threadIdx.z + blockDim.z * blockIdx.x) < NumElements) {
    float* GlobA = &A[(threadIdx.z + blockDim.z * blockIdx.x)][ExtraOffsetA + 0 + threadIdx.y * 56];
    if (threadIdx.x < 56) {
      GlobA[threadIdx.x] = 0.0;
    }
  }
}
void initialize_m56_56_p_8eacdf4(float ** A, int ExtraOffsetA, unsigned NumElements) {
  dim3 Block(64, 9, 1);
  dim3 Grid((NumElements + 1 - 1) / 1, 1, 1);
  kernel_initialize_m56_56_p_8eacdf4<<<Grid,Block>>>(A, ExtraOffsetA, NumElements);
  CHECK_ERR;
}
void scopyAddScale_m56_56_n9_56_pp_88237f4(float Scale, const float ** A, int ExtraOffsetA, float ** B, int ExtraOffsetB, unsigned NumElements) {
  initialize_m56_56_p_8eacdf4(B, ExtraOffsetB, NumElements);
  dim3 Block(64, 9, 1);
  dim3 Grid((NumElements + 1 - 1) / 1, 1, 1);
  kernel_scopyAddScale_m56_56_n9_56_pp_88237f4<<<Grid,Block>>>(Scale, A, ExtraOffsetA, B, ExtraOffsetB, NumElements);
  CHECK_ERR;
}
__global__ void
    __launch_bounds__(576)
kernel_scopyAddScale_m35_56_n9_56_pp_07aaceb(float Scale, const float ** A, int ExtraOffsetA, float ** B, int ExtraOffsetB, unsigned NumElements) {
  if ((threadIdx.z + blockDim.z * blockIdx.x) < NumElements) {
    const float* GlobMatA = &A[(threadIdx.z + blockDim.z * blockIdx.x)][ExtraOffsetA + 0 + threadIdx.y * 40];
    float* GlobMatB = &B[(threadIdx.z + blockDim.z * blockIdx.x)][ExtraOffsetB + 0 + threadIdx.y * 56];
    if (threadIdx.x < 35) {
      GlobMatB[threadIdx.x] += Scale * GlobMatA[threadIdx.x];
    }
  }
}
void scopyAddScale_m35_56_n9_56_pp_07aaceb(float Scale, const float ** A, int ExtraOffsetA, float ** B, int ExtraOffsetB, unsigned NumElements) {
  dim3 Block(64, 9, 1);
  dim3 Grid((NumElements + 1 - 1) / 1, 1, 1);
  kernel_scopyAddScale_m35_56_n9_56_pp_07aaceb<<<Grid,Block>>>(Scale, A, ExtraOffsetA, B, ExtraOffsetB, NumElements);
  CHECK_ERR;
}
__global__ void
    __launch_bounds__(288)
kernel_scopyAddScale_m20_56_n9_56_pp_6eb2445(float Scale, const float ** A, int ExtraOffsetA, float ** B, int ExtraOffsetB, unsigned NumElements) {
  if ((threadIdx.z + blockDim.z * blockIdx.x) < NumElements) {
    const float* GlobMatA = &A[(threadIdx.z + blockDim.z * blockIdx.x)][ExtraOffsetA + 0 + threadIdx.y * 24];
    float* GlobMatB = &B[(threadIdx.z + blockDim.z * blockIdx.x)][ExtraOffsetB + 0 + threadIdx.y * 56];
    if (threadIdx.x < 20) {
      GlobMatB[threadIdx.x] += Scale * GlobMatA[threadIdx.x];
    }
  }
}
void scopyAddScale_m20_56_n9_56_pp_6eb2445(float Scale, const float ** A, int ExtraOffsetA, float ** B, int ExtraOffsetB, unsigned NumElements) {
  dim3 Block(32, 9, 1);
  dim3 Grid((NumElements + 1 - 1) / 1, 1, 1);
  kernel_scopyAddScale_m20_56_n9_56_pp_6eb2445<<<Grid,Block>>>(Scale, A, ExtraOffsetA, B, ExtraOffsetB, NumElements);
  CHECK_ERR;
}
__global__ void
    __launch_bounds__(288)
kernel_scopyAddScale_m10_56_n9_56_pp_afcc802(float Scale, const float ** A, int ExtraOffsetA, float ** B, int ExtraOffsetB, unsigned NumElements) {
  if ((threadIdx.z + blockDim.z * blockIdx.x) < NumElements) {
    const float* GlobMatA = &A[(threadIdx.z + blockDim.z * blockIdx.x)][ExtraOffsetA + 0 + threadIdx.y * 16];
    float* GlobMatB = &B[(threadIdx.z + blockDim.z * blockIdx.x)][ExtraOffsetB + 0 + threadIdx.y * 56];
    if (threadIdx.x < 10) {
      GlobMatB[threadIdx.x] += Scale * GlobMatA[threadIdx.x];
    }
  }
}
void scopyAddScale_m10_56_n9_56_pp_afcc802(float Scale, const float ** A, int ExtraOffsetA, float ** B, int ExtraOffsetB, unsigned NumElements) {
  dim3 Block(32, 9, 1);
  dim3 Grid((NumElements + 1 - 1) / 1, 1, 1);
  kernel_scopyAddScale_m10_56_n9_56_pp_afcc802<<<Grid,Block>>>(Scale, A, ExtraOffsetA, B, ExtraOffsetB, NumElements);
  CHECK_ERR;
}
__global__ void
    __launch_bounds__(288)
kernel_scopyAddScale_m4_56_n9_56_pp_c707c5f(float Scale, const float ** A, int ExtraOffsetA, float ** B, int ExtraOffsetB, unsigned NumElements) {
  if ((threadIdx.z + blockDim.z * blockIdx.x) < NumElements) {
    const float* GlobMatA = &A[(threadIdx.z + blockDim.z * blockIdx.x)][ExtraOffsetA + 0 + threadIdx.y * 8];
    float* GlobMatB = &B[(threadIdx.z + blockDim.z * blockIdx.x)][ExtraOffsetB + 0 + threadIdx.y * 56];
    if (threadIdx.x < 4) {
      GlobMatB[threadIdx.x] += Scale * GlobMatA[threadIdx.x];
    }
  }
}
void scopyAddScale_m4_56_n9_56_pp_c707c5f(float Scale, const float ** A, int ExtraOffsetA, float ** B, int ExtraOffsetB, unsigned NumElements) {
  dim3 Block(32, 9, 1);
  dim3 Grid((NumElements + 1 - 1) / 1, 1, 1);
  kernel_scopyAddScale_m4_56_n9_56_pp_c707c5f<<<Grid,Block>>>(Scale, A, ExtraOffsetA, B, ExtraOffsetB, NumElements);
  CHECK_ERR;
}
__global__ void
    __launch_bounds__(288)
kernel_scopyAddScale_m1_56_n9_56_pp_96cca66(float Scale, const float ** A, int ExtraOffsetA, float ** B, int ExtraOffsetB, unsigned NumElements) {
  if ((threadIdx.z + blockDim.z * blockIdx.x) < NumElements) {
    const float* GlobMatA = &A[(threadIdx.z + blockDim.z * blockIdx.x)][ExtraOffsetA + 0 + threadIdx.y * 8];
    float* GlobMatB = &B[(threadIdx.z + blockDim.z * blockIdx.x)][ExtraOffsetB + 0 + threadIdx.y * 56];
    if (threadIdx.x < 1) {
      GlobMatB[threadIdx.x] += Scale * GlobMatA[threadIdx.x];
    }
  }
}
void scopyAddScale_m1_56_n9_56_pp_96cca66(float Scale, const float ** A, int ExtraOffsetA, float ** B, int ExtraOffsetB, unsigned NumElements) {
  dim3 Block(32, 9, 1);
  dim3 Grid((NumElements + 1 - 1) / 1, 1, 1);
  kernel_scopyAddScale_m1_56_n9_56_pp_96cca66<<<Grid,Block>>>(Scale, A, ExtraOffsetA, B, ExtraOffsetB, NumElements);
  CHECK_ERR;
}
__global__ void
    __launch_bounds__(64)
kernel_sgemm_NT_NT_m35_40_n9_56_k53_nps_960fd76(const float * __restrict__ A, int ExtraOffsetA, const float ** B, int ExtraOffsetB, float * C, int ExtraOffsetC, unsigned NumElements) {
  if ((threadIdx.y + blockDim.y * blockIdx.x) < NumElements) {
    const float* GlobMatA = &A[40 + ExtraOffsetA];
    const float* GlobMatB = &B[(threadIdx.y + blockDim.y * blockIdx.x)][1 + ExtraOffsetB];
    float* GlobMatC = &C[(threadIdx.y + blockDim.y * blockIdx.x) * 360 + 0 + ExtraOffsetC];
    __shared__ float Scratch[501];
    float* ShrMatB = &Scratch[threadIdx.y * 501];

    // using ExtendedPatchLoader
#pragma unroll
    for (int i = 0; i < 7; ++i) {
      ShrMatB[threadIdx.x + i * 64] = GlobMatB[threadIdx.x + i * 64];
    }
    if (threadIdx.x < 53) {
      ShrMatB[threadIdx.x + 448] = GlobMatB[threadIdx.x + 448];
    }
    __syncthreads();

    if (threadIdx.x < 35) {

      float Results[9] = {0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f};
      float Value;

      for (int k = 0; k < 53; ++k) {
        Value = GlobMatA[threadIdx.x + 40 * k];

#pragma unroll
        for (int n = 0; n < 9; ++n) {
          Results[n] += Value * ShrMatB[k + 56 * n];
        }
      }

#pragma unroll
      for (int n = 0; n < 9; ++n) {
        GlobMatC[threadIdx.x + 40 * n] = Results[n];
      }
    }
  }
}
void sgemm_NT_NT_m35_40_n9_56_k53_nps_960fd76(const float * __restrict__ A, int ExtraOffsetA, const float ** B, int ExtraOffsetB, float * C, int ExtraOffsetC, unsigned NumElements) {
  dim3 Block(64, 1, 1);
  dim3 Grid((NumElements + 1 - 1) / 1, 1, 1);
  kernel_sgemm_NT_NT_m35_40_n9_56_k53_nps_960fd76<<<Grid,Block>>>(A, ExtraOffsetA, B, ExtraOffsetB, C, ExtraOffsetC, NumElements);
  CHECK_ERR;
}
__global__ void
    __launch_bounds__(64)
kernel_sgemm_NT_NT_m35_40_n9_9_k9_spp_36989d0(const float * A, int ExtraOffsetA, const float ** B, int ExtraOffsetB, float ** C, int ExtraOffsetC, unsigned NumElements) {
  if ((threadIdx.y + blockDim.y * blockIdx.x) < NumElements) {
    const float* GlobMatA = &A[(threadIdx.y + blockDim.y * blockIdx.x) * 360 + 0 + ExtraOffsetA];
    const float* GlobMatB = &B[(threadIdx.y + blockDim.y * blockIdx.x)][0 + ExtraOffsetB];
    float* GlobMatC = &C[(threadIdx.y + blockDim.y * blockIdx.x)][0 + ExtraOffsetC];
    __shared__ float Scratch[81];
    float* ShrMatB = &Scratch[threadIdx.y * 81];

    // using ExtendedPatchLoader
    ShrMatB[threadIdx.x + 0] = GlobMatB[threadIdx.x + 0];
    if (threadIdx.x < 17) {
      ShrMatB[threadIdx.x + 64] = GlobMatB[threadIdx.x + 64];
    }
    __syncthreads();

    if (threadIdx.x < 35) {

      float Results[9] = {0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f};
      float Value;

      for (int k = 0; k < 9; ++k) {
        Value = GlobMatA[threadIdx.x + 40 * k];

#pragma unroll
        for (int n = 0; n < 9; ++n) {
          Results[n] += Value * ShrMatB[k + 9 * n];
        }
      }

#pragma unroll
      for (int n = 0; n < 9; ++n) {
        GlobMatC[threadIdx.x + 40 * n] = Results[n];
      }
    }
  }
}
void sgemm_NT_NT_m35_40_n9_9_k9_spp_36989d0(const float * A, int ExtraOffsetA, const float ** B, int ExtraOffsetB, float ** C, int ExtraOffsetC, unsigned NumElements) {
  dim3 Block(64, 1, 1);
  dim3 Grid((NumElements + 1 - 1) / 1, 1, 1);
  kernel_sgemm_NT_NT_m35_40_n9_9_k9_spp_36989d0<<<Grid,Block>>>(A, ExtraOffsetA, B, ExtraOffsetB, C, ExtraOffsetC, NumElements);
  CHECK_ERR;
}
__global__ void
    __launch_bounds__(64)
kernel_sgemm_NT_NT_m35_40_n9_56_k54_nps_2269f5e(const float * __restrict__ A, int ExtraOffsetA, const float ** B, int ExtraOffsetB, float * C, int ExtraOffsetC, unsigned NumElements) {
  if ((threadIdx.y + blockDim.y * blockIdx.x) < NumElements) {
    const float* GlobMatA = &A[40 + ExtraOffsetA];
    const float* GlobMatB = &B[(threadIdx.y + blockDim.y * blockIdx.x)][1 + ExtraOffsetB];
    float* GlobMatC = &C[(threadIdx.y + blockDim.y * blockIdx.x) * 360 + 0 + ExtraOffsetC];
    __shared__ float Scratch[502];
    float* ShrMatB = &Scratch[threadIdx.y * 502];

    // using ExtendedPatchLoader
#pragma unroll
    for (int i = 0; i < 7; ++i) {
      ShrMatB[threadIdx.x + i * 64] = GlobMatB[threadIdx.x + i * 64];
    }
    if (threadIdx.x < 54) {
      ShrMatB[threadIdx.x + 448] = GlobMatB[threadIdx.x + 448];
    }
    __syncthreads();

    if (threadIdx.x < 35) {

      float Results[9] = {0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f};
      float Value;

      for (int k = 0; k < 54; ++k) {
        Value = GlobMatA[threadIdx.x + 40 * k];

#pragma unroll
        for (int n = 0; n < 9; ++n) {
          Results[n] += Value * ShrMatB[k + 56 * n];
        }
      }

#pragma unroll
      for (int n = 0; n < 9; ++n) {
        GlobMatC[threadIdx.x + 40 * n] = Results[n];
      }
    }
  }
}
void sgemm_NT_NT_m35_40_n9_56_k54_nps_2269f5e(const float * __restrict__ A, int ExtraOffsetA, const float ** B, int ExtraOffsetB, float * C, int ExtraOffsetC, unsigned NumElements) {
  dim3 Block(64, 1, 1);
  dim3 Grid((NumElements + 1 - 1) / 1, 1, 1);
  kernel_sgemm_NT_NT_m35_40_n9_56_k54_nps_2269f5e<<<Grid,Block>>>(A, ExtraOffsetA, B, ExtraOffsetB, C, ExtraOffsetC, NumElements);
  CHECK_ERR;
}
__global__ void
    __launch_bounds__(64)
kernel_sgemm_NT_NT_m35_40_n9_9_k9_spp_9d7d789(const float * A, int ExtraOffsetA, const float ** B, int ExtraOffsetB, float ** C, int ExtraOffsetC, unsigned NumElements) {
  if ((threadIdx.y + blockDim.y * blockIdx.x) < NumElements) {
    const float* GlobMatA = &A[(threadIdx.y + blockDim.y * blockIdx.x) * 360 + 0 + ExtraOffsetA];
    const float* GlobMatB = &B[(threadIdx.y + blockDim.y * blockIdx.x)][0 + ExtraOffsetB];
    float* GlobMatC = &C[(threadIdx.y + blockDim.y * blockIdx.x)][0 + ExtraOffsetC];
    __shared__ float Scratch[81];
    float* ShrMatB = &Scratch[threadIdx.y * 81];

    // using ExtendedPatchLoader
    ShrMatB[threadIdx.x + 0] = GlobMatB[threadIdx.x + 0];
    if (threadIdx.x < 17) {
      ShrMatB[threadIdx.x + 64] = GlobMatB[threadIdx.x + 64];
    }
    __syncthreads();

    if (threadIdx.x < 35) {

      float Results[9] = {0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f};
      float Value;

      for (int k = 0; k < 9; ++k) {
        Value = GlobMatA[threadIdx.x + 40 * k];

#pragma unroll
        for (int n = 0; n < 9; ++n) {
          Results[n] += Value * ShrMatB[k + 9 * n];
        }
      }

#pragma unroll
      for (int n = 0; n < 9; ++n) {
        GlobMatC[threadIdx.x + 40 * n] = Results[n] + GlobMatC[threadIdx.x + 40 * n];
      }
    }
  }
}
void sgemm_NT_NT_m35_40_n9_9_k9_spp_9d7d789(const float * A, int ExtraOffsetA, const float ** B, int ExtraOffsetB, float ** C, int ExtraOffsetC, unsigned NumElements) {
  dim3 Block(64, 1, 1);
  dim3 Grid((NumElements + 1 - 1) / 1, 1, 1);
  kernel_sgemm_NT_NT_m35_40_n9_9_k9_spp_9d7d789<<<Grid,Block>>>(A, ExtraOffsetA, B, ExtraOffsetB, C, ExtraOffsetC, NumElements);
  CHECK_ERR;
}
__global__ void
    __launch_bounds__(64)
kernel_sgemm_NT_NT_m35_40_n9_56_k55_nps_5d71d3c(const float * __restrict__ A, int ExtraOffsetA, const float ** B, int ExtraOffsetB, float * C, int ExtraOffsetC, unsigned NumElements) {
  if ((threadIdx.y + blockDim.y * blockIdx.x) < NumElements) {
    const float* GlobMatA = &A[40 + ExtraOffsetA];
    const float* GlobMatB = &B[(threadIdx.y + blockDim.y * blockIdx.x)][1 + ExtraOffsetB];
    float* GlobMatC = &C[(threadIdx.y + blockDim.y * blockIdx.x) * 360 + 0 + ExtraOffsetC];
    __shared__ float Scratch[503];
    float* ShrMatB = &Scratch[threadIdx.y * 503];

    // using ExtendedPatchLoader
#pragma unroll
    for (int i = 0; i < 7; ++i) {
      ShrMatB[threadIdx.x + i * 64] = GlobMatB[threadIdx.x + i * 64];
    }
    if (threadIdx.x < 55) {
      ShrMatB[threadIdx.x + 448] = GlobMatB[threadIdx.x + 448];
    }
    __syncthreads();

    if (threadIdx.x < 35) {

      float Results[9] = {0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f};
      float Value;

      for (int k = 0; k < 55; ++k) {
        Value = GlobMatA[threadIdx.x + 40 * k];

#pragma unroll
        for (int n = 0; n < 9; ++n) {
          Results[n] += Value * ShrMatB[k + 56 * n];
        }
      }

#pragma unroll
      for (int n = 0; n < 9; ++n) {
        GlobMatC[threadIdx.x + 40 * n] = Results[n];
      }
    }
  }
}
void sgemm_NT_NT_m35_40_n9_56_k55_nps_5d71d3c(const float * __restrict__ A, int ExtraOffsetA, const float ** B, int ExtraOffsetB, float * C, int ExtraOffsetC, unsigned NumElements) {
  dim3 Block(64, 1, 1);
  dim3 Grid((NumElements + 1 - 1) / 1, 1, 1);
  kernel_sgemm_NT_NT_m35_40_n9_56_k55_nps_5d71d3c<<<Grid,Block>>>(A, ExtraOffsetA, B, ExtraOffsetB, C, ExtraOffsetC, NumElements);
  CHECK_ERR;
}
__global__ void
    __launch_bounds__(64)
kernel_sgemm_NT_NT_m20_40_n9_40_k32_nps_c33eb71(const float * __restrict__ A, int ExtraOffsetA, const float ** B, int ExtraOffsetB, float * C, int ExtraOffsetC, unsigned NumElements) {
  if ((threadIdx.y + blockDim.y * blockIdx.x) < NumElements) {
    const float* GlobMatA = &A[40 + ExtraOffsetA];
    const float* GlobMatB = &B[(threadIdx.y + blockDim.y * blockIdx.x)][1 + ExtraOffsetB];
    float* GlobMatC = &C[(threadIdx.y + blockDim.y * blockIdx.x) * 216 + 0 + ExtraOffsetC];
    __shared__ float Scratch[576];
    float* ShrMatB = &Scratch[threadIdx.y * 288];

    // using ExactPatchLoader
#pragma unroll
    for (int i = 0; i < 9; ++i) {
      ShrMatB[threadIdx.x + 0 + i * 32] = GlobMatB[threadIdx.x + 0 + i * 40];
    }
    __syncthreads();

    if (threadIdx.x < 20) {

      float Results[9] = {0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f};
      float Value;

      for (int k = 0; k < 32; ++k) {
        Value = GlobMatA[threadIdx.x + 40 * k];

#pragma unroll
        for (int n = 0; n < 9; ++n) {
          Results[n] += Value * ShrMatB[k + 32 * n];
        }
      }

#pragma unroll
      for (int n = 0; n < 9; ++n) {
        GlobMatC[threadIdx.x + 24 * n] = Results[n];
      }
    }
  }
}
void sgemm_NT_NT_m20_40_n9_40_k32_nps_c33eb71(const float * __restrict__ A, int ExtraOffsetA, const float ** B, int ExtraOffsetB, float * C, int ExtraOffsetC, unsigned NumElements) {
  dim3 Block(32, 2, 1);
  dim3 Grid((NumElements + 2 - 1) / 2, 1, 1);
  kernel_sgemm_NT_NT_m20_40_n9_40_k32_nps_c33eb71<<<Grid,Block>>>(A, ExtraOffsetA, B, ExtraOffsetB, C, ExtraOffsetC, NumElements);
  CHECK_ERR;
}
__global__ void
    __launch_bounds__(64)
kernel_sgemm_NT_NT_m20_24_n9_9_k9_spp_0ba88dd(const float * A, int ExtraOffsetA, const float ** B, int ExtraOffsetB, float ** C, int ExtraOffsetC, unsigned NumElements) {
  if ((threadIdx.y + blockDim.y * blockIdx.x) < NumElements) {
    const float* GlobMatA = &A[(threadIdx.y + blockDim.y * blockIdx.x) * 216 + 0 + ExtraOffsetA];
    const float* GlobMatB = &B[(threadIdx.y + blockDim.y * blockIdx.x)][0 + ExtraOffsetB];
    float* GlobMatC = &C[(threadIdx.y + blockDim.y * blockIdx.x)][0 + ExtraOffsetC];
    __shared__ float Scratch[162];
    float* ShrMatB = &Scratch[threadIdx.y * 81];

    // using ExtendedPatchLoader
    ShrMatB[threadIdx.x + 0] = GlobMatB[threadIdx.x + 0];
    ShrMatB[threadIdx.x + 32] = GlobMatB[threadIdx.x + 32];
    if (threadIdx.x < 17) {
      ShrMatB[threadIdx.x + 64] = GlobMatB[threadIdx.x + 64];
    }
    __syncthreads();

    if (threadIdx.x < 20) {

      float Results[9] = {0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f};
      float Value;

      for (int k = 0; k < 9; ++k) {
        Value = GlobMatA[threadIdx.x + 24 * k];

#pragma unroll
        for (int n = 0; n < 9; ++n) {
          Results[n] += Value * ShrMatB[k + 9 * n];
        }
      }

#pragma unroll
      for (int n = 0; n < 9; ++n) {
        GlobMatC[threadIdx.x + 24 * n] = Results[n];
      }
    }
  }
}
void sgemm_NT_NT_m20_24_n9_9_k9_spp_0ba88dd(const float * A, int ExtraOffsetA, const float ** B, int ExtraOffsetB, float ** C, int ExtraOffsetC, unsigned NumElements) {
  dim3 Block(32, 2, 1);
  dim3 Grid((NumElements + 2 - 1) / 2, 1, 1);
  kernel_sgemm_NT_NT_m20_24_n9_9_k9_spp_0ba88dd<<<Grid,Block>>>(A, ExtraOffsetA, B, ExtraOffsetB, C, ExtraOffsetC, NumElements);
  CHECK_ERR;
}
__global__ void
    __launch_bounds__(64)
kernel_sgemm_NT_NT_m20_40_n9_40_k33_nps_c9310c2(const float * __restrict__ A, int ExtraOffsetA, const float ** B, int ExtraOffsetB, float * C, int ExtraOffsetC, unsigned NumElements) {
  if ((threadIdx.y + blockDim.y * blockIdx.x) < NumElements) {
    const float* GlobMatA = &A[40 + ExtraOffsetA];
    const float* GlobMatB = &B[(threadIdx.y + blockDim.y * blockIdx.x)][1 + ExtraOffsetB];
    float* GlobMatC = &C[(threadIdx.y + blockDim.y * blockIdx.x) * 216 + 0 + ExtraOffsetC];
    __shared__ float Scratch[706];
    float* ShrMatB = &Scratch[threadIdx.y * 353];

    // using ExtendedPatchLoader
#pragma unroll
    for (int i = 0; i < 11; ++i) {
      ShrMatB[threadIdx.x + i * 32] = GlobMatB[threadIdx.x + i * 32];
    }
    if (threadIdx.x < 1) {
      ShrMatB[threadIdx.x + 352] = GlobMatB[threadIdx.x + 352];
    }
    __syncthreads();

    if (threadIdx.x < 20) {

      float Results[9] = {0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f};
      float Value;

      for (int k = 0; k < 33; ++k) {
        Value = GlobMatA[threadIdx.x + 40 * k];

#pragma unroll
        for (int n = 0; n < 9; ++n) {
          Results[n] += Value * ShrMatB[k + 40 * n];
        }
      }

#pragma unroll
      for (int n = 0; n < 9; ++n) {
        GlobMatC[threadIdx.x + 24 * n] = Results[n];
      }
    }
  }
}
void sgemm_NT_NT_m20_40_n9_40_k33_nps_c9310c2(const float * __restrict__ A, int ExtraOffsetA, const float ** B, int ExtraOffsetB, float * C, int ExtraOffsetC, unsigned NumElements) {
  dim3 Block(32, 2, 1);
  dim3 Grid((NumElements + 2 - 1) / 2, 1, 1);
  kernel_sgemm_NT_NT_m20_40_n9_40_k33_nps_c9310c2<<<Grid,Block>>>(A, ExtraOffsetA, B, ExtraOffsetB, C, ExtraOffsetC, NumElements);
  CHECK_ERR;
}
__global__ void
    __launch_bounds__(64)
kernel_sgemm_NT_NT_m20_24_n9_9_k9_spp_88ab2d4(const float * A, int ExtraOffsetA, const float ** B, int ExtraOffsetB, float ** C, int ExtraOffsetC, unsigned NumElements) {
  if ((threadIdx.y + blockDim.y * blockIdx.x) < NumElements) {
    const float* GlobMatA = &A[(threadIdx.y + blockDim.y * blockIdx.x) * 216 + 0 + ExtraOffsetA];
    const float* GlobMatB = &B[(threadIdx.y + blockDim.y * blockIdx.x)][0 + ExtraOffsetB];
    float* GlobMatC = &C[(threadIdx.y + blockDim.y * blockIdx.x)][0 + ExtraOffsetC];
    __shared__ float Scratch[162];
    float* ShrMatB = &Scratch[threadIdx.y * 81];

    // using ExtendedPatchLoader
    ShrMatB[threadIdx.x + 0] = GlobMatB[threadIdx.x + 0];
    ShrMatB[threadIdx.x + 32] = GlobMatB[threadIdx.x + 32];
    if (threadIdx.x < 17) {
      ShrMatB[threadIdx.x + 64] = GlobMatB[threadIdx.x + 64];
    }
    __syncthreads();

    if (threadIdx.x < 20) {

      float Results[9] = {0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f};
      float Value;

      for (int k = 0; k < 9; ++k) {
        Value = GlobMatA[threadIdx.x + 24 * k];

#pragma unroll
        for (int n = 0; n < 9; ++n) {
          Results[n] += Value * ShrMatB[k + 9 * n];
        }
      }

#pragma unroll
      for (int n = 0; n < 9; ++n) {
        GlobMatC[threadIdx.x + 24 * n] = Results[n] + GlobMatC[threadIdx.x + 24 * n];
      }
    }
  }
}
void sgemm_NT_NT_m20_24_n9_9_k9_spp_88ab2d4(const float * A, int ExtraOffsetA, const float ** B, int ExtraOffsetB, float ** C, int ExtraOffsetC, unsigned NumElements) {
  dim3 Block(32, 2, 1);
  dim3 Grid((NumElements + 2 - 1) / 2, 1, 1);
  kernel_sgemm_NT_NT_m20_24_n9_9_k9_spp_88ab2d4<<<Grid,Block>>>(A, ExtraOffsetA, B, ExtraOffsetB, C, ExtraOffsetC, NumElements);
  CHECK_ERR;
}
__global__ void
    __launch_bounds__(64)
kernel_sgemm_NT_NT_m20_40_n9_40_k34_nps_4344487(const float * __restrict__ A, int ExtraOffsetA, const float ** B, int ExtraOffsetB, float * C, int ExtraOffsetC, unsigned NumElements) {
  if ((threadIdx.y + blockDim.y * blockIdx.x) < NumElements) {
    const float* GlobMatA = &A[40 + ExtraOffsetA];
    const float* GlobMatB = &B[(threadIdx.y + blockDim.y * blockIdx.x)][1 + ExtraOffsetB];
    float* GlobMatC = &C[(threadIdx.y + blockDim.y * blockIdx.x) * 216 + 0 + ExtraOffsetC];
    __shared__ float Scratch[708];
    float* ShrMatB = &Scratch[threadIdx.y * 354];

    // using ExtendedPatchLoader
#pragma unroll
    for (int i = 0; i < 11; ++i) {
      ShrMatB[threadIdx.x + i * 32] = GlobMatB[threadIdx.x + i * 32];
    }
    if (threadIdx.x < 2) {
      ShrMatB[threadIdx.x + 352] = GlobMatB[threadIdx.x + 352];
    }
    __syncthreads();

    if (threadIdx.x < 20) {

      float Results[9] = {0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f};
      float Value;

      for (int k = 0; k < 34; ++k) {
        Value = GlobMatA[threadIdx.x + 40 * k];

#pragma unroll
        for (int n = 0; n < 9; ++n) {
          Results[n] += Value * ShrMatB[k + 40 * n];
        }
      }

#pragma unroll
      for (int n = 0; n < 9; ++n) {
        GlobMatC[threadIdx.x + 24 * n] = Results[n];
      }
    }
  }
}
void sgemm_NT_NT_m20_40_n9_40_k34_nps_4344487(const float * __restrict__ A, int ExtraOffsetA, const float ** B, int ExtraOffsetB, float * C, int ExtraOffsetC, unsigned NumElements) {
  dim3 Block(32, 2, 1);
  dim3 Grid((NumElements + 2 - 1) / 2, 1, 1);
  kernel_sgemm_NT_NT_m20_40_n9_40_k34_nps_4344487<<<Grid,Block>>>(A, ExtraOffsetA, B, ExtraOffsetB, C, ExtraOffsetC, NumElements);
  CHECK_ERR;
}
__global__ void
    __launch_bounds__(96)
kernel_sgemm_NT_NT_m10_40_n9_24_k17_nps_41e96ee(const float * __restrict__ A, int ExtraOffsetA, const float ** B, int ExtraOffsetB, float * C, int ExtraOffsetC, unsigned NumElements) {
  if ((threadIdx.y + blockDim.y * blockIdx.x) < NumElements) {
    const float* GlobMatA = &A[40 + ExtraOffsetA];
    const float* GlobMatB = &B[(threadIdx.y + blockDim.y * blockIdx.x)][1 + ExtraOffsetB];
    float* GlobMatC = &C[(threadIdx.y + blockDim.y * blockIdx.x) * 144 + 0 + ExtraOffsetC];
    __shared__ float Scratch[627];
    float* ShrMatB = &Scratch[threadIdx.y * 209];

    // using ExtendedPatchLoader
#pragma unroll
    for (int i = 0; i < 6; ++i) {
      ShrMatB[threadIdx.x + i * 32] = GlobMatB[threadIdx.x + i * 32];
    }
    if (threadIdx.x < 17) {
      ShrMatB[threadIdx.x + 192] = GlobMatB[threadIdx.x + 192];
    }
    __syncthreads();

    if (threadIdx.x < 10) {

      float Results[9] = {0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f};
      float Value;

      for (int k = 0; k < 17; ++k) {
        Value = GlobMatA[threadIdx.x + 40 * k];

#pragma unroll
        for (int n = 0; n < 9; ++n) {
          Results[n] += Value * ShrMatB[k + 24 * n];
        }
      }

#pragma unroll
      for (int n = 0; n < 9; ++n) {
        GlobMatC[threadIdx.x + 16 * n] = Results[n];
      }
    }
  }
}
void sgemm_NT_NT_m10_40_n9_24_k17_nps_41e96ee(const float * __restrict__ A, int ExtraOffsetA, const float ** B, int ExtraOffsetB, float * C, int ExtraOffsetC, unsigned NumElements) {
  dim3 Block(32, 3, 1);
  dim3 Grid((NumElements + 3 - 1) / 3, 1, 1);
  kernel_sgemm_NT_NT_m10_40_n9_24_k17_nps_41e96ee<<<Grid,Block>>>(A, ExtraOffsetA, B, ExtraOffsetB, C, ExtraOffsetC, NumElements);
  CHECK_ERR;
}
__global__ void
    __launch_bounds__(96)
kernel_sgemm_NT_NT_m10_16_n9_9_k9_spp_afbfa38(const float * A, int ExtraOffsetA, const float ** B, int ExtraOffsetB, float ** C, int ExtraOffsetC, unsigned NumElements) {
  if ((threadIdx.y + blockDim.y * blockIdx.x) < NumElements) {
    const float* GlobMatA = &A[(threadIdx.y + blockDim.y * blockIdx.x) * 144 + 0 + ExtraOffsetA];
    const float* GlobMatB = &B[(threadIdx.y + blockDim.y * blockIdx.x)][0 + ExtraOffsetB];
    float* GlobMatC = &C[(threadIdx.y + blockDim.y * blockIdx.x)][0 + ExtraOffsetC];
    __shared__ float Scratch[243];
    float* ShrMatB = &Scratch[threadIdx.y * 81];

    // using ExtendedPatchLoader
    ShrMatB[threadIdx.x + 0] = GlobMatB[threadIdx.x + 0];
    ShrMatB[threadIdx.x + 32] = GlobMatB[threadIdx.x + 32];
    if (threadIdx.x < 17) {
      ShrMatB[threadIdx.x + 64] = GlobMatB[threadIdx.x + 64];
    }
    __syncthreads();

    if (threadIdx.x < 10) {

      float Results[9] = {0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f};
      float Value;

      for (int k = 0; k < 9; ++k) {
        Value = GlobMatA[threadIdx.x + 16 * k];

#pragma unroll
        for (int n = 0; n < 9; ++n) {
          Results[n] += Value * ShrMatB[k + 9 * n];
        }
      }

#pragma unroll
      for (int n = 0; n < 9; ++n) {
        GlobMatC[threadIdx.x + 16 * n] = Results[n];
      }
    }
  }
}
void sgemm_NT_NT_m10_16_n9_9_k9_spp_afbfa38(const float * A, int ExtraOffsetA, const float ** B, int ExtraOffsetB, float ** C, int ExtraOffsetC, unsigned NumElements) {
  dim3 Block(32, 3, 1);
  dim3 Grid((NumElements + 3 - 1) / 3, 1, 1);
  kernel_sgemm_NT_NT_m10_16_n9_9_k9_spp_afbfa38<<<Grid,Block>>>(A, ExtraOffsetA, B, ExtraOffsetB, C, ExtraOffsetC, NumElements);
  CHECK_ERR;
}
__global__ void
    __launch_bounds__(96)
kernel_sgemm_NT_NT_m10_40_n9_24_k18_nps_20ec594(const float * __restrict__ A, int ExtraOffsetA, const float ** B, int ExtraOffsetB, float * C, int ExtraOffsetC, unsigned NumElements) {
  if ((threadIdx.y + blockDim.y * blockIdx.x) < NumElements) {
    const float* GlobMatA = &A[40 + ExtraOffsetA];
    const float* GlobMatB = &B[(threadIdx.y + blockDim.y * blockIdx.x)][1 + ExtraOffsetB];
    float* GlobMatC = &C[(threadIdx.y + blockDim.y * blockIdx.x) * 144 + 0 + ExtraOffsetC];
    __shared__ float Scratch[630];
    float* ShrMatB = &Scratch[threadIdx.y * 210];

    // using ExtendedPatchLoader
#pragma unroll
    for (int i = 0; i < 6; ++i) {
      ShrMatB[threadIdx.x + i * 32] = GlobMatB[threadIdx.x + i * 32];
    }
    if (threadIdx.x < 18) {
      ShrMatB[threadIdx.x + 192] = GlobMatB[threadIdx.x + 192];
    }
    __syncthreads();

    if (threadIdx.x < 10) {

      float Results[9] = {0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f};
      float Value;

      for (int k = 0; k < 18; ++k) {
        Value = GlobMatA[threadIdx.x + 40 * k];

#pragma unroll
        for (int n = 0; n < 9; ++n) {
          Results[n] += Value * ShrMatB[k + 24 * n];
        }
      }

#pragma unroll
      for (int n = 0; n < 9; ++n) {
        GlobMatC[threadIdx.x + 16 * n] = Results[n];
      }
    }
  }
}
void sgemm_NT_NT_m10_40_n9_24_k18_nps_20ec594(const float * __restrict__ A, int ExtraOffsetA, const float ** B, int ExtraOffsetB, float * C, int ExtraOffsetC, unsigned NumElements) {
  dim3 Block(32, 3, 1);
  dim3 Grid((NumElements + 3 - 1) / 3, 1, 1);
  kernel_sgemm_NT_NT_m10_40_n9_24_k18_nps_20ec594<<<Grid,Block>>>(A, ExtraOffsetA, B, ExtraOffsetB, C, ExtraOffsetC, NumElements);
  CHECK_ERR;
}
__global__ void
    __launch_bounds__(96)
kernel_sgemm_NT_NT_m10_16_n9_9_k9_spp_fe81cdb(const float * A, int ExtraOffsetA, const float ** B, int ExtraOffsetB, float ** C, int ExtraOffsetC, unsigned NumElements) {
  if ((threadIdx.y + blockDim.y * blockIdx.x) < NumElements) {
    const float* GlobMatA = &A[(threadIdx.y + blockDim.y * blockIdx.x) * 144 + 0 + ExtraOffsetA];
    const float* GlobMatB = &B[(threadIdx.y + blockDim.y * blockIdx.x)][0 + ExtraOffsetB];
    float* GlobMatC = &C[(threadIdx.y + blockDim.y * blockIdx.x)][0 + ExtraOffsetC];
    __shared__ float Scratch[243];
    float* ShrMatB = &Scratch[threadIdx.y * 81];

    // using ExtendedPatchLoader
    ShrMatB[threadIdx.x + 0] = GlobMatB[threadIdx.x + 0];
    ShrMatB[threadIdx.x + 32] = GlobMatB[threadIdx.x + 32];
    if (threadIdx.x < 17) {
      ShrMatB[threadIdx.x + 64] = GlobMatB[threadIdx.x + 64];
    }
    __syncthreads();

    if (threadIdx.x < 10) {

      float Results[9] = {0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f};
      float Value;

      for (int k = 0; k < 9; ++k) {
        Value = GlobMatA[threadIdx.x + 16 * k];

#pragma unroll
        for (int n = 0; n < 9; ++n) {
          Results[n] += Value * ShrMatB[k + 9 * n];
        }
      }

#pragma unroll
      for (int n = 0; n < 9; ++n) {
        GlobMatC[threadIdx.x + 16 * n] = Results[n] + GlobMatC[threadIdx.x + 16 * n];
      }
    }
  }
}
void sgemm_NT_NT_m10_16_n9_9_k9_spp_fe81cdb(const float * A, int ExtraOffsetA, const float ** B, int ExtraOffsetB, float ** C, int ExtraOffsetC, unsigned NumElements) {
  dim3 Block(32, 3, 1);
  dim3 Grid((NumElements + 3 - 1) / 3, 1, 1);
  kernel_sgemm_NT_NT_m10_16_n9_9_k9_spp_fe81cdb<<<Grid,Block>>>(A, ExtraOffsetA, B, ExtraOffsetB, C, ExtraOffsetC, NumElements);
  CHECK_ERR;
}
__global__ void
    __launch_bounds__(96)
kernel_sgemm_NT_NT_m10_40_n9_24_k19_nps_277801b(const float * __restrict__ A, int ExtraOffsetA, const float ** B, int ExtraOffsetB, float * C, int ExtraOffsetC, unsigned NumElements) {
  if ((threadIdx.y + blockDim.y * blockIdx.x) < NumElements) {
    const float* GlobMatA = &A[40 + ExtraOffsetA];
    const float* GlobMatB = &B[(threadIdx.y + blockDim.y * blockIdx.x)][1 + ExtraOffsetB];
    float* GlobMatC = &C[(threadIdx.y + blockDim.y * blockIdx.x) * 144 + 0 + ExtraOffsetC];
    __shared__ float Scratch[633];
    float* ShrMatB = &Scratch[threadIdx.y * 211];

    // using ExtendedPatchLoader
#pragma unroll
    for (int i = 0; i < 6; ++i) {
      ShrMatB[threadIdx.x + i * 32] = GlobMatB[threadIdx.x + i * 32];
    }
    if (threadIdx.x < 19) {
      ShrMatB[threadIdx.x + 192] = GlobMatB[threadIdx.x + 192];
    }
    __syncthreads();

    if (threadIdx.x < 10) {

      float Results[9] = {0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f};
      float Value;

      for (int k = 0; k < 19; ++k) {
        Value = GlobMatA[threadIdx.x + 40 * k];

#pragma unroll
        for (int n = 0; n < 9; ++n) {
          Results[n] += Value * ShrMatB[k + 24 * n];
        }
      }

#pragma unroll
      for (int n = 0; n < 9; ++n) {
        GlobMatC[threadIdx.x + 16 * n] = Results[n];
      }
    }
  }
}
void sgemm_NT_NT_m10_40_n9_24_k19_nps_277801b(const float * __restrict__ A, int ExtraOffsetA, const float ** B, int ExtraOffsetB, float * C, int ExtraOffsetC, unsigned NumElements) {
  dim3 Block(32, 3, 1);
  dim3 Grid((NumElements + 3 - 1) / 3, 1, 1);
  kernel_sgemm_NT_NT_m10_40_n9_24_k19_nps_277801b<<<Grid,Block>>>(A, ExtraOffsetA, B, ExtraOffsetB, C, ExtraOffsetC, NumElements);
  CHECK_ERR;
}
__global__ void
    __launch_bounds__(96)
kernel_sgemm_NT_NT_m4_40_n9_16_k7_nps_625feb0(const float * __restrict__ A, int ExtraOffsetA, const float ** B, int ExtraOffsetB, float * C, int ExtraOffsetC, unsigned NumElements) {
  if ((threadIdx.y + blockDim.y * blockIdx.x) < NumElements) {
    const float* GlobMatA = &A[40 + ExtraOffsetA];
    const float* GlobMatB = &B[(threadIdx.y + blockDim.y * blockIdx.x)][1 + ExtraOffsetB];
    float* GlobMatC = &C[(threadIdx.y + blockDim.y * blockIdx.x) * 72 + 0 + ExtraOffsetC];
    __shared__ float Scratch[405];
    float* ShrMatB = &Scratch[threadIdx.y * 135];

    // using ExtendedPatchLoader
#pragma unroll
    for (int i = 0; i < 4; ++i) {
      ShrMatB[threadIdx.x + i * 32] = GlobMatB[threadIdx.x + i * 32];
    }
    if (threadIdx.x < 7) {
      ShrMatB[threadIdx.x + 128] = GlobMatB[threadIdx.x + 128];
    }
    __syncthreads();

    if (threadIdx.x < 4) {

      float Results[9] = {0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f};
      float Value;

      for (int k = 0; k < 7; ++k) {
        Value = GlobMatA[threadIdx.x + 40 * k];

#pragma unroll
        for (int n = 0; n < 9; ++n) {
          Results[n] += Value * ShrMatB[k + 16 * n];
        }
      }

#pragma unroll
      for (int n = 0; n < 9; ++n) {
        GlobMatC[threadIdx.x + 8 * n] = Results[n];
      }
    }
  }
}
void sgemm_NT_NT_m4_40_n9_16_k7_nps_625feb0(const float * __restrict__ A, int ExtraOffsetA, const float ** B, int ExtraOffsetB, float * C, int ExtraOffsetC, unsigned NumElements) {
  dim3 Block(32, 3, 1);
  dim3 Grid((NumElements + 3 - 1) / 3, 1, 1);
  kernel_sgemm_NT_NT_m4_40_n9_16_k7_nps_625feb0<<<Grid,Block>>>(A, ExtraOffsetA, B, ExtraOffsetB, C, ExtraOffsetC, NumElements);
  CHECK_ERR;
}
__global__ void
    __launch_bounds__(96)
kernel_sgemm_NT_NT_m4_8_n9_9_k9_spp_67deb58(const float * A, int ExtraOffsetA, const float ** B, int ExtraOffsetB, float ** C, int ExtraOffsetC, unsigned NumElements) {
  if ((threadIdx.y + blockDim.y * blockIdx.x) < NumElements) {
    const float* GlobMatA = &A[(threadIdx.y + blockDim.y * blockIdx.x) * 72 + 0 + ExtraOffsetA];
    const float* GlobMatB = &B[(threadIdx.y + blockDim.y * blockIdx.x)][0 + ExtraOffsetB];
    float* GlobMatC = &C[(threadIdx.y + blockDim.y * blockIdx.x)][0 + ExtraOffsetC];
    __shared__ float Scratch[243];
    float* ShrMatB = &Scratch[threadIdx.y * 81];

    // using ExtendedPatchLoader
    ShrMatB[threadIdx.x + 0] = GlobMatB[threadIdx.x + 0];
    ShrMatB[threadIdx.x + 32] = GlobMatB[threadIdx.x + 32];
    if (threadIdx.x < 17) {
      ShrMatB[threadIdx.x + 64] = GlobMatB[threadIdx.x + 64];
    }
    __syncthreads();

    if (threadIdx.x < 4) {

      float Results[9] = {0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f};
      float Value;

      for (int k = 0; k < 9; ++k) {
        Value = GlobMatA[threadIdx.x + 8 * k];

#pragma unroll
        for (int n = 0; n < 9; ++n) {
          Results[n] += Value * ShrMatB[k + 9 * n];
        }
      }

#pragma unroll
      for (int n = 0; n < 9; ++n) {
        GlobMatC[threadIdx.x + 8 * n] = Results[n];
      }
    }
  }
}
void sgemm_NT_NT_m4_8_n9_9_k9_spp_67deb58(const float * A, int ExtraOffsetA, const float ** B, int ExtraOffsetB, float ** C, int ExtraOffsetC, unsigned NumElements) {
  dim3 Block(32, 3, 1);
  dim3 Grid((NumElements + 3 - 1) / 3, 1, 1);
  kernel_sgemm_NT_NT_m4_8_n9_9_k9_spp_67deb58<<<Grid,Block>>>(A, ExtraOffsetA, B, ExtraOffsetB, C, ExtraOffsetC, NumElements);
  CHECK_ERR;
}
__global__ void
    __launch_bounds__(96)
kernel_sgemm_NT_NT_m4_40_n9_16_k8_nps_8f843c2(const float * __restrict__ A, int ExtraOffsetA, const float ** B, int ExtraOffsetB, float * C, int ExtraOffsetC, unsigned NumElements) {
  if ((threadIdx.y + blockDim.y * blockIdx.x) < NumElements) {
    const float* GlobMatA = &A[40 + ExtraOffsetA];
    const float* GlobMatB = &B[(threadIdx.y + blockDim.y * blockIdx.x)][1 + ExtraOffsetB];
    float* GlobMatC = &C[(threadIdx.y + blockDim.y * blockIdx.x) * 72 + 0 + ExtraOffsetC];
    __shared__ float Scratch[408];
    float* ShrMatB = &Scratch[threadIdx.y * 136];

    // using ExtendedPatchLoader
#pragma unroll
    for (int i = 0; i < 4; ++i) {
      ShrMatB[threadIdx.x + i * 32] = GlobMatB[threadIdx.x + i * 32];
    }
    if (threadIdx.x < 8) {
      ShrMatB[threadIdx.x + 128] = GlobMatB[threadIdx.x + 128];
    }
    __syncthreads();

    if (threadIdx.x < 4) {

      float Results[9] = {0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f};
      float Value;

      for (int k = 0; k < 8; ++k) {
        Value = GlobMatA[threadIdx.x + 40 * k];

#pragma unroll
        for (int n = 0; n < 9; ++n) {
          Results[n] += Value * ShrMatB[k + 16 * n];
        }
      }

#pragma unroll
      for (int n = 0; n < 9; ++n) {
        GlobMatC[threadIdx.x + 8 * n] = Results[n];
      }
    }
  }
}
void sgemm_NT_NT_m4_40_n9_16_k8_nps_8f843c2(const float * __restrict__ A, int ExtraOffsetA, const float ** B, int ExtraOffsetB, float * C, int ExtraOffsetC, unsigned NumElements) {
  dim3 Block(32, 3, 1);
  dim3 Grid((NumElements + 3 - 1) / 3, 1, 1);
  kernel_sgemm_NT_NT_m4_40_n9_16_k8_nps_8f843c2<<<Grid,Block>>>(A, ExtraOffsetA, B, ExtraOffsetB, C, ExtraOffsetC, NumElements);
  CHECK_ERR;
}
__global__ void
    __launch_bounds__(96)
kernel_sgemm_NT_NT_m4_8_n9_9_k9_spp_784eae3(const float * A, int ExtraOffsetA, const float ** B, int ExtraOffsetB, float ** C, int ExtraOffsetC, unsigned NumElements) {
  if ((threadIdx.y + blockDim.y * blockIdx.x) < NumElements) {
    const float* GlobMatA = &A[(threadIdx.y + blockDim.y * blockIdx.x) * 72 + 0 + ExtraOffsetA];
    const float* GlobMatB = &B[(threadIdx.y + blockDim.y * blockIdx.x)][0 + ExtraOffsetB];
    float* GlobMatC = &C[(threadIdx.y + blockDim.y * blockIdx.x)][0 + ExtraOffsetC];
    __shared__ float Scratch[243];
    float* ShrMatB = &Scratch[threadIdx.y * 81];

    // using ExtendedPatchLoader
    ShrMatB[threadIdx.x + 0] = GlobMatB[threadIdx.x + 0];
    ShrMatB[threadIdx.x + 32] = GlobMatB[threadIdx.x + 32];
    if (threadIdx.x < 17) {
      ShrMatB[threadIdx.x + 64] = GlobMatB[threadIdx.x + 64];
    }
    __syncthreads();

    if (threadIdx.x < 4) {

      float Results[9] = {0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f};
      float Value;

      for (int k = 0; k < 9; ++k) {
        Value = GlobMatA[threadIdx.x + 8 * k];

#pragma unroll
        for (int n = 0; n < 9; ++n) {
          Results[n] += Value * ShrMatB[k + 9 * n];
        }
      }

#pragma unroll
      for (int n = 0; n < 9; ++n) {
        GlobMatC[threadIdx.x + 8 * n] = Results[n] + GlobMatC[threadIdx.x + 8 * n];
      }
    }
  }
}
void sgemm_NT_NT_m4_8_n9_9_k9_spp_784eae3(const float * A, int ExtraOffsetA, const float ** B, int ExtraOffsetB, float ** C, int ExtraOffsetC, unsigned NumElements) {
  dim3 Block(32, 3, 1);
  dim3 Grid((NumElements + 3 - 1) / 3, 1, 1);
  kernel_sgemm_NT_NT_m4_8_n9_9_k9_spp_784eae3<<<Grid,Block>>>(A, ExtraOffsetA, B, ExtraOffsetB, C, ExtraOffsetC, NumElements);
  CHECK_ERR;
}
__global__ void
    __launch_bounds__(96)
kernel_sgemm_NT_NT_m4_40_n9_16_k9_nps_d74565c(const float * __restrict__ A, int ExtraOffsetA, const float ** B, int ExtraOffsetB, float * C, int ExtraOffsetC, unsigned NumElements) {
  if ((threadIdx.y + blockDim.y * blockIdx.x) < NumElements) {
    const float* GlobMatA = &A[40 + ExtraOffsetA];
    const float* GlobMatB = &B[(threadIdx.y + blockDim.y * blockIdx.x)][1 + ExtraOffsetB];
    float* GlobMatC = &C[(threadIdx.y + blockDim.y * blockIdx.x) * 72 + 0 + ExtraOffsetC];
    __shared__ float Scratch[411];
    float* ShrMatB = &Scratch[threadIdx.y * 137];

    // using ExtendedPatchLoader
#pragma unroll
    for (int i = 0; i < 4; ++i) {
      ShrMatB[threadIdx.x + i * 32] = GlobMatB[threadIdx.x + i * 32];
    }
    if (threadIdx.x < 9) {
      ShrMatB[threadIdx.x + 128] = GlobMatB[threadIdx.x + 128];
    }
    __syncthreads();

    if (threadIdx.x < 4) {

      float Results[9] = {0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f};
      float Value;

      for (int k = 0; k < 9; ++k) {
        Value = GlobMatA[threadIdx.x + 40 * k];

#pragma unroll
        for (int n = 0; n < 9; ++n) {
          Results[n] += Value * ShrMatB[k + 16 * n];
        }
      }

#pragma unroll
      for (int n = 0; n < 9; ++n) {
        GlobMatC[threadIdx.x + 8 * n] = Results[n];
      }
    }
  }
}
void sgemm_NT_NT_m4_40_n9_16_k9_nps_d74565c(const float * __restrict__ A, int ExtraOffsetA, const float ** B, int ExtraOffsetB, float * C, int ExtraOffsetC, unsigned NumElements) {
  dim3 Block(32, 3, 1);
  dim3 Grid((NumElements + 3 - 1) / 3, 1, 1);
  kernel_sgemm_NT_NT_m4_40_n9_16_k9_nps_d74565c<<<Grid,Block>>>(A, ExtraOffsetA, B, ExtraOffsetB, C, ExtraOffsetC, NumElements);
  CHECK_ERR;
}
__global__ void
    __launch_bounds__(96)
kernel_sgemm_NT_NT_m1_40_n9_8_k1_nps_2533254(const float * __restrict__ A, int ExtraOffsetA, const float ** B, int ExtraOffsetB, float * C, int ExtraOffsetC, unsigned NumElements) {
  if ((threadIdx.y + blockDim.y * blockIdx.x) < NumElements) {
    const float* GlobMatA = &A[40 + ExtraOffsetA];
    const float* GlobMatB = &B[(threadIdx.y + blockDim.y * blockIdx.x)][1 + ExtraOffsetB];
    float* GlobMatC = &C[(threadIdx.y + blockDim.y * blockIdx.x) * 72 + 0 + ExtraOffsetC];
    __shared__ float Scratch[195];
    float* ShrMatB = &Scratch[threadIdx.y * 65];

    // using ExtendedPatchLoader
    ShrMatB[threadIdx.x + 0] = GlobMatB[threadIdx.x + 0];
    ShrMatB[threadIdx.x + 32] = GlobMatB[threadIdx.x + 32];
    if (threadIdx.x < 1) {
      ShrMatB[threadIdx.x + 64] = GlobMatB[threadIdx.x + 64];
    }
    __syncthreads();

    if (threadIdx.x < 1) {

      float Results[9] = {0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f};
      float Value;

      for (int k = 0; k < 1; ++k) {
        Value = GlobMatA[threadIdx.x + 40 * k];

#pragma unroll
        for (int n = 0; n < 9; ++n) {
          Results[n] += Value * ShrMatB[k + 8 * n];
        }
      }

#pragma unroll
      for (int n = 0; n < 9; ++n) {
        GlobMatC[threadIdx.x + 8 * n] = Results[n];
      }
    }
  }
}
void sgemm_NT_NT_m1_40_n9_8_k1_nps_2533254(const float * __restrict__ A, int ExtraOffsetA, const float ** B, int ExtraOffsetB, float * C, int ExtraOffsetC, unsigned NumElements) {
  dim3 Block(32, 3, 1);
  dim3 Grid((NumElements + 3 - 1) / 3, 1, 1);
  kernel_sgemm_NT_NT_m1_40_n9_8_k1_nps_2533254<<<Grid,Block>>>(A, ExtraOffsetA, B, ExtraOffsetB, C, ExtraOffsetC, NumElements);
  CHECK_ERR;
}
__global__ void
    __launch_bounds__(96)
kernel_sgemm_NT_NT_m1_8_n9_9_k9_spp_5c45fb8(const float * A, int ExtraOffsetA, const float ** B, int ExtraOffsetB, float ** C, int ExtraOffsetC, unsigned NumElements) {
  if ((threadIdx.y + blockDim.y * blockIdx.x) < NumElements) {
    const float* GlobMatA = &A[(threadIdx.y + blockDim.y * blockIdx.x) * 72 + 0 + ExtraOffsetA];
    const float* GlobMatB = &B[(threadIdx.y + blockDim.y * blockIdx.x)][0 + ExtraOffsetB];
    float* GlobMatC = &C[(threadIdx.y + blockDim.y * blockIdx.x)][0 + ExtraOffsetC];
    __shared__ float Scratch[243];
    float* ShrMatB = &Scratch[threadIdx.y * 81];

    // using ExtendedPatchLoader
    ShrMatB[threadIdx.x + 0] = GlobMatB[threadIdx.x + 0];
    ShrMatB[threadIdx.x + 32] = GlobMatB[threadIdx.x + 32];
    if (threadIdx.x < 17) {
      ShrMatB[threadIdx.x + 64] = GlobMatB[threadIdx.x + 64];
    }
    __syncthreads();

    if (threadIdx.x < 1) {

      float Results[9] = {0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f};
      float Value;

      for (int k = 0; k < 9; ++k) {
        Value = GlobMatA[threadIdx.x + 8 * k];

#pragma unroll
        for (int n = 0; n < 9; ++n) {
          Results[n] += Value * ShrMatB[k + 9 * n];
        }
      }

#pragma unroll
      for (int n = 0; n < 9; ++n) {
        GlobMatC[threadIdx.x + 8 * n] = Results[n];
      }
    }
  }
}
void sgemm_NT_NT_m1_8_n9_9_k9_spp_5c45fb8(const float * A, int ExtraOffsetA, const float ** B, int ExtraOffsetB, float ** C, int ExtraOffsetC, unsigned NumElements) {
  dim3 Block(32, 3, 1);
  dim3 Grid((NumElements + 3 - 1) / 3, 1, 1);
  kernel_sgemm_NT_NT_m1_8_n9_9_k9_spp_5c45fb8<<<Grid,Block>>>(A, ExtraOffsetA, B, ExtraOffsetB, C, ExtraOffsetC, NumElements);
  CHECK_ERR;
}
__global__ void
    __launch_bounds__(96)
kernel_sgemm_NT_NT_m1_40_n9_8_k2_nps_a20e3a8(const float * __restrict__ A, int ExtraOffsetA, const float ** B, int ExtraOffsetB, float * C, int ExtraOffsetC, unsigned NumElements) {
  if ((threadIdx.y + blockDim.y * blockIdx.x) < NumElements) {
    const float* GlobMatA = &A[40 + ExtraOffsetA];
    const float* GlobMatB = &B[(threadIdx.y + blockDim.y * blockIdx.x)][1 + ExtraOffsetB];
    float* GlobMatC = &C[(threadIdx.y + blockDim.y * blockIdx.x) * 72 + 0 + ExtraOffsetC];
    __shared__ float Scratch[198];
    float* ShrMatB = &Scratch[threadIdx.y * 66];

    // using ExtendedPatchLoader
    ShrMatB[threadIdx.x + 0] = GlobMatB[threadIdx.x + 0];
    ShrMatB[threadIdx.x + 32] = GlobMatB[threadIdx.x + 32];
    if (threadIdx.x < 2) {
      ShrMatB[threadIdx.x + 64] = GlobMatB[threadIdx.x + 64];
    }
    __syncthreads();

    if (threadIdx.x < 1) {

      float Results[9] = {0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f};
      float Value;

      for (int k = 0; k < 2; ++k) {
        Value = GlobMatA[threadIdx.x + 40 * k];

#pragma unroll
        for (int n = 0; n < 9; ++n) {
          Results[n] += Value * ShrMatB[k + 8 * n];
        }
      }

#pragma unroll
      for (int n = 0; n < 9; ++n) {
        GlobMatC[threadIdx.x + 8 * n] = Results[n];
      }
    }
  }
}
void sgemm_NT_NT_m1_40_n9_8_k2_nps_a20e3a8(const float * __restrict__ A, int ExtraOffsetA, const float ** B, int ExtraOffsetB, float * C, int ExtraOffsetC, unsigned NumElements) {
  dim3 Block(32, 3, 1);
  dim3 Grid((NumElements + 3 - 1) / 3, 1, 1);
  kernel_sgemm_NT_NT_m1_40_n9_8_k2_nps_a20e3a8<<<Grid,Block>>>(A, ExtraOffsetA, B, ExtraOffsetB, C, ExtraOffsetC, NumElements);
  CHECK_ERR;
}
__global__ void
    __launch_bounds__(96)
kernel_sgemm_NT_NT_m1_8_n9_9_k9_spp_7cd43db(const float * A, int ExtraOffsetA, const float ** B, int ExtraOffsetB, float ** C, int ExtraOffsetC, unsigned NumElements) {
  if ((threadIdx.y + blockDim.y * blockIdx.x) < NumElements) {
    const float* GlobMatA = &A[(threadIdx.y + blockDim.y * blockIdx.x) * 72 + 0 + ExtraOffsetA];
    const float* GlobMatB = &B[(threadIdx.y + blockDim.y * blockIdx.x)][0 + ExtraOffsetB];
    float* GlobMatC = &C[(threadIdx.y + blockDim.y * blockIdx.x)][0 + ExtraOffsetC];
    __shared__ float Scratch[243];
    float* ShrMatB = &Scratch[threadIdx.y * 81];

    // using ExtendedPatchLoader
    ShrMatB[threadIdx.x + 0] = GlobMatB[threadIdx.x + 0];
    ShrMatB[threadIdx.x + 32] = GlobMatB[threadIdx.x + 32];
    if (threadIdx.x < 17) {
      ShrMatB[threadIdx.x + 64] = GlobMatB[threadIdx.x + 64];
    }
    __syncthreads();

    if (threadIdx.x < 1) {

      float Results[9] = {0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f};
      float Value;

      for (int k = 0; k < 9; ++k) {
        Value = GlobMatA[threadIdx.x + 8 * k];

#pragma unroll
        for (int n = 0; n < 9; ++n) {
          Results[n] += Value * ShrMatB[k + 9 * n];
        }
      }

#pragma unroll
      for (int n = 0; n < 9; ++n) {
        GlobMatC[threadIdx.x + 8 * n] = Results[n] + GlobMatC[threadIdx.x + 8 * n];
      }
    }
  }
}
void sgemm_NT_NT_m1_8_n9_9_k9_spp_7cd43db(const float * A, int ExtraOffsetA, const float ** B, int ExtraOffsetB, float ** C, int ExtraOffsetC, unsigned NumElements) {
  dim3 Block(32, 3, 1);
  dim3 Grid((NumElements + 3 - 1) / 3, 1, 1);
  kernel_sgemm_NT_NT_m1_8_n9_9_k9_spp_7cd43db<<<Grid,Block>>>(A, ExtraOffsetA, B, ExtraOffsetB, C, ExtraOffsetC, NumElements);
  CHECK_ERR;
}
__global__ void
    __launch_bounds__(96)
kernel_sgemm_NT_NT_m1_40_n9_8_k3_nps_cd5eaf9(const float * __restrict__ A, int ExtraOffsetA, const float ** B, int ExtraOffsetB, float * C, int ExtraOffsetC, unsigned NumElements) {
  if ((threadIdx.y + blockDim.y * blockIdx.x) < NumElements) {
    const float* GlobMatA = &A[40 + ExtraOffsetA];
    const float* GlobMatB = &B[(threadIdx.y + blockDim.y * blockIdx.x)][1 + ExtraOffsetB];
    float* GlobMatC = &C[(threadIdx.y + blockDim.y * blockIdx.x) * 72 + 0 + ExtraOffsetC];
    __shared__ float Scratch[201];
    float* ShrMatB = &Scratch[threadIdx.y * 67];

    // using ExtendedPatchLoader
    ShrMatB[threadIdx.x + 0] = GlobMatB[threadIdx.x + 0];
    ShrMatB[threadIdx.x + 32] = GlobMatB[threadIdx.x + 32];
    if (threadIdx.x < 3) {
      ShrMatB[threadIdx.x + 64] = GlobMatB[threadIdx.x + 64];
    }
    __syncthreads();

    if (threadIdx.x < 1) {

      float Results[9] = {0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f};
      float Value;

      for (int k = 0; k < 3; ++k) {
        Value = GlobMatA[threadIdx.x + 40 * k];

#pragma unroll
        for (int n = 0; n < 9; ++n) {
          Results[n] += Value * ShrMatB[k + 8 * n];
        }
      }

#pragma unroll
      for (int n = 0; n < 9; ++n) {
        GlobMatC[threadIdx.x + 8 * n] = Results[n];
      }
    }
  }
}
void sgemm_NT_NT_m1_40_n9_8_k3_nps_cd5eaf9(const float * __restrict__ A, int ExtraOffsetA, const float ** B, int ExtraOffsetB, float * C, int ExtraOffsetC, unsigned NumElements) {
  dim3 Block(32, 3, 1);
  dim3 Grid((NumElements + 3 - 1) / 3, 1, 1);
  kernel_sgemm_NT_NT_m1_40_n9_8_k3_nps_cd5eaf9<<<Grid,Block>>>(A, ExtraOffsetA, B, ExtraOffsetB, C, ExtraOffsetC, NumElements);
  CHECK_ERR;
}
__global__ void
    __launch_bounds__(64)
kernel_sgemm_NT_NT_m49_56_n9_56_k56_nps_4ae838d(const float * __restrict__ A, int ExtraOffsetA, const float ** B, int ExtraOffsetB, float * C, int ExtraOffsetC, unsigned NumElements) {
  if ((threadIdx.y + blockDim.y * blockIdx.x) < NumElements) {
    const float* GlobMatA = &A[0 + ExtraOffsetA];
    const float* GlobMatB = &B[(threadIdx.y + blockDim.y * blockIdx.x)][0 + ExtraOffsetB];
    float* GlobMatC = &C[(threadIdx.y + blockDim.y * blockIdx.x) * 504 + 0 + ExtraOffsetC];
    __shared__ float Scratch[504];
    float* ShrMatB = &Scratch[threadIdx.y * 504];

    // using ExtendedPatchLoader
#pragma unroll
    for (int i = 0; i < 7; ++i) {
      ShrMatB[threadIdx.x + i * 64] = GlobMatB[threadIdx.x + i * 64];
    }
    if (threadIdx.x < 56) {
      ShrMatB[threadIdx.x + 448] = GlobMatB[threadIdx.x + 448];
    }
    __syncthreads();

    if (threadIdx.x < 49) {

      float Results[9] = {0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f};
      float Value;

      for (int k = 0; k < 56; ++k) {
        Value = GlobMatA[threadIdx.x + 56 * k];

#pragma unroll
        for (int n = 0; n < 9; ++n) {
          Results[n] += Value * ShrMatB[k + 56 * n];
        }
      }

#pragma unroll
      for (int n = 0; n < 9; ++n) {
        GlobMatC[threadIdx.x + 56 * n] = Results[n];
      }
    }
  }
}
void sgemm_NT_NT_m49_56_n9_56_k56_nps_4ae838d(const float * __restrict__ A, int ExtraOffsetA, const float ** B, int ExtraOffsetB, float * C, int ExtraOffsetC, unsigned NumElements) {
  dim3 Block(64, 1, 1);
  dim3 Grid((NumElements + 1 - 1) / 1, 1, 1);
  kernel_sgemm_NT_NT_m49_56_n9_56_k56_nps_4ae838d<<<Grid,Block>>>(A, ExtraOffsetA, B, ExtraOffsetB, C, ExtraOffsetC, NumElements);
  CHECK_ERR;
}
__global__ void
    __launch_bounds__(64)
kernel_sgemm_NT_NT_m49_56_n9_9_k9_spp_539b719(const float * A, int ExtraOffsetA, const float ** B, int ExtraOffsetB, float ** C, int ExtraOffsetC, unsigned NumElements) {
  if ((threadIdx.y + blockDim.y * blockIdx.x) < NumElements) {
    const float* GlobMatA = &A[(threadIdx.y + blockDim.y * blockIdx.x) * 504 + 0 + ExtraOffsetA];
    const float* GlobMatB = &B[(threadIdx.y + blockDim.y * blockIdx.x)][0 + ExtraOffsetB];
    float* GlobMatC = &C[(threadIdx.y + blockDim.y * blockIdx.x)][0 + ExtraOffsetC];
    __shared__ float Scratch[81];
    float* ShrMatB = &Scratch[threadIdx.y * 81];

    // using ExtendedPatchLoader
    ShrMatB[threadIdx.x + 0] = GlobMatB[threadIdx.x + 0];
    if (threadIdx.x < 17) {
      ShrMatB[threadIdx.x + 64] = GlobMatB[threadIdx.x + 64];
    }
    __syncthreads();

    if (threadIdx.x < 49) {

      float Results[9] = {0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f};
      float Value;

      for (int k = 0; k < 9; ++k) {
        Value = GlobMatA[threadIdx.x + 56 * k];

#pragma unroll
        for (int n = 0; n < 9; ++n) {
          Results[n] += Value * ShrMatB[k + 9 * n];
        }
      }

#pragma unroll
      for (int n = 0; n < 9; ++n) {
        GlobMatC[threadIdx.x + 56 * n] = Results[n];
      }
    }
  }
}
void sgemm_NT_NT_m49_56_n9_9_k9_spp_539b719(const float * A, int ExtraOffsetA, const float ** B, int ExtraOffsetB, float ** C, int ExtraOffsetC, unsigned NumElements) {
  dim3 Block(64, 1, 1);
  dim3 Grid((NumElements + 1 - 1) / 1, 1, 1);
  kernel_sgemm_NT_NT_m49_56_n9_9_k9_spp_539b719<<<Grid,Block>>>(A, ExtraOffsetA, B, ExtraOffsetB, C, ExtraOffsetC, NumElements);
  CHECK_ERR;
}
__global__ void
    __launch_bounds__(64)
kernel_sgemm_NT_NT_m49_56_n9_9_k9_spp_082a417(const float * A, int ExtraOffsetA, const float ** B, int ExtraOffsetB, float ** C, int ExtraOffsetC, unsigned NumElements) {
  if ((threadIdx.y + blockDim.y * blockIdx.x) < NumElements) {
    const float* GlobMatA = &A[(threadIdx.y + blockDim.y * blockIdx.x) * 504 + 0 + ExtraOffsetA];
    const float* GlobMatB = &B[(threadIdx.y + blockDim.y * blockIdx.x)][0 + ExtraOffsetB];
    float* GlobMatC = &C[(threadIdx.y + blockDim.y * blockIdx.x)][0 + ExtraOffsetC];
    __shared__ float Scratch[81];
    float* ShrMatB = &Scratch[threadIdx.y * 81];

    // using ExtendedPatchLoader
    ShrMatB[threadIdx.x + 0] = GlobMatB[threadIdx.x + 0];
    if (threadIdx.x < 17) {
      ShrMatB[threadIdx.x + 64] = GlobMatB[threadIdx.x + 64];
    }
    __syncthreads();

    if (threadIdx.x < 49) {

      float Results[9] = {0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f};
      float Value;

      for (int k = 0; k < 9; ++k) {
        Value = GlobMatA[threadIdx.x + 56 * k];

#pragma unroll
        for (int n = 0; n < 9; ++n) {
          Results[n] += Value * ShrMatB[k + 9 * n];
        }
      }

#pragma unroll
      for (int n = 0; n < 9; ++n) {
        GlobMatC[threadIdx.x + 56 * n] = Results[n] + GlobMatC[threadIdx.x + 56 * n];
      }
    }
  }
}
void sgemm_NT_NT_m49_56_n9_9_k9_spp_082a417(const float * A, int ExtraOffsetA, const float ** B, int ExtraOffsetB, float ** C, int ExtraOffsetC, unsigned NumElements) {
  dim3 Block(64, 1, 1);
  dim3 Grid((NumElements + 1 - 1) / 1, 1, 1);
  kernel_sgemm_NT_NT_m49_56_n9_9_k9_spp_082a417<<<Grid,Block>>>(A, ExtraOffsetA, B, ExtraOffsetB, C, ExtraOffsetC, NumElements);
  CHECK_ERR;
}
__global__ void
    __launch_bounds__(64)
kernel_sgemm_NT_NT_m49_56_n9_9_k9_pps_834df01(const float ** A, int ExtraOffsetA, const float ** B, int ExtraOffsetB, float * C, int ExtraOffsetC, unsigned NumElements) {
  if ((threadIdx.y + blockDim.y * blockIdx.x) < NumElements) {
    const float* GlobMatA = &A[(threadIdx.y + blockDim.y * blockIdx.x)][0 + ExtraOffsetA];
    const float* GlobMatB = &B[(threadIdx.y + blockDim.y * blockIdx.x)][0 + ExtraOffsetB];
    float* GlobMatC = &C[(threadIdx.y + blockDim.y * blockIdx.x) * 504 + 0 + ExtraOffsetC];
    __shared__ float Scratch[81];
    float* ShrMatB = &Scratch[threadIdx.y * 81];

    // using ExtendedPatchLoader
    ShrMatB[threadIdx.x + 0] = GlobMatB[threadIdx.x + 0];
    if (threadIdx.x < 17) {
      ShrMatB[threadIdx.x + 64] = GlobMatB[threadIdx.x + 64];
    }
    __syncthreads();

    if (threadIdx.x < 49) {

      float Results[9] = {0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f};
      float Value;

      for (int k = 0; k < 9; ++k) {
        Value = GlobMatA[threadIdx.x + 56 * k];

#pragma unroll
        for (int n = 0; n < 9; ++n) {
          Results[n] += Value * ShrMatB[k + 9 * n];
        }
      }

#pragma unroll
      for (int n = 0; n < 9; ++n) {
        GlobMatC[threadIdx.x + 56 * n] = Results[n];
      }
    }
  }
}
void sgemm_NT_NT_m49_56_n9_9_k9_pps_834df01(const float ** A, int ExtraOffsetA, const float ** B, int ExtraOffsetB, float * C, int ExtraOffsetC, unsigned NumElements) {
  dim3 Block(64, 1, 1);
  dim3 Grid((NumElements + 1 - 1) / 1, 1, 1);
  kernel_sgemm_NT_NT_m49_56_n9_9_k9_pps_834df01<<<Grid,Block>>>(A, ExtraOffsetA, B, ExtraOffsetB, C, ExtraOffsetC, NumElements);
  CHECK_ERR;
}
__global__ void
    __launch_bounds__(64)
kernel_sgemm_NT_NT_m56_56_n9_56_k49_nsp_572e194(const float * __restrict__ A, int ExtraOffsetA, const float * B, int ExtraOffsetB, float ** C, int ExtraOffsetC, unsigned NumElements) {
  if ((threadIdx.y + blockDim.y * blockIdx.x) < NumElements) {
    const float* GlobMatA = &A[0 + ExtraOffsetA];
    const float* GlobMatB = &B[(threadIdx.y + blockDim.y * blockIdx.x) * 504 + 0 + ExtraOffsetB];
    float* GlobMatC = &C[(threadIdx.y + blockDim.y * blockIdx.x)][0 + ExtraOffsetC];
    __shared__ float Scratch[497];
    float* ShrMatB = &Scratch[threadIdx.y * 497];

    // using ExtendedPatchLoader
#pragma unroll
    for (int i = 0; i < 7; ++i) {
      ShrMatB[threadIdx.x + i * 64] = GlobMatB[threadIdx.x + i * 64];
    }
    if (threadIdx.x < 49) {
      ShrMatB[threadIdx.x + 448] = GlobMatB[threadIdx.x + 448];
    }
    __syncthreads();

    if (threadIdx.x < 56) {

      float Results[9] = {0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f};
      float Value;

      for (int k = 0; k < 49; ++k) {
        Value = GlobMatA[threadIdx.x + 56 * k];

#pragma unroll
        for (int n = 0; n < 9; ++n) {
          Results[n] += Value * ShrMatB[k + 56 * n];
        }
      }

#pragma unroll
      for (int n = 0; n < 9; ++n) {
        GlobMatC[threadIdx.x + 56 * n] = Results[n] + GlobMatC[threadIdx.x + 56 * n];
      }
    }
  }
}
void sgemm_NT_NT_m56_56_n9_56_k49_nsp_572e194(const float * __restrict__ A, int ExtraOffsetA, const float * B, int ExtraOffsetB, float ** C, int ExtraOffsetC, unsigned NumElements) {
  dim3 Block(64, 1, 1);
  dim3 Grid((NumElements + 1 - 1) / 1, 1, 1);
  kernel_sgemm_NT_NT_m56_56_n9_56_k49_nsp_572e194<<<Grid,Block>>>(A, ExtraOffsetA, B, ExtraOffsetB, C, ExtraOffsetC, NumElements);
  CHECK_ERR;
}
