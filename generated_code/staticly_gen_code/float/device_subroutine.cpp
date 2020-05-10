#include "gemmgen_aux.h"
__global__ void
    __launch_bounds__(64)
kernel_sgemm_NT_NT_m35_64_n9_9_k9_pps_91656fc(const float ** A, int ExtraOffsetA, const float ** B, int ExtraOffsetB, float * C, int ExtraOffsetC, unsigned NumElements) {
  if ((threadIdx.y + blockDim.y * blockIdx.x) < NumElements) {
    const float* GlobMatA = &A[(threadIdx.y + blockDim.y * blockIdx.x)][0 + ExtraOffsetA];
    const float* GlobMatB = &B[(threadIdx.y + blockDim.y * blockIdx.x)][0 + ExtraOffsetB];
    float* GlobMatC = &C[(threadIdx.y + blockDim.y * blockIdx.x) * 432 + 0 + ExtraOffsetC];
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
        Value = GlobMatA[threadIdx.x + 64 * k];

#pragma unroll
        for (int n = 0; n < 9; ++n) {
          Results[n] += Value * ShrMatB[k + 9 * n];
        }
      }

#pragma unroll
      for (int n = 0; n < 9; ++n) {
        GlobMatC[threadIdx.x + 48 * n] = Results[n];
      }
    }
  }
}
void sgemm_NT_NT_m35_64_n9_9_k9_pps_91656fc(const float ** A, int ExtraOffsetA, const float ** B, int ExtraOffsetB, float * C, int ExtraOffsetC, unsigned NumElements) {
  dim3 Block(64, 1, 1);
  dim3 Grid((NumElements + 1 - 1) / 1, 1, 1);
  kernel_sgemm_NT_NT_m35_64_n9_9_k9_pps_91656fc<<<Grid,Block>>>(A, ExtraOffsetA, B, ExtraOffsetB, C, ExtraOffsetC, NumElements);
  CHECK_ERR;
}
__global__ void
    __launch_bounds__(64)
kernel_sgemm_NT_NT_m53_64_n9_48_k35_nsp_49df72e(const float * __restrict__ A, int ExtraOffsetA, const float * B, int ExtraOffsetB, float ** C, int ExtraOffsetC, unsigned NumElements) {
  if ((threadIdx.y + blockDim.y * blockIdx.x) < NumElements) {
    const float* GlobMatA = &A[1 + ExtraOffsetA];
    const float* GlobMatB = &B[(threadIdx.y + blockDim.y * blockIdx.x) * 432 + 0 + ExtraOffsetB];
    float* GlobMatC = &C[(threadIdx.y + blockDim.y * blockIdx.x)][1 + ExtraOffsetC];
    __shared__ float Scratch[419];
    float* ShrMatB = &Scratch[threadIdx.y * 419];

    // using ExtendedPatchLoader
#pragma unroll
    for (int i = 0; i < 6; ++i) {
      ShrMatB[threadIdx.x + i * 64] = GlobMatB[threadIdx.x + i * 64];
    }
    if (threadIdx.x < 35) {
      ShrMatB[threadIdx.x + 384] = GlobMatB[threadIdx.x + 384];
    }
    __syncthreads();

    if (threadIdx.x < 53) {

      float Results[9] = {0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f};
      float Value;

      for (int k = 0; k < 35; ++k) {
        Value = GlobMatA[threadIdx.x + 64 * k];

#pragma unroll
        for (int n = 0; n < 9; ++n) {
          Results[n] += Value * ShrMatB[k + 48 * n];
        }
      }

#pragma unroll
      for (int n = 0; n < 9; ++n) {
        GlobMatC[threadIdx.x + 64 * n] = Results[n] + GlobMatC[threadIdx.x + 64 * n];
      }
    }
  }
}
void sgemm_NT_NT_m53_64_n9_48_k35_nsp_49df72e(const float * __restrict__ A, int ExtraOffsetA, const float * B, int ExtraOffsetB, float ** C, int ExtraOffsetC, unsigned NumElements) {
  dim3 Block(64, 1, 1);
  dim3 Grid((NumElements + 1 - 1) / 1, 1, 1);
  kernel_sgemm_NT_NT_m53_64_n9_48_k35_nsp_49df72e<<<Grid,Block>>>(A, ExtraOffsetA, B, ExtraOffsetB, C, ExtraOffsetC, NumElements);
  CHECK_ERR;
}
__global__ void
    __launch_bounds__(64)
kernel_sgemm_NT_NT_m54_64_n9_48_k35_nsp_13f0270(const float * __restrict__ A, int ExtraOffsetA, const float * B, int ExtraOffsetB, float ** C, int ExtraOffsetC, unsigned NumElements) {
  if ((threadIdx.y + blockDim.y * blockIdx.x) < NumElements) {
    const float* GlobMatA = &A[1 + ExtraOffsetA];
    const float* GlobMatB = &B[(threadIdx.y + blockDim.y * blockIdx.x) * 432 + 0 + ExtraOffsetB];
    float* GlobMatC = &C[(threadIdx.y + blockDim.y * blockIdx.x)][1 + ExtraOffsetC];
    __shared__ float Scratch[419];
    float* ShrMatB = &Scratch[threadIdx.y * 419];

    // using ExtendedPatchLoader
#pragma unroll
    for (int i = 0; i < 6; ++i) {
      ShrMatB[threadIdx.x + i * 64] = GlobMatB[threadIdx.x + i * 64];
    }
    if (threadIdx.x < 35) {
      ShrMatB[threadIdx.x + 384] = GlobMatB[threadIdx.x + 384];
    }
    __syncthreads();

    if (threadIdx.x < 54) {

      float Results[9] = {0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f};
      float Value;

      for (int k = 0; k < 35; ++k) {
        Value = GlobMatA[threadIdx.x + 64 * k];

#pragma unroll
        for (int n = 0; n < 9; ++n) {
          Results[n] += Value * ShrMatB[k + 48 * n];
        }
      }

#pragma unroll
      for (int n = 0; n < 9; ++n) {
        GlobMatC[threadIdx.x + 64 * n] = Results[n] + GlobMatC[threadIdx.x + 64 * n];
      }
    }
  }
}
void sgemm_NT_NT_m54_64_n9_48_k35_nsp_13f0270(const float * __restrict__ A, int ExtraOffsetA, const float * B, int ExtraOffsetB, float ** C, int ExtraOffsetC, unsigned NumElements) {
  dim3 Block(64, 1, 1);
  dim3 Grid((NumElements + 1 - 1) / 1, 1, 1);
  kernel_sgemm_NT_NT_m54_64_n9_48_k35_nsp_13f0270<<<Grid,Block>>>(A, ExtraOffsetA, B, ExtraOffsetB, C, ExtraOffsetC, NumElements);
  CHECK_ERR;
}
__global__ void
    __launch_bounds__(64)
kernel_sgemm_NT_NT_m55_64_n9_48_k35_nsp_6daa8d7(const float * __restrict__ A, int ExtraOffsetA, const float * B, int ExtraOffsetB, float ** C, int ExtraOffsetC, unsigned NumElements) {
  if ((threadIdx.y + blockDim.y * blockIdx.x) < NumElements) {
    const float* GlobMatA = &A[1 + ExtraOffsetA];
    const float* GlobMatB = &B[(threadIdx.y + blockDim.y * blockIdx.x) * 432 + 0 + ExtraOffsetB];
    float* GlobMatC = &C[(threadIdx.y + blockDim.y * blockIdx.x)][1 + ExtraOffsetC];
    __shared__ float Scratch[419];
    float* ShrMatB = &Scratch[threadIdx.y * 419];

    // using ExtendedPatchLoader
#pragma unroll
    for (int i = 0; i < 6; ++i) {
      ShrMatB[threadIdx.x + i * 64] = GlobMatB[threadIdx.x + i * 64];
    }
    if (threadIdx.x < 35) {
      ShrMatB[threadIdx.x + 384] = GlobMatB[threadIdx.x + 384];
    }
    __syncthreads();

    if (threadIdx.x < 55) {

      float Results[9] = {0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f};
      float Value;

      for (int k = 0; k < 35; ++k) {
        Value = GlobMatA[threadIdx.x + 64 * k];

#pragma unroll
        for (int n = 0; n < 9; ++n) {
          Results[n] += Value * ShrMatB[k + 48 * n];
        }
      }

#pragma unroll
      for (int n = 0; n < 9; ++n) {
        GlobMatC[threadIdx.x + 64 * n] = Results[n] + GlobMatC[threadIdx.x + 64 * n];
      }
    }
  }
}
void sgemm_NT_NT_m55_64_n9_48_k35_nsp_6daa8d7(const float * __restrict__ A, int ExtraOffsetA, const float * B, int ExtraOffsetB, float ** C, int ExtraOffsetC, unsigned NumElements) {
  dim3 Block(64, 1, 1);
  dim3 Grid((NumElements + 1 - 1) / 1, 1, 1);
  kernel_sgemm_NT_NT_m55_64_n9_48_k35_nsp_6daa8d7<<<Grid,Block>>>(A, ExtraOffsetA, B, ExtraOffsetB, C, ExtraOffsetC, NumElements);
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
kernel_sgemm_NT_NT_m56_64_n6_64_k56_npp_968d18f(const float * __restrict__ A, int ExtraOffsetA, const float ** B, int ExtraOffsetB, float ** C, int ExtraOffsetC, unsigned NumElements) {
  if ((threadIdx.y + blockDim.y * blockIdx.x) < NumElements) {
    const float* GlobMatA = &A[0 + ExtraOffsetA];
    const float* GlobMatB = &B[(threadIdx.y + blockDim.y * blockIdx.x)][0 + ExtraOffsetB];
    float* GlobMatC = &C[(threadIdx.y + blockDim.y * blockIdx.x)][0 + ExtraOffsetC];
    __shared__ float Scratch[376];
    float* ShrMatB = &Scratch[threadIdx.y * 376];

    // using ExtendedPatchLoader
#pragma unroll
    for (int i = 0; i < 5; ++i) {
      ShrMatB[threadIdx.x + i * 64] = GlobMatB[threadIdx.x + i * 64];
    }
    if (threadIdx.x < 56) {
      ShrMatB[threadIdx.x + 320] = GlobMatB[threadIdx.x + 320];
    }
    __syncthreads();

    if (threadIdx.x < 56) {

      float Results[6] = {0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f};
      float Value;

      for (int k = 0; k < 56; ++k) {
        Value = GlobMatA[threadIdx.x + 64 * k];

#pragma unroll
        for (int n = 0; n < 6; ++n) {
          Results[n] += Value * ShrMatB[k + 64 * n];
        }
      }

#pragma unroll
      for (int n = 0; n < 6; ++n) {
        GlobMatC[threadIdx.x + 64 * n] = Results[n];
      }
    }
  }
}
void sgemm_NT_NT_m56_64_n6_64_k56_npp_968d18f(const float * __restrict__ A, int ExtraOffsetA, const float ** B, int ExtraOffsetB, float ** C, int ExtraOffsetC, unsigned NumElements) {
  dim3 Block(64, 1, 1);
  dim3 Grid((NumElements + 1 - 1) / 1, 1, 1);
  kernel_sgemm_NT_NT_m56_64_n6_64_k56_npp_968d18f<<<Grid,Block>>>(A, ExtraOffsetA, B, ExtraOffsetB, C, ExtraOffsetC, NumElements);
  CHECK_ERR;
}
__global__ void
    __launch_bounds__(64)
kernel_sgemm_NT_NT_m56_64_n6_64_k56_npp_3250429(const float * __restrict__ A, int ExtraOffsetA, const float ** B, int ExtraOffsetB, float ** C, int ExtraOffsetC, unsigned NumElements) {
  if ((threadIdx.y + blockDim.y * blockIdx.x) < NumElements) {
    const float* GlobMatA = &A[0 + ExtraOffsetA];
    const float* GlobMatB = &B[(threadIdx.y + blockDim.y * blockIdx.x)][0 + ExtraOffsetB];
    float* GlobMatC = &C[(threadIdx.y + blockDim.y * blockIdx.x)][0 + ExtraOffsetC];
    __shared__ float Scratch[376];
    float* ShrMatB = &Scratch[threadIdx.y * 376];

    // using ExtendedPatchLoader
#pragma unroll
    for (int i = 0; i < 5; ++i) {
      ShrMatB[threadIdx.x + i * 64] = GlobMatB[threadIdx.x + i * 64];
    }
    if (threadIdx.x < 56) {
      ShrMatB[threadIdx.x + 320] = GlobMatB[threadIdx.x + 320];
    }
    __syncthreads();

    if (threadIdx.x < 56) {

      float Results[6] = {0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f};
      float Value;

      for (int k = 0; k < 56; ++k) {
        Value = GlobMatA[threadIdx.x + 64 * k];

#pragma unroll
        for (int n = 0; n < 6; ++n) {
          Results[n] += Value * ShrMatB[k + 64 * n];
        }
      }

#pragma unroll
      for (int n = 0; n < 6; ++n) {
        GlobMatC[threadIdx.x + 64 * n] = Results[n] + GlobMatC[threadIdx.x + 64 * n];
      }
    }
  }
}
void sgemm_NT_NT_m56_64_n6_64_k56_npp_3250429(const float * __restrict__ A, int ExtraOffsetA, const float ** B, int ExtraOffsetB, float ** C, int ExtraOffsetC, unsigned NumElements) {
  dim3 Block(64, 1, 1);
  dim3 Grid((NumElements + 1 - 1) / 1, 1, 1);
  kernel_sgemm_NT_NT_m56_64_n6_64_k56_npp_3250429<<<Grid,Block>>>(A, ExtraOffsetA, B, ExtraOffsetB, C, ExtraOffsetC, NumElements);
  CHECK_ERR;
}
__global__ void
    __launch_bounds__(32)
kernel_sgemm_NT_NT_m21_32_n9_64_k56_nps_3c76782(const float * __restrict__ A, int ExtraOffsetA, const float ** B, int ExtraOffsetB, float * C, int ExtraOffsetC, unsigned NumElements) {
  if ((threadIdx.y + blockDim.y * blockIdx.x) < NumElements) {
    const float* GlobMatA = &A[0 + ExtraOffsetA];
    const float* GlobMatB = &B[(threadIdx.y + blockDim.y * blockIdx.x)][0 + ExtraOffsetB];
    float* GlobMatC = &C[(threadIdx.y + blockDim.y * blockIdx.x) * 288 + 0 + ExtraOffsetC];
    __shared__ float Scratch[568];
    float* ShrMatB = &Scratch[threadIdx.y * 568];

    // using ExtendedPatchLoader
#pragma unroll
    for (int i = 0; i < 17; ++i) {
      ShrMatB[threadIdx.x + i * 32] = GlobMatB[threadIdx.x + i * 32];
    }
    if (threadIdx.x < 24) {
      ShrMatB[threadIdx.x + 544] = GlobMatB[threadIdx.x + 544];
    }
    __syncthreads();

    if (threadIdx.x < 21) {

      float Results[9] = {0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f};
      float Value;

      for (int k = 0; k < 56; ++k) {
        Value = GlobMatA[threadIdx.x + 32 * k];

#pragma unroll
        for (int n = 0; n < 9; ++n) {
          Results[n] += Value * ShrMatB[k + 64 * n];
        }
      }

#pragma unroll
      for (int n = 0; n < 9; ++n) {
        GlobMatC[threadIdx.x + 32 * n] = Results[n];
      }
    }
  }
}
void sgemm_NT_NT_m21_32_n9_64_k56_nps_3c76782(const float * __restrict__ A, int ExtraOffsetA, const float ** B, int ExtraOffsetB, float * C, int ExtraOffsetC, unsigned NumElements) {
  dim3 Block(32, 1, 1);
  dim3 Grid((NumElements + 1 - 1) / 1, 1, 1);
  kernel_sgemm_NT_NT_m21_32_n9_64_k56_nps_3c76782<<<Grid,Block>>>(A, ExtraOffsetA, B, ExtraOffsetB, C, ExtraOffsetC, NumElements);
  CHECK_ERR;
}
__global__ void
    __launch_bounds__(64)
kernel_sgemm_NT_NT_m21_32_n9_9_k9_sps_97c1866(const float * A, int ExtraOffsetA, const float ** B, int ExtraOffsetB, float * C, int ExtraOffsetC, unsigned NumElements) {
  if ((threadIdx.y + blockDim.y * blockIdx.x) < NumElements) {
    const float* GlobMatA = &A[(threadIdx.y + blockDim.y * blockIdx.x) * 288 + 0 + ExtraOffsetA];
    const float* GlobMatB = &B[(threadIdx.y + blockDim.y * blockIdx.x)][0 + ExtraOffsetB];
    float* GlobMatC = &C[(threadIdx.y + blockDim.y * blockIdx.x) * 288 + 0 + ExtraOffsetC];
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
        Value = GlobMatA[threadIdx.x + 32 * k];

#pragma unroll
        for (int n = 0; n < 9; ++n) {
          Results[n] += Value * ShrMatB[k + 9 * n];
        }
      }

#pragma unroll
      for (int n = 0; n < 9; ++n) {
        GlobMatC[threadIdx.x + 32 * n] = Results[n];
      }
    }
  }
}
void sgemm_NT_NT_m21_32_n9_9_k9_sps_97c1866(const float * A, int ExtraOffsetA, const float ** B, int ExtraOffsetB, float * C, int ExtraOffsetC, unsigned NumElements) {
  dim3 Block(32, 2, 1);
  dim3 Grid((NumElements + 2 - 1) / 2, 1, 1);
  kernel_sgemm_NT_NT_m21_32_n9_9_k9_sps_97c1866<<<Grid,Block>>>(A, ExtraOffsetA, B, ExtraOffsetB, C, ExtraOffsetC, NumElements);
  CHECK_ERR;
}
__global__ void
    __launch_bounds__(64)
kernel_sgemm_NT_NT_m56_64_n9_32_k21_nsp_0e604b1(const float * __restrict__ A, int ExtraOffsetA, const float * B, int ExtraOffsetB, float ** C, int ExtraOffsetC, unsigned NumElements) {
  if ((threadIdx.y + blockDim.y * blockIdx.x) < NumElements) {
    const float* GlobMatA = &A[0 + ExtraOffsetA];
    const float* GlobMatB = &B[(threadIdx.y + blockDim.y * blockIdx.x) * 288 + 0 + ExtraOffsetB];
    float* GlobMatC = &C[(threadIdx.y + blockDim.y * blockIdx.x)][0 + ExtraOffsetC];
    __shared__ float Scratch[277];
    float* ShrMatB = &Scratch[threadIdx.y * 277];

    // using ExtendedPatchLoader
#pragma unroll
    for (int i = 0; i < 4; ++i) {
      ShrMatB[threadIdx.x + i * 64] = GlobMatB[threadIdx.x + i * 64];
    }
    if (threadIdx.x < 21) {
      ShrMatB[threadIdx.x + 256] = GlobMatB[threadIdx.x + 256];
    }
    __syncthreads();

    if (threadIdx.x < 56) {

      float Results[9] = {0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f};
      float Value;

      for (int k = 0; k < 21; ++k) {
        Value = GlobMatA[threadIdx.x + 64 * k];

#pragma unroll
        for (int n = 0; n < 9; ++n) {
          Results[n] += Value * ShrMatB[k + 32 * n];
        }
      }

#pragma unroll
      for (int n = 0; n < 9; ++n) {
        GlobMatC[threadIdx.x + 64 * n] = Results[n] + GlobMatC[threadIdx.x + 64 * n];
      }
    }
  }
}
void sgemm_NT_NT_m56_64_n9_32_k21_nsp_0e604b1(const float * __restrict__ A, int ExtraOffsetA, const float * B, int ExtraOffsetB, float ** C, int ExtraOffsetC, unsigned NumElements) {
  dim3 Block(64, 1, 1);
  dim3 Grid((NumElements + 1 - 1) / 1, 1, 1);
  kernel_sgemm_NT_NT_m56_64_n9_32_k21_nsp_0e604b1<<<Grid,Block>>>(A, ExtraOffsetA, B, ExtraOffsetB, C, ExtraOffsetC, NumElements);
  CHECK_ERR;
}
__global__ void
    __launch_bounds__(64)
kernel_sgemm_NT_NT_m21_32_n9_32_k21_nss_3b271a8(const float * __restrict__ A, int ExtraOffsetA, const float * B, int ExtraOffsetB, float * C, int ExtraOffsetC, unsigned NumElements) {
  if ((threadIdx.y + blockDim.y * blockIdx.x) < NumElements) {
    const float* GlobMatA = &A[0 + ExtraOffsetA];
    const float* GlobMatB = &B[(threadIdx.y + blockDim.y * blockIdx.x) * 288 + 0 + ExtraOffsetB];
    float* GlobMatC = &C[(threadIdx.y + blockDim.y * blockIdx.x) * 288 + 0 + ExtraOffsetC];
    __shared__ float Scratch[554];
    float* ShrMatB = &Scratch[threadIdx.y * 277];

    // using ExtendedPatchLoader
#pragma unroll
    for (int i = 0; i < 8; ++i) {
      ShrMatB[threadIdx.x + i * 32] = GlobMatB[threadIdx.x + i * 32];
    }
    if (threadIdx.x < 21) {
      ShrMatB[threadIdx.x + 256] = GlobMatB[threadIdx.x + 256];
    }
    __syncthreads();

    if (threadIdx.x < 21) {

      float Results[9] = {0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f};
      float Value;

      for (int k = 0; k < 21; ++k) {
        Value = GlobMatA[threadIdx.x + 32 * k];

#pragma unroll
        for (int n = 0; n < 9; ++n) {
          Results[n] += Value * ShrMatB[k + 32 * n];
        }
      }

#pragma unroll
      for (int n = 0; n < 9; ++n) {
        GlobMatC[threadIdx.x + 32 * n] = Results[n];
      }
    }
  }
}
void sgemm_NT_NT_m21_32_n9_32_k21_nss_3b271a8(const float * __restrict__ A, int ExtraOffsetA, const float * B, int ExtraOffsetB, float * C, int ExtraOffsetC, unsigned NumElements) {
  dim3 Block(32, 2, 1);
  dim3 Grid((NumElements + 2 - 1) / 2, 1, 1);
  kernel_sgemm_NT_NT_m21_32_n9_32_k21_nss_3b271a8<<<Grid,Block>>>(A, ExtraOffsetA, B, ExtraOffsetB, C, ExtraOffsetC, NumElements);
  CHECK_ERR;
}
__global__ void
    __launch_bounds__(576)
kernel_scopyAddScale_m56_64_n9_64_pp_1869608(float Scale, const float ** A, int ExtraOffsetA, float ** B, int ExtraOffsetB, unsigned NumElements) {
  if ((threadIdx.z + blockDim.z * blockIdx.x) < NumElements) {
    const float* GlobMatA = &A[(threadIdx.z + blockDim.z * blockIdx.x)][ExtraOffsetA + 0 + threadIdx.y * 64];
    float* GlobMatB = &B[(threadIdx.z + blockDim.z * blockIdx.x)][ExtraOffsetB + 0 + threadIdx.y * 64];
    if (threadIdx.x < 56) {
      GlobMatB[threadIdx.x] = Scale * GlobMatA[threadIdx.x];
    }
  }
}
void scopyAddScale_m56_64_n9_64_pp_1869608(float Scale, const float ** A, int ExtraOffsetA, float ** B, int ExtraOffsetB, unsigned NumElements) {
  dim3 Block(64, 9, 1);
  dim3 Grid((NumElements + 1 - 1) / 1, 1, 1);
  kernel_scopyAddScale_m56_64_n9_64_pp_1869608<<<Grid,Block>>>(Scale, A, ExtraOffsetA, B, ExtraOffsetB, NumElements);
  CHECK_ERR;
}
__global__ void
    __launch_bounds__(576)
kernel_scopyAddScale_m35_64_n9_64_pp_8aadec7(float Scale, const float ** A, int ExtraOffsetA, float ** B, int ExtraOffsetB, unsigned NumElements) {
  if ((threadIdx.z + blockDim.z * blockIdx.x) < NumElements) {
    const float* GlobMatA = &A[(threadIdx.z + blockDim.z * blockIdx.x)][ExtraOffsetA + 0 + threadIdx.y * 48];
    float* GlobMatB = &B[(threadIdx.z + blockDim.z * blockIdx.x)][ExtraOffsetB + 0 + threadIdx.y * 64];
    if (threadIdx.x < 35) {
      GlobMatB[threadIdx.x] += Scale * GlobMatA[threadIdx.x];
    }
  }
}
void scopyAddScale_m35_64_n9_64_pp_8aadec7(float Scale, const float ** A, int ExtraOffsetA, float ** B, int ExtraOffsetB, unsigned NumElements) {
  dim3 Block(64, 9, 1);
  dim3 Grid((NumElements + 1 - 1) / 1, 1, 1);
  kernel_scopyAddScale_m35_64_n9_64_pp_8aadec7<<<Grid,Block>>>(Scale, A, ExtraOffsetA, B, ExtraOffsetB, NumElements);
  CHECK_ERR;
}
__global__ void
    __launch_bounds__(288)
kernel_scopyAddScale_m20_64_n9_64_pp_f68fa48(float Scale, const float ** A, int ExtraOffsetA, float ** B, int ExtraOffsetB, unsigned NumElements) {
  if ((threadIdx.z + blockDim.z * blockIdx.x) < NumElements) {
    const float* GlobMatA = &A[(threadIdx.z + blockDim.z * blockIdx.x)][ExtraOffsetA + 0 + threadIdx.y * 32];
    float* GlobMatB = &B[(threadIdx.z + blockDim.z * blockIdx.x)][ExtraOffsetB + 0 + threadIdx.y * 64];
    if (threadIdx.x < 20) {
      GlobMatB[threadIdx.x] += Scale * GlobMatA[threadIdx.x];
    }
  }
}
void scopyAddScale_m20_64_n9_64_pp_f68fa48(float Scale, const float ** A, int ExtraOffsetA, float ** B, int ExtraOffsetB, unsigned NumElements) {
  dim3 Block(32, 9, 1);
  dim3 Grid((NumElements + 1 - 1) / 1, 1, 1);
  kernel_scopyAddScale_m20_64_n9_64_pp_f68fa48<<<Grid,Block>>>(Scale, A, ExtraOffsetA, B, ExtraOffsetB, NumElements);
  CHECK_ERR;
}
__global__ void
    __launch_bounds__(288)
kernel_scopyAddScale_m10_64_n9_64_pp_5f5bfb7(float Scale, const float ** A, int ExtraOffsetA, float ** B, int ExtraOffsetB, unsigned NumElements) {
  if ((threadIdx.z + blockDim.z * blockIdx.x) < NumElements) {
    const float* GlobMatA = &A[(threadIdx.z + blockDim.z * blockIdx.x)][ExtraOffsetA + 0 + threadIdx.y * 16];
    float* GlobMatB = &B[(threadIdx.z + blockDim.z * blockIdx.x)][ExtraOffsetB + 0 + threadIdx.y * 64];
    if (threadIdx.x < 10) {
      GlobMatB[threadIdx.x] += Scale * GlobMatA[threadIdx.x];
    }
  }
}
void scopyAddScale_m10_64_n9_64_pp_5f5bfb7(float Scale, const float ** A, int ExtraOffsetA, float ** B, int ExtraOffsetB, unsigned NumElements) {
  dim3 Block(32, 9, 1);
  dim3 Grid((NumElements + 1 - 1) / 1, 1, 1);
  kernel_scopyAddScale_m10_64_n9_64_pp_5f5bfb7<<<Grid,Block>>>(Scale, A, ExtraOffsetA, B, ExtraOffsetB, NumElements);
  CHECK_ERR;
}
__global__ void
    __launch_bounds__(288)
kernel_scopyAddScale_m4_64_n9_64_pp_8b6c4ec(float Scale, const float ** A, int ExtraOffsetA, float ** B, int ExtraOffsetB, unsigned NumElements) {
  if ((threadIdx.z + blockDim.z * blockIdx.x) < NumElements) {
    const float* GlobMatA = &A[(threadIdx.z + blockDim.z * blockIdx.x)][ExtraOffsetA + 0 + threadIdx.y * 16];
    float* GlobMatB = &B[(threadIdx.z + blockDim.z * blockIdx.x)][ExtraOffsetB + 0 + threadIdx.y * 64];
    if (threadIdx.x < 4) {
      GlobMatB[threadIdx.x] += Scale * GlobMatA[threadIdx.x];
    }
  }
}
void scopyAddScale_m4_64_n9_64_pp_8b6c4ec(float Scale, const float ** A, int ExtraOffsetA, float ** B, int ExtraOffsetB, unsigned NumElements) {
  dim3 Block(32, 9, 1);
  dim3 Grid((NumElements + 1 - 1) / 1, 1, 1);
  kernel_scopyAddScale_m4_64_n9_64_pp_8b6c4ec<<<Grid,Block>>>(Scale, A, ExtraOffsetA, B, ExtraOffsetB, NumElements);
  CHECK_ERR;
}
__global__ void
    __launch_bounds__(288)
kernel_scopyAddScale_m1_64_n9_64_pp_ddfc1ea(float Scale, const float ** A, int ExtraOffsetA, float ** B, int ExtraOffsetB, unsigned NumElements) {
  if ((threadIdx.z + blockDim.z * blockIdx.x) < NumElements) {
    const float* GlobMatA = &A[(threadIdx.z + blockDim.z * blockIdx.x)][ExtraOffsetA + 0 + threadIdx.y * 16];
    float* GlobMatB = &B[(threadIdx.z + blockDim.z * blockIdx.x)][ExtraOffsetB + 0 + threadIdx.y * 64];
    if (threadIdx.x < 1) {
      GlobMatB[threadIdx.x] += Scale * GlobMatA[threadIdx.x];
    }
  }
}
void scopyAddScale_m1_64_n9_64_pp_ddfc1ea(float Scale, const float ** A, int ExtraOffsetA, float ** B, int ExtraOffsetB, unsigned NumElements) {
  dim3 Block(32, 9, 1);
  dim3 Grid((NumElements + 1 - 1) / 1, 1, 1);
  kernel_scopyAddScale_m1_64_n9_64_pp_ddfc1ea<<<Grid,Block>>>(Scale, A, ExtraOffsetA, B, ExtraOffsetB, NumElements);
  CHECK_ERR;
}
__global__ void
    __launch_bounds__(64)
kernel_sgemm_NT_NT_m35_48_n9_64_k53_nps_7539a56(const float * __restrict__ A, int ExtraOffsetA, const float ** B, int ExtraOffsetB, float * C, int ExtraOffsetC, unsigned NumElements) {
  if ((threadIdx.y + blockDim.y * blockIdx.x) < NumElements) {
    const float* GlobMatA = &A[48 + ExtraOffsetA];
    const float* GlobMatB = &B[(threadIdx.y + blockDim.y * blockIdx.x)][1 + ExtraOffsetB];
    float* GlobMatC = &C[(threadIdx.y + blockDim.y * blockIdx.x) * 432 + 0 + ExtraOffsetC];
    __shared__ float Scratch[565];
    float* ShrMatB = &Scratch[threadIdx.y * 565];

    // using ExtendedPatchLoader
#pragma unroll
    for (int i = 0; i < 8; ++i) {
      ShrMatB[threadIdx.x + i * 64] = GlobMatB[threadIdx.x + i * 64];
    }
    if (threadIdx.x < 53) {
      ShrMatB[threadIdx.x + 512] = GlobMatB[threadIdx.x + 512];
    }
    __syncthreads();

    if (threadIdx.x < 35) {

      float Results[9] = {0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f};
      float Value;

      for (int k = 0; k < 53; ++k) {
        Value = GlobMatA[threadIdx.x + 48 * k];

#pragma unroll
        for (int n = 0; n < 9; ++n) {
          Results[n] += Value * ShrMatB[k + 64 * n];
        }
      }

#pragma unroll
      for (int n = 0; n < 9; ++n) {
        GlobMatC[threadIdx.x + 48 * n] = Results[n];
      }
    }
  }
}
void sgemm_NT_NT_m35_48_n9_64_k53_nps_7539a56(const float * __restrict__ A, int ExtraOffsetA, const float ** B, int ExtraOffsetB, float * C, int ExtraOffsetC, unsigned NumElements) {
  dim3 Block(64, 1, 1);
  dim3 Grid((NumElements + 1 - 1) / 1, 1, 1);
  kernel_sgemm_NT_NT_m35_48_n9_64_k53_nps_7539a56<<<Grid,Block>>>(A, ExtraOffsetA, B, ExtraOffsetB, C, ExtraOffsetC, NumElements);
  CHECK_ERR;
}
__global__ void
    __launch_bounds__(64)
kernel_sgemm_NT_NT_m35_48_n9_9_k9_spp_bcf4f26(const float * A, int ExtraOffsetA, const float ** B, int ExtraOffsetB, float ** C, int ExtraOffsetC, unsigned NumElements) {
  if ((threadIdx.y + blockDim.y * blockIdx.x) < NumElements) {
    const float* GlobMatA = &A[(threadIdx.y + blockDim.y * blockIdx.x) * 432 + 0 + ExtraOffsetA];
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
        Value = GlobMatA[threadIdx.x + 48 * k];

#pragma unroll
        for (int n = 0; n < 9; ++n) {
          Results[n] += Value * ShrMatB[k + 9 * n];
        }
      }

#pragma unroll
      for (int n = 0; n < 9; ++n) {
        GlobMatC[threadIdx.x + 48 * n] = Results[n];
      }
    }
  }
}
void sgemm_NT_NT_m35_48_n9_9_k9_spp_bcf4f26(const float * A, int ExtraOffsetA, const float ** B, int ExtraOffsetB, float ** C, int ExtraOffsetC, unsigned NumElements) {
  dim3 Block(64, 1, 1);
  dim3 Grid((NumElements + 1 - 1) / 1, 1, 1);
  kernel_sgemm_NT_NT_m35_48_n9_9_k9_spp_bcf4f26<<<Grid,Block>>>(A, ExtraOffsetA, B, ExtraOffsetB, C, ExtraOffsetC, NumElements);
  CHECK_ERR;
}
__global__ void
    __launch_bounds__(64)
kernel_sgemm_NT_NT_m35_48_n9_64_k54_nps_bfcf3fd(const float * __restrict__ A, int ExtraOffsetA, const float ** B, int ExtraOffsetB, float * C, int ExtraOffsetC, unsigned NumElements) {
  if ((threadIdx.y + blockDim.y * blockIdx.x) < NumElements) {
    const float* GlobMatA = &A[48 + ExtraOffsetA];
    const float* GlobMatB = &B[(threadIdx.y + blockDim.y * blockIdx.x)][1 + ExtraOffsetB];
    float* GlobMatC = &C[(threadIdx.y + blockDim.y * blockIdx.x) * 432 + 0 + ExtraOffsetC];
    __shared__ float Scratch[566];
    float* ShrMatB = &Scratch[threadIdx.y * 566];

    // using ExtendedPatchLoader
#pragma unroll
    for (int i = 0; i < 8; ++i) {
      ShrMatB[threadIdx.x + i * 64] = GlobMatB[threadIdx.x + i * 64];
    }
    if (threadIdx.x < 54) {
      ShrMatB[threadIdx.x + 512] = GlobMatB[threadIdx.x + 512];
    }
    __syncthreads();

    if (threadIdx.x < 35) {

      float Results[9] = {0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f};
      float Value;

      for (int k = 0; k < 54; ++k) {
        Value = GlobMatA[threadIdx.x + 48 * k];

#pragma unroll
        for (int n = 0; n < 9; ++n) {
          Results[n] += Value * ShrMatB[k + 64 * n];
        }
      }

#pragma unroll
      for (int n = 0; n < 9; ++n) {
        GlobMatC[threadIdx.x + 48 * n] = Results[n];
      }
    }
  }
}
void sgemm_NT_NT_m35_48_n9_64_k54_nps_bfcf3fd(const float * __restrict__ A, int ExtraOffsetA, const float ** B, int ExtraOffsetB, float * C, int ExtraOffsetC, unsigned NumElements) {
  dim3 Block(64, 1, 1);
  dim3 Grid((NumElements + 1 - 1) / 1, 1, 1);
  kernel_sgemm_NT_NT_m35_48_n9_64_k54_nps_bfcf3fd<<<Grid,Block>>>(A, ExtraOffsetA, B, ExtraOffsetB, C, ExtraOffsetC, NumElements);
  CHECK_ERR;
}
__global__ void
    __launch_bounds__(64)
kernel_sgemm_NT_NT_m35_48_n9_9_k9_spp_d135d69(const float * A, int ExtraOffsetA, const float ** B, int ExtraOffsetB, float ** C, int ExtraOffsetC, unsigned NumElements) {
  if ((threadIdx.y + blockDim.y * blockIdx.x) < NumElements) {
    const float* GlobMatA = &A[(threadIdx.y + blockDim.y * blockIdx.x) * 432 + 0 + ExtraOffsetA];
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
        Value = GlobMatA[threadIdx.x + 48 * k];

#pragma unroll
        for (int n = 0; n < 9; ++n) {
          Results[n] += Value * ShrMatB[k + 9 * n];
        }
      }

#pragma unroll
      for (int n = 0; n < 9; ++n) {
        GlobMatC[threadIdx.x + 48 * n] = Results[n] + GlobMatC[threadIdx.x + 48 * n];
      }
    }
  }
}
void sgemm_NT_NT_m35_48_n9_9_k9_spp_d135d69(const float * A, int ExtraOffsetA, const float ** B, int ExtraOffsetB, float ** C, int ExtraOffsetC, unsigned NumElements) {
  dim3 Block(64, 1, 1);
  dim3 Grid((NumElements + 1 - 1) / 1, 1, 1);
  kernel_sgemm_NT_NT_m35_48_n9_9_k9_spp_d135d69<<<Grid,Block>>>(A, ExtraOffsetA, B, ExtraOffsetB, C, ExtraOffsetC, NumElements);
  CHECK_ERR;
}
__global__ void
    __launch_bounds__(64)
kernel_sgemm_NT_NT_m35_48_n9_64_k55_nps_72bd187(const float * __restrict__ A, int ExtraOffsetA, const float ** B, int ExtraOffsetB, float * C, int ExtraOffsetC, unsigned NumElements) {
  if ((threadIdx.y + blockDim.y * blockIdx.x) < NumElements) {
    const float* GlobMatA = &A[48 + ExtraOffsetA];
    const float* GlobMatB = &B[(threadIdx.y + blockDim.y * blockIdx.x)][1 + ExtraOffsetB];
    float* GlobMatC = &C[(threadIdx.y + blockDim.y * blockIdx.x) * 432 + 0 + ExtraOffsetC];
    __shared__ float Scratch[567];
    float* ShrMatB = &Scratch[threadIdx.y * 567];

    // using ExtendedPatchLoader
#pragma unroll
    for (int i = 0; i < 8; ++i) {
      ShrMatB[threadIdx.x + i * 64] = GlobMatB[threadIdx.x + i * 64];
    }
    if (threadIdx.x < 55) {
      ShrMatB[threadIdx.x + 512] = GlobMatB[threadIdx.x + 512];
    }
    __syncthreads();

    if (threadIdx.x < 35) {

      float Results[9] = {0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f};
      float Value;

      for (int k = 0; k < 55; ++k) {
        Value = GlobMatA[threadIdx.x + 48 * k];

#pragma unroll
        for (int n = 0; n < 9; ++n) {
          Results[n] += Value * ShrMatB[k + 64 * n];
        }
      }

#pragma unroll
      for (int n = 0; n < 9; ++n) {
        GlobMatC[threadIdx.x + 48 * n] = Results[n];
      }
    }
  }
}
void sgemm_NT_NT_m35_48_n9_64_k55_nps_72bd187(const float * __restrict__ A, int ExtraOffsetA, const float ** B, int ExtraOffsetB, float * C, int ExtraOffsetC, unsigned NumElements) {
  dim3 Block(64, 1, 1);
  dim3 Grid((NumElements + 1 - 1) / 1, 1, 1);
  kernel_sgemm_NT_NT_m35_48_n9_64_k55_nps_72bd187<<<Grid,Block>>>(A, ExtraOffsetA, B, ExtraOffsetB, C, ExtraOffsetC, NumElements);
  CHECK_ERR;
}
__global__ void
    __launch_bounds__(64)
kernel_sgemm_NT_NT_m20_48_n9_48_k32_nps_c879df6(const float * __restrict__ A, int ExtraOffsetA, const float ** B, int ExtraOffsetB, float * C, int ExtraOffsetC, unsigned NumElements) {
  if ((threadIdx.y + blockDim.y * blockIdx.x) < NumElements) {
    const float* GlobMatA = &A[48 + ExtraOffsetA];
    const float* GlobMatB = &B[(threadIdx.y + blockDim.y * blockIdx.x)][1 + ExtraOffsetB];
    float* GlobMatC = &C[(threadIdx.y + blockDim.y * blockIdx.x) * 288 + 0 + ExtraOffsetC];
    __shared__ float Scratch[576];
    float* ShrMatB = &Scratch[threadIdx.y * 288];

    // using ExactPatchLoader
#pragma unroll
    for (int i = 0; i < 9; ++i) {
      ShrMatB[threadIdx.x + 0 + i * 32] = GlobMatB[threadIdx.x + 0 + i * 48];
    }
    __syncthreads();

    if (threadIdx.x < 20) {

      float Results[9] = {0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f};
      float Value;

      for (int k = 0; k < 32; ++k) {
        Value = GlobMatA[threadIdx.x + 48 * k];

#pragma unroll
        for (int n = 0; n < 9; ++n) {
          Results[n] += Value * ShrMatB[k + 32 * n];
        }
      }

#pragma unroll
      for (int n = 0; n < 9; ++n) {
        GlobMatC[threadIdx.x + 32 * n] = Results[n];
      }
    }
  }
}
void sgemm_NT_NT_m20_48_n9_48_k32_nps_c879df6(const float * __restrict__ A, int ExtraOffsetA, const float ** B, int ExtraOffsetB, float * C, int ExtraOffsetC, unsigned NumElements) {
  dim3 Block(32, 2, 1);
  dim3 Grid((NumElements + 2 - 1) / 2, 1, 1);
  kernel_sgemm_NT_NT_m20_48_n9_48_k32_nps_c879df6<<<Grid,Block>>>(A, ExtraOffsetA, B, ExtraOffsetB, C, ExtraOffsetC, NumElements);
  CHECK_ERR;
}
__global__ void
    __launch_bounds__(64)
kernel_sgemm_NT_NT_m20_32_n9_9_k9_spp_7cafb26(const float * A, int ExtraOffsetA, const float ** B, int ExtraOffsetB, float ** C, int ExtraOffsetC, unsigned NumElements) {
  if ((threadIdx.y + blockDim.y * blockIdx.x) < NumElements) {
    const float* GlobMatA = &A[(threadIdx.y + blockDim.y * blockIdx.x) * 288 + 0 + ExtraOffsetA];
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
        Value = GlobMatA[threadIdx.x + 32 * k];

#pragma unroll
        for (int n = 0; n < 9; ++n) {
          Results[n] += Value * ShrMatB[k + 9 * n];
        }
      }

#pragma unroll
      for (int n = 0; n < 9; ++n) {
        GlobMatC[threadIdx.x + 32 * n] = Results[n];
      }
    }
  }
}
void sgemm_NT_NT_m20_32_n9_9_k9_spp_7cafb26(const float * A, int ExtraOffsetA, const float ** B, int ExtraOffsetB, float ** C, int ExtraOffsetC, unsigned NumElements) {
  dim3 Block(32, 2, 1);
  dim3 Grid((NumElements + 2 - 1) / 2, 1, 1);
  kernel_sgemm_NT_NT_m20_32_n9_9_k9_spp_7cafb26<<<Grid,Block>>>(A, ExtraOffsetA, B, ExtraOffsetB, C, ExtraOffsetC, NumElements);
  CHECK_ERR;
}
__global__ void
    __launch_bounds__(32)
kernel_sgemm_NT_NT_m20_48_n9_48_k33_nps_34e9903(const float * __restrict__ A, int ExtraOffsetA, const float ** B, int ExtraOffsetB, float * C, int ExtraOffsetC, unsigned NumElements) {
  if ((threadIdx.y + blockDim.y * blockIdx.x) < NumElements) {
    const float* GlobMatA = &A[48 + ExtraOffsetA];
    const float* GlobMatB = &B[(threadIdx.y + blockDim.y * blockIdx.x)][1 + ExtraOffsetB];
    float* GlobMatC = &C[(threadIdx.y + blockDim.y * blockIdx.x) * 288 + 0 + ExtraOffsetC];
    __shared__ float Scratch[417];
    float* ShrMatB = &Scratch[threadIdx.y * 417];

    // using ExtendedPatchLoader
#pragma unroll
    for (int i = 0; i < 13; ++i) {
      ShrMatB[threadIdx.x + i * 32] = GlobMatB[threadIdx.x + i * 32];
    }
    if (threadIdx.x < 1) {
      ShrMatB[threadIdx.x + 416] = GlobMatB[threadIdx.x + 416];
    }
    __syncthreads();

    if (threadIdx.x < 20) {

      float Results[9] = {0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f};
      float Value;

      for (int k = 0; k < 33; ++k) {
        Value = GlobMatA[threadIdx.x + 48 * k];

#pragma unroll
        for (int n = 0; n < 9; ++n) {
          Results[n] += Value * ShrMatB[k + 48 * n];
        }
      }

#pragma unroll
      for (int n = 0; n < 9; ++n) {
        GlobMatC[threadIdx.x + 32 * n] = Results[n];
      }
    }
  }
}
void sgemm_NT_NT_m20_48_n9_48_k33_nps_34e9903(const float * __restrict__ A, int ExtraOffsetA, const float ** B, int ExtraOffsetB, float * C, int ExtraOffsetC, unsigned NumElements) {
  dim3 Block(32, 1, 1);
  dim3 Grid((NumElements + 1 - 1) / 1, 1, 1);
  kernel_sgemm_NT_NT_m20_48_n9_48_k33_nps_34e9903<<<Grid,Block>>>(A, ExtraOffsetA, B, ExtraOffsetB, C, ExtraOffsetC, NumElements);
  CHECK_ERR;
}
__global__ void
    __launch_bounds__(64)
kernel_sgemm_NT_NT_m20_32_n9_9_k9_spp_51f73de(const float * A, int ExtraOffsetA, const float ** B, int ExtraOffsetB, float ** C, int ExtraOffsetC, unsigned NumElements) {
  if ((threadIdx.y + blockDim.y * blockIdx.x) < NumElements) {
    const float* GlobMatA = &A[(threadIdx.y + blockDim.y * blockIdx.x) * 288 + 0 + ExtraOffsetA];
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
        Value = GlobMatA[threadIdx.x + 32 * k];

#pragma unroll
        for (int n = 0; n < 9; ++n) {
          Results[n] += Value * ShrMatB[k + 9 * n];
        }
      }

#pragma unroll
      for (int n = 0; n < 9; ++n) {
        GlobMatC[threadIdx.x + 32 * n] = Results[n] + GlobMatC[threadIdx.x + 32 * n];
      }
    }
  }
}
void sgemm_NT_NT_m20_32_n9_9_k9_spp_51f73de(const float * A, int ExtraOffsetA, const float ** B, int ExtraOffsetB, float ** C, int ExtraOffsetC, unsigned NumElements) {
  dim3 Block(32, 2, 1);
  dim3 Grid((NumElements + 2 - 1) / 2, 1, 1);
  kernel_sgemm_NT_NT_m20_32_n9_9_k9_spp_51f73de<<<Grid,Block>>>(A, ExtraOffsetA, B, ExtraOffsetB, C, ExtraOffsetC, NumElements);
  CHECK_ERR;
}
__global__ void
    __launch_bounds__(32)
kernel_sgemm_NT_NT_m20_48_n9_48_k34_nps_c44f1e3(const float * __restrict__ A, int ExtraOffsetA, const float ** B, int ExtraOffsetB, float * C, int ExtraOffsetC, unsigned NumElements) {
  if ((threadIdx.y + blockDim.y * blockIdx.x) < NumElements) {
    const float* GlobMatA = &A[48 + ExtraOffsetA];
    const float* GlobMatB = &B[(threadIdx.y + blockDim.y * blockIdx.x)][1 + ExtraOffsetB];
    float* GlobMatC = &C[(threadIdx.y + blockDim.y * blockIdx.x) * 288 + 0 + ExtraOffsetC];
    __shared__ float Scratch[418];
    float* ShrMatB = &Scratch[threadIdx.y * 418];

    // using ExtendedPatchLoader
#pragma unroll
    for (int i = 0; i < 13; ++i) {
      ShrMatB[threadIdx.x + i * 32] = GlobMatB[threadIdx.x + i * 32];
    }
    if (threadIdx.x < 2) {
      ShrMatB[threadIdx.x + 416] = GlobMatB[threadIdx.x + 416];
    }
    __syncthreads();

    if (threadIdx.x < 20) {

      float Results[9] = {0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f};
      float Value;

      for (int k = 0; k < 34; ++k) {
        Value = GlobMatA[threadIdx.x + 48 * k];

#pragma unroll
        for (int n = 0; n < 9; ++n) {
          Results[n] += Value * ShrMatB[k + 48 * n];
        }
      }

#pragma unroll
      for (int n = 0; n < 9; ++n) {
        GlobMatC[threadIdx.x + 32 * n] = Results[n];
      }
    }
  }
}
void sgemm_NT_NT_m20_48_n9_48_k34_nps_c44f1e3(const float * __restrict__ A, int ExtraOffsetA, const float ** B, int ExtraOffsetB, float * C, int ExtraOffsetC, unsigned NumElements) {
  dim3 Block(32, 1, 1);
  dim3 Grid((NumElements + 1 - 1) / 1, 1, 1);
  kernel_sgemm_NT_NT_m20_48_n9_48_k34_nps_c44f1e3<<<Grid,Block>>>(A, ExtraOffsetA, B, ExtraOffsetB, C, ExtraOffsetC, NumElements);
  CHECK_ERR;
}
__global__ void
    __launch_bounds__(64)
kernel_sgemm_NT_NT_m10_48_n9_32_k17_nps_eeac906(const float * __restrict__ A, int ExtraOffsetA, const float ** B, int ExtraOffsetB, float * C, int ExtraOffsetC, unsigned NumElements) {
  if ((threadIdx.y + blockDim.y * blockIdx.x) < NumElements) {
    const float* GlobMatA = &A[48 + ExtraOffsetA];
    const float* GlobMatB = &B[(threadIdx.y + blockDim.y * blockIdx.x)][1 + ExtraOffsetB];
    float* GlobMatC = &C[(threadIdx.y + blockDim.y * blockIdx.x) * 144 + 0 + ExtraOffsetC];
    __shared__ float Scratch[546];
    float* ShrMatB = &Scratch[threadIdx.y * 273];

    // using ExtendedPatchLoader
#pragma unroll
    for (int i = 0; i < 8; ++i) {
      ShrMatB[threadIdx.x + i * 32] = GlobMatB[threadIdx.x + i * 32];
    }
    if (threadIdx.x < 17) {
      ShrMatB[threadIdx.x + 256] = GlobMatB[threadIdx.x + 256];
    }
    __syncthreads();

    if (threadIdx.x < 10) {

      float Results[9] = {0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f};
      float Value;

      for (int k = 0; k < 17; ++k) {
        Value = GlobMatA[threadIdx.x + 48 * k];

#pragma unroll
        for (int n = 0; n < 9; ++n) {
          Results[n] += Value * ShrMatB[k + 32 * n];
        }
      }

#pragma unroll
      for (int n = 0; n < 9; ++n) {
        GlobMatC[threadIdx.x + 16 * n] = Results[n];
      }
    }
  }
}
void sgemm_NT_NT_m10_48_n9_32_k17_nps_eeac906(const float * __restrict__ A, int ExtraOffsetA, const float ** B, int ExtraOffsetB, float * C, int ExtraOffsetC, unsigned NumElements) {
  dim3 Block(32, 2, 1);
  dim3 Grid((NumElements + 2 - 1) / 2, 1, 1);
  kernel_sgemm_NT_NT_m10_48_n9_32_k17_nps_eeac906<<<Grid,Block>>>(A, ExtraOffsetA, B, ExtraOffsetB, C, ExtraOffsetC, NumElements);
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
    __launch_bounds__(64)
kernel_sgemm_NT_NT_m10_48_n9_32_k18_nps_2125e35(const float * __restrict__ A, int ExtraOffsetA, const float ** B, int ExtraOffsetB, float * C, int ExtraOffsetC, unsigned NumElements) {
  if ((threadIdx.y + blockDim.y * blockIdx.x) < NumElements) {
    const float* GlobMatA = &A[48 + ExtraOffsetA];
    const float* GlobMatB = &B[(threadIdx.y + blockDim.y * blockIdx.x)][1 + ExtraOffsetB];
    float* GlobMatC = &C[(threadIdx.y + blockDim.y * blockIdx.x) * 144 + 0 + ExtraOffsetC];
    __shared__ float Scratch[548];
    float* ShrMatB = &Scratch[threadIdx.y * 274];

    // using ExtendedPatchLoader
#pragma unroll
    for (int i = 0; i < 8; ++i) {
      ShrMatB[threadIdx.x + i * 32] = GlobMatB[threadIdx.x + i * 32];
    }
    if (threadIdx.x < 18) {
      ShrMatB[threadIdx.x + 256] = GlobMatB[threadIdx.x + 256];
    }
    __syncthreads();

    if (threadIdx.x < 10) {

      float Results[9] = {0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f};
      float Value;

      for (int k = 0; k < 18; ++k) {
        Value = GlobMatA[threadIdx.x + 48 * k];

#pragma unroll
        for (int n = 0; n < 9; ++n) {
          Results[n] += Value * ShrMatB[k + 32 * n];
        }
      }

#pragma unroll
      for (int n = 0; n < 9; ++n) {
        GlobMatC[threadIdx.x + 16 * n] = Results[n];
      }
    }
  }
}
void sgemm_NT_NT_m10_48_n9_32_k18_nps_2125e35(const float * __restrict__ A, int ExtraOffsetA, const float ** B, int ExtraOffsetB, float * C, int ExtraOffsetC, unsigned NumElements) {
  dim3 Block(32, 2, 1);
  dim3 Grid((NumElements + 2 - 1) / 2, 1, 1);
  kernel_sgemm_NT_NT_m10_48_n9_32_k18_nps_2125e35<<<Grid,Block>>>(A, ExtraOffsetA, B, ExtraOffsetB, C, ExtraOffsetC, NumElements);
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
    __launch_bounds__(64)
kernel_sgemm_NT_NT_m10_48_n9_32_k19_nps_63bce85(const float * __restrict__ A, int ExtraOffsetA, const float ** B, int ExtraOffsetB, float * C, int ExtraOffsetC, unsigned NumElements) {
  if ((threadIdx.y + blockDim.y * blockIdx.x) < NumElements) {
    const float* GlobMatA = &A[48 + ExtraOffsetA];
    const float* GlobMatB = &B[(threadIdx.y + blockDim.y * blockIdx.x)][1 + ExtraOffsetB];
    float* GlobMatC = &C[(threadIdx.y + blockDim.y * blockIdx.x) * 144 + 0 + ExtraOffsetC];
    __shared__ float Scratch[550];
    float* ShrMatB = &Scratch[threadIdx.y * 275];

    // using ExtendedPatchLoader
#pragma unroll
    for (int i = 0; i < 8; ++i) {
      ShrMatB[threadIdx.x + i * 32] = GlobMatB[threadIdx.x + i * 32];
    }
    if (threadIdx.x < 19) {
      ShrMatB[threadIdx.x + 256] = GlobMatB[threadIdx.x + 256];
    }
    __syncthreads();

    if (threadIdx.x < 10) {

      float Results[9] = {0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f};
      float Value;

      for (int k = 0; k < 19; ++k) {
        Value = GlobMatA[threadIdx.x + 48 * k];

#pragma unroll
        for (int n = 0; n < 9; ++n) {
          Results[n] += Value * ShrMatB[k + 32 * n];
        }
      }

#pragma unroll
      for (int n = 0; n < 9; ++n) {
        GlobMatC[threadIdx.x + 16 * n] = Results[n];
      }
    }
  }
}
void sgemm_NT_NT_m10_48_n9_32_k19_nps_63bce85(const float * __restrict__ A, int ExtraOffsetA, const float ** B, int ExtraOffsetB, float * C, int ExtraOffsetC, unsigned NumElements) {
  dim3 Block(32, 2, 1);
  dim3 Grid((NumElements + 2 - 1) / 2, 1, 1);
  kernel_sgemm_NT_NT_m10_48_n9_32_k19_nps_63bce85<<<Grid,Block>>>(A, ExtraOffsetA, B, ExtraOffsetB, C, ExtraOffsetC, NumElements);
  CHECK_ERR;
}
__global__ void
    __launch_bounds__(96)
kernel_sgemm_NT_NT_m4_48_n9_16_k7_nps_dff4b53(const float * __restrict__ A, int ExtraOffsetA, const float ** B, int ExtraOffsetB, float * C, int ExtraOffsetC, unsigned NumElements) {
  if ((threadIdx.y + blockDim.y * blockIdx.x) < NumElements) {
    const float* GlobMatA = &A[48 + ExtraOffsetA];
    const float* GlobMatB = &B[(threadIdx.y + blockDim.y * blockIdx.x)][1 + ExtraOffsetB];
    float* GlobMatC = &C[(threadIdx.y + blockDim.y * blockIdx.x) * 144 + 0 + ExtraOffsetC];
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
        Value = GlobMatA[threadIdx.x + 48 * k];

#pragma unroll
        for (int n = 0; n < 9; ++n) {
          Results[n] += Value * ShrMatB[k + 16 * n];
        }
      }

#pragma unroll
      for (int n = 0; n < 9; ++n) {
        GlobMatC[threadIdx.x + 16 * n] = Results[n];
      }
    }
  }
}
void sgemm_NT_NT_m4_48_n9_16_k7_nps_dff4b53(const float * __restrict__ A, int ExtraOffsetA, const float ** B, int ExtraOffsetB, float * C, int ExtraOffsetC, unsigned NumElements) {
  dim3 Block(32, 3, 1);
  dim3 Grid((NumElements + 3 - 1) / 3, 1, 1);
  kernel_sgemm_NT_NT_m4_48_n9_16_k7_nps_dff4b53<<<Grid,Block>>>(A, ExtraOffsetA, B, ExtraOffsetB, C, ExtraOffsetC, NumElements);
  CHECK_ERR;
}
__global__ void
    __launch_bounds__(96)
kernel_sgemm_NT_NT_m4_16_n9_9_k9_spp_b787603(const float * A, int ExtraOffsetA, const float ** B, int ExtraOffsetB, float ** C, int ExtraOffsetC, unsigned NumElements) {
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

    if (threadIdx.x < 4) {

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
void sgemm_NT_NT_m4_16_n9_9_k9_spp_b787603(const float * A, int ExtraOffsetA, const float ** B, int ExtraOffsetB, float ** C, int ExtraOffsetC, unsigned NumElements) {
  dim3 Block(32, 3, 1);
  dim3 Grid((NumElements + 3 - 1) / 3, 1, 1);
  kernel_sgemm_NT_NT_m4_16_n9_9_k9_spp_b787603<<<Grid,Block>>>(A, ExtraOffsetA, B, ExtraOffsetB, C, ExtraOffsetC, NumElements);
  CHECK_ERR;
}
__global__ void
    __launch_bounds__(96)
kernel_sgemm_NT_NT_m4_48_n9_16_k8_nps_7c6bd3e(const float * __restrict__ A, int ExtraOffsetA, const float ** B, int ExtraOffsetB, float * C, int ExtraOffsetC, unsigned NumElements) {
  if ((threadIdx.y + blockDim.y * blockIdx.x) < NumElements) {
    const float* GlobMatA = &A[48 + ExtraOffsetA];
    const float* GlobMatB = &B[(threadIdx.y + blockDim.y * blockIdx.x)][1 + ExtraOffsetB];
    float* GlobMatC = &C[(threadIdx.y + blockDim.y * blockIdx.x) * 144 + 0 + ExtraOffsetC];
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
        Value = GlobMatA[threadIdx.x + 48 * k];

#pragma unroll
        for (int n = 0; n < 9; ++n) {
          Results[n] += Value * ShrMatB[k + 16 * n];
        }
      }

#pragma unroll
      for (int n = 0; n < 9; ++n) {
        GlobMatC[threadIdx.x + 16 * n] = Results[n];
      }
    }
  }
}
void sgemm_NT_NT_m4_48_n9_16_k8_nps_7c6bd3e(const float * __restrict__ A, int ExtraOffsetA, const float ** B, int ExtraOffsetB, float * C, int ExtraOffsetC, unsigned NumElements) {
  dim3 Block(32, 3, 1);
  dim3 Grid((NumElements + 3 - 1) / 3, 1, 1);
  kernel_sgemm_NT_NT_m4_48_n9_16_k8_nps_7c6bd3e<<<Grid,Block>>>(A, ExtraOffsetA, B, ExtraOffsetB, C, ExtraOffsetC, NumElements);
  CHECK_ERR;
}
__global__ void
    __launch_bounds__(96)
kernel_sgemm_NT_NT_m4_16_n9_9_k9_spp_fbdc4fe(const float * A, int ExtraOffsetA, const float ** B, int ExtraOffsetB, float ** C, int ExtraOffsetC, unsigned NumElements) {
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

    if (threadIdx.x < 4) {

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
void sgemm_NT_NT_m4_16_n9_9_k9_spp_fbdc4fe(const float * A, int ExtraOffsetA, const float ** B, int ExtraOffsetB, float ** C, int ExtraOffsetC, unsigned NumElements) {
  dim3 Block(32, 3, 1);
  dim3 Grid((NumElements + 3 - 1) / 3, 1, 1);
  kernel_sgemm_NT_NT_m4_16_n9_9_k9_spp_fbdc4fe<<<Grid,Block>>>(A, ExtraOffsetA, B, ExtraOffsetB, C, ExtraOffsetC, NumElements);
  CHECK_ERR;
}
__global__ void
    __launch_bounds__(96)
kernel_sgemm_NT_NT_m4_48_n9_16_k9_nps_83bc4d7(const float * __restrict__ A, int ExtraOffsetA, const float ** B, int ExtraOffsetB, float * C, int ExtraOffsetC, unsigned NumElements) {
  if ((threadIdx.y + blockDim.y * blockIdx.x) < NumElements) {
    const float* GlobMatA = &A[48 + ExtraOffsetA];
    const float* GlobMatB = &B[(threadIdx.y + blockDim.y * blockIdx.x)][1 + ExtraOffsetB];
    float* GlobMatC = &C[(threadIdx.y + blockDim.y * blockIdx.x) * 144 + 0 + ExtraOffsetC];
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
        Value = GlobMatA[threadIdx.x + 48 * k];

#pragma unroll
        for (int n = 0; n < 9; ++n) {
          Results[n] += Value * ShrMatB[k + 16 * n];
        }
      }

#pragma unroll
      for (int n = 0; n < 9; ++n) {
        GlobMatC[threadIdx.x + 16 * n] = Results[n];
      }
    }
  }
}
void sgemm_NT_NT_m4_48_n9_16_k9_nps_83bc4d7(const float * __restrict__ A, int ExtraOffsetA, const float ** B, int ExtraOffsetB, float * C, int ExtraOffsetC, unsigned NumElements) {
  dim3 Block(32, 3, 1);
  dim3 Grid((NumElements + 3 - 1) / 3, 1, 1);
  kernel_sgemm_NT_NT_m4_48_n9_16_k9_nps_83bc4d7<<<Grid,Block>>>(A, ExtraOffsetA, B, ExtraOffsetB, C, ExtraOffsetC, NumElements);
  CHECK_ERR;
}
__global__ void
    __launch_bounds__(96)
kernel_sgemm_NT_NT_m1_48_n9_16_k1_nps_7a8cc39(const float * __restrict__ A, int ExtraOffsetA, const float ** B, int ExtraOffsetB, float * C, int ExtraOffsetC, unsigned NumElements) {
  if ((threadIdx.y + blockDim.y * blockIdx.x) < NumElements) {
    const float* GlobMatA = &A[48 + ExtraOffsetA];
    const float* GlobMatB = &B[(threadIdx.y + blockDim.y * blockIdx.x)][1 + ExtraOffsetB];
    float* GlobMatC = &C[(threadIdx.y + blockDim.y * blockIdx.x) * 144 + 0 + ExtraOffsetC];
    __shared__ float Scratch[387];
    float* ShrMatB = &Scratch[threadIdx.y * 129];

    // using ExtendedPatchLoader
#pragma unroll
    for (int i = 0; i < 4; ++i) {
      ShrMatB[threadIdx.x + i * 32] = GlobMatB[threadIdx.x + i * 32];
    }
    if (threadIdx.x < 1) {
      ShrMatB[threadIdx.x + 128] = GlobMatB[threadIdx.x + 128];
    }
    __syncthreads();

    if (threadIdx.x < 1) {

      float Results[9] = {0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f};
      float Value;

      for (int k = 0; k < 1; ++k) {
        Value = GlobMatA[threadIdx.x + 48 * k];

#pragma unroll
        for (int n = 0; n < 9; ++n) {
          Results[n] += Value * ShrMatB[k + 16 * n];
        }
      }

#pragma unroll
      for (int n = 0; n < 9; ++n) {
        GlobMatC[threadIdx.x + 16 * n] = Results[n];
      }
    }
  }
}
void sgemm_NT_NT_m1_48_n9_16_k1_nps_7a8cc39(const float * __restrict__ A, int ExtraOffsetA, const float ** B, int ExtraOffsetB, float * C, int ExtraOffsetC, unsigned NumElements) {
  dim3 Block(32, 3, 1);
  dim3 Grid((NumElements + 3 - 1) / 3, 1, 1);
  kernel_sgemm_NT_NT_m1_48_n9_16_k1_nps_7a8cc39<<<Grid,Block>>>(A, ExtraOffsetA, B, ExtraOffsetB, C, ExtraOffsetC, NumElements);
  CHECK_ERR;
}
__global__ void
    __launch_bounds__(96)
kernel_sgemm_NT_NT_m1_16_n9_9_k9_spp_8da6d94(const float * A, int ExtraOffsetA, const float ** B, int ExtraOffsetB, float ** C, int ExtraOffsetC, unsigned NumElements) {
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

    if (threadIdx.x < 1) {

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
void sgemm_NT_NT_m1_16_n9_9_k9_spp_8da6d94(const float * A, int ExtraOffsetA, const float ** B, int ExtraOffsetB, float ** C, int ExtraOffsetC, unsigned NumElements) {
  dim3 Block(32, 3, 1);
  dim3 Grid((NumElements + 3 - 1) / 3, 1, 1);
  kernel_sgemm_NT_NT_m1_16_n9_9_k9_spp_8da6d94<<<Grid,Block>>>(A, ExtraOffsetA, B, ExtraOffsetB, C, ExtraOffsetC, NumElements);
  CHECK_ERR;
}
__global__ void
    __launch_bounds__(96)
kernel_sgemm_NT_NT_m1_48_n9_16_k2_nps_04350e1(const float * __restrict__ A, int ExtraOffsetA, const float ** B, int ExtraOffsetB, float * C, int ExtraOffsetC, unsigned NumElements) {
  if ((threadIdx.y + blockDim.y * blockIdx.x) < NumElements) {
    const float* GlobMatA = &A[48 + ExtraOffsetA];
    const float* GlobMatB = &B[(threadIdx.y + blockDim.y * blockIdx.x)][1 + ExtraOffsetB];
    float* GlobMatC = &C[(threadIdx.y + blockDim.y * blockIdx.x) * 144 + 0 + ExtraOffsetC];
    __shared__ float Scratch[390];
    float* ShrMatB = &Scratch[threadIdx.y * 130];

    // using ExtendedPatchLoader
#pragma unroll
    for (int i = 0; i < 4; ++i) {
      ShrMatB[threadIdx.x + i * 32] = GlobMatB[threadIdx.x + i * 32];
    }
    if (threadIdx.x < 2) {
      ShrMatB[threadIdx.x + 128] = GlobMatB[threadIdx.x + 128];
    }
    __syncthreads();

    if (threadIdx.x < 1) {

      float Results[9] = {0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f};
      float Value;

      for (int k = 0; k < 2; ++k) {
        Value = GlobMatA[threadIdx.x + 48 * k];

#pragma unroll
        for (int n = 0; n < 9; ++n) {
          Results[n] += Value * ShrMatB[k + 16 * n];
        }
      }

#pragma unroll
      for (int n = 0; n < 9; ++n) {
        GlobMatC[threadIdx.x + 16 * n] = Results[n];
      }
    }
  }
}
void sgemm_NT_NT_m1_48_n9_16_k2_nps_04350e1(const float * __restrict__ A, int ExtraOffsetA, const float ** B, int ExtraOffsetB, float * C, int ExtraOffsetC, unsigned NumElements) {
  dim3 Block(32, 3, 1);
  dim3 Grid((NumElements + 3 - 1) / 3, 1, 1);
  kernel_sgemm_NT_NT_m1_48_n9_16_k2_nps_04350e1<<<Grid,Block>>>(A, ExtraOffsetA, B, ExtraOffsetB, C, ExtraOffsetC, NumElements);
  CHECK_ERR;
}
__global__ void
    __launch_bounds__(96)
kernel_sgemm_NT_NT_m1_16_n9_9_k9_spp_fe6bf6f(const float * A, int ExtraOffsetA, const float ** B, int ExtraOffsetB, float ** C, int ExtraOffsetC, unsigned NumElements) {
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

    if (threadIdx.x < 1) {

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
void sgemm_NT_NT_m1_16_n9_9_k9_spp_fe6bf6f(const float * A, int ExtraOffsetA, const float ** B, int ExtraOffsetB, float ** C, int ExtraOffsetC, unsigned NumElements) {
  dim3 Block(32, 3, 1);
  dim3 Grid((NumElements + 3 - 1) / 3, 1, 1);
  kernel_sgemm_NT_NT_m1_16_n9_9_k9_spp_fe6bf6f<<<Grid,Block>>>(A, ExtraOffsetA, B, ExtraOffsetB, C, ExtraOffsetC, NumElements);
  CHECK_ERR;
}
__global__ void
    __launch_bounds__(96)
kernel_sgemm_NT_NT_m1_48_n9_16_k3_nps_c033cdc(const float * __restrict__ A, int ExtraOffsetA, const float ** B, int ExtraOffsetB, float * C, int ExtraOffsetC, unsigned NumElements) {
  if ((threadIdx.y + blockDim.y * blockIdx.x) < NumElements) {
    const float* GlobMatA = &A[48 + ExtraOffsetA];
    const float* GlobMatB = &B[(threadIdx.y + blockDim.y * blockIdx.x)][1 + ExtraOffsetB];
    float* GlobMatC = &C[(threadIdx.y + blockDim.y * blockIdx.x) * 144 + 0 + ExtraOffsetC];
    __shared__ float Scratch[393];
    float* ShrMatB = &Scratch[threadIdx.y * 131];

    // using ExtendedPatchLoader
#pragma unroll
    for (int i = 0; i < 4; ++i) {
      ShrMatB[threadIdx.x + i * 32] = GlobMatB[threadIdx.x + i * 32];
    }
    if (threadIdx.x < 3) {
      ShrMatB[threadIdx.x + 128] = GlobMatB[threadIdx.x + 128];
    }
    __syncthreads();

    if (threadIdx.x < 1) {

      float Results[9] = {0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f};
      float Value;

      for (int k = 0; k < 3; ++k) {
        Value = GlobMatA[threadIdx.x + 48 * k];

#pragma unroll
        for (int n = 0; n < 9; ++n) {
          Results[n] += Value * ShrMatB[k + 16 * n];
        }
      }

#pragma unroll
      for (int n = 0; n < 9; ++n) {
        GlobMatC[threadIdx.x + 16 * n] = Results[n];
      }
    }
  }
}
void sgemm_NT_NT_m1_48_n9_16_k3_nps_c033cdc(const float * __restrict__ A, int ExtraOffsetA, const float ** B, int ExtraOffsetB, float * C, int ExtraOffsetC, unsigned NumElements) {
  dim3 Block(32, 3, 1);
  dim3 Grid((NumElements + 3 - 1) / 3, 1, 1);
  kernel_sgemm_NT_NT_m1_48_n9_16_k3_nps_c033cdc<<<Grid,Block>>>(A, ExtraOffsetA, B, ExtraOffsetB, C, ExtraOffsetC, NumElements);
  CHECK_ERR;
}
__global__ void
    __launch_bounds__(64)
kernel_sgemm_NT_NT_m49_64_n9_64_k56_nps_29706f6(const float * __restrict__ A, int ExtraOffsetA, const float ** B, int ExtraOffsetB, float * C, int ExtraOffsetC, unsigned NumElements) {
  if ((threadIdx.y + blockDim.y * blockIdx.x) < NumElements) {
    const float* GlobMatA = &A[0 + ExtraOffsetA];
    const float* GlobMatB = &B[(threadIdx.y + blockDim.y * blockIdx.x)][0 + ExtraOffsetB];
    float* GlobMatC = &C[(threadIdx.y + blockDim.y * blockIdx.x) * 576 + 0 + ExtraOffsetC];
    __shared__ float Scratch[568];
    float* ShrMatB = &Scratch[threadIdx.y * 568];

    // using ExtendedPatchLoader
#pragma unroll
    for (int i = 0; i < 8; ++i) {
      ShrMatB[threadIdx.x + i * 64] = GlobMatB[threadIdx.x + i * 64];
    }
    if (threadIdx.x < 56) {
      ShrMatB[threadIdx.x + 512] = GlobMatB[threadIdx.x + 512];
    }
    __syncthreads();

    if (threadIdx.x < 49) {

      float Results[9] = {0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f};
      float Value;

      for (int k = 0; k < 56; ++k) {
        Value = GlobMatA[threadIdx.x + 64 * k];

#pragma unroll
        for (int n = 0; n < 9; ++n) {
          Results[n] += Value * ShrMatB[k + 64 * n];
        }
      }

#pragma unroll
      for (int n = 0; n < 9; ++n) {
        GlobMatC[threadIdx.x + 64 * n] = Results[n];
      }
    }
  }
}
void sgemm_NT_NT_m49_64_n9_64_k56_nps_29706f6(const float * __restrict__ A, int ExtraOffsetA, const float ** B, int ExtraOffsetB, float * C, int ExtraOffsetC, unsigned NumElements) {
  dim3 Block(64, 1, 1);
  dim3 Grid((NumElements + 1 - 1) / 1, 1, 1);
  kernel_sgemm_NT_NT_m49_64_n9_64_k56_nps_29706f6<<<Grid,Block>>>(A, ExtraOffsetA, B, ExtraOffsetB, C, ExtraOffsetC, NumElements);
  CHECK_ERR;
}
__global__ void
    __launch_bounds__(64)
kernel_sgemm_NT_NT_m49_64_n9_9_k9_spp_43ea6a7(const float * A, int ExtraOffsetA, const float ** B, int ExtraOffsetB, float ** C, int ExtraOffsetC, unsigned NumElements) {
  if ((threadIdx.y + blockDim.y * blockIdx.x) < NumElements) {
    const float* GlobMatA = &A[(threadIdx.y + blockDim.y * blockIdx.x) * 576 + 0 + ExtraOffsetA];
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
        Value = GlobMatA[threadIdx.x + 64 * k];

#pragma unroll
        for (int n = 0; n < 9; ++n) {
          Results[n] += Value * ShrMatB[k + 9 * n];
        }
      }

#pragma unroll
      for (int n = 0; n < 9; ++n) {
        GlobMatC[threadIdx.x + 64 * n] = Results[n];
      }
    }
  }
}
void sgemm_NT_NT_m49_64_n9_9_k9_spp_43ea6a7(const float * A, int ExtraOffsetA, const float ** B, int ExtraOffsetB, float ** C, int ExtraOffsetC, unsigned NumElements) {
  dim3 Block(64, 1, 1);
  dim3 Grid((NumElements + 1 - 1) / 1, 1, 1);
  kernel_sgemm_NT_NT_m49_64_n9_9_k9_spp_43ea6a7<<<Grid,Block>>>(A, ExtraOffsetA, B, ExtraOffsetB, C, ExtraOffsetC, NumElements);
  CHECK_ERR;
}
__global__ void
    __launch_bounds__(64)
kernel_sgemm_NT_NT_m49_64_n9_9_k9_spp_b4ee274(const float * A, int ExtraOffsetA, const float ** B, int ExtraOffsetB, float ** C, int ExtraOffsetC, unsigned NumElements) {
  if ((threadIdx.y + blockDim.y * blockIdx.x) < NumElements) {
    const float* GlobMatA = &A[(threadIdx.y + blockDim.y * blockIdx.x) * 576 + 0 + ExtraOffsetA];
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
        Value = GlobMatA[threadIdx.x + 64 * k];

#pragma unroll
        for (int n = 0; n < 9; ++n) {
          Results[n] += Value * ShrMatB[k + 9 * n];
        }
      }

#pragma unroll
      for (int n = 0; n < 9; ++n) {
        GlobMatC[threadIdx.x + 64 * n] = Results[n] + GlobMatC[threadIdx.x + 64 * n];
      }
    }
  }
}
void sgemm_NT_NT_m49_64_n9_9_k9_spp_b4ee274(const float * A, int ExtraOffsetA, const float ** B, int ExtraOffsetB, float ** C, int ExtraOffsetC, unsigned NumElements) {
  dim3 Block(64, 1, 1);
  dim3 Grid((NumElements + 1 - 1) / 1, 1, 1);
  kernel_sgemm_NT_NT_m49_64_n9_9_k9_spp_b4ee274<<<Grid,Block>>>(A, ExtraOffsetA, B, ExtraOffsetB, C, ExtraOffsetC, NumElements);
  CHECK_ERR;
}
__global__ void
    __launch_bounds__(64)
kernel_sgemm_NT_NT_m49_64_n9_9_k9_pps_4046c83(const float ** A, int ExtraOffsetA, const float ** B, int ExtraOffsetB, float * C, int ExtraOffsetC, unsigned NumElements) {
  if ((threadIdx.y + blockDim.y * blockIdx.x) < NumElements) {
    const float* GlobMatA = &A[(threadIdx.y + blockDim.y * blockIdx.x)][0 + ExtraOffsetA];
    const float* GlobMatB = &B[(threadIdx.y + blockDim.y * blockIdx.x)][0 + ExtraOffsetB];
    float* GlobMatC = &C[(threadIdx.y + blockDim.y * blockIdx.x) * 576 + 0 + ExtraOffsetC];
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
        Value = GlobMatA[threadIdx.x + 64 * k];

#pragma unroll
        for (int n = 0; n < 9; ++n) {
          Results[n] += Value * ShrMatB[k + 9 * n];
        }
      }

#pragma unroll
      for (int n = 0; n < 9; ++n) {
        GlobMatC[threadIdx.x + 64 * n] = Results[n];
      }
    }
  }
}
void sgemm_NT_NT_m49_64_n9_9_k9_pps_4046c83(const float ** A, int ExtraOffsetA, const float ** B, int ExtraOffsetB, float * C, int ExtraOffsetC, unsigned NumElements) {
  dim3 Block(64, 1, 1);
  dim3 Grid((NumElements + 1 - 1) / 1, 1, 1);
  kernel_sgemm_NT_NT_m49_64_n9_9_k9_pps_4046c83<<<Grid,Block>>>(A, ExtraOffsetA, B, ExtraOffsetB, C, ExtraOffsetC, NumElements);
  CHECK_ERR;
}
__global__ void
    __launch_bounds__(64)
kernel_sgemm_NT_NT_m56_64_n9_64_k49_nsp_01c641c(const float * __restrict__ A, int ExtraOffsetA, const float * B, int ExtraOffsetB, float ** C, int ExtraOffsetC, unsigned NumElements) {
  if ((threadIdx.y + blockDim.y * blockIdx.x) < NumElements) {
    const float* GlobMatA = &A[0 + ExtraOffsetA];
    const float* GlobMatB = &B[(threadIdx.y + blockDim.y * blockIdx.x) * 576 + 0 + ExtraOffsetB];
    float* GlobMatC = &C[(threadIdx.y + blockDim.y * blockIdx.x)][0 + ExtraOffsetC];
    __shared__ float Scratch[561];
    float* ShrMatB = &Scratch[threadIdx.y * 561];

    // using ExtendedPatchLoader
#pragma unroll
    for (int i = 0; i < 8; ++i) {
      ShrMatB[threadIdx.x + i * 64] = GlobMatB[threadIdx.x + i * 64];
    }
    if (threadIdx.x < 49) {
      ShrMatB[threadIdx.x + 512] = GlobMatB[threadIdx.x + 512];
    }
    __syncthreads();

    if (threadIdx.x < 56) {

      float Results[9] = {0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f};
      float Value;

      for (int k = 0; k < 49; ++k) {
        Value = GlobMatA[threadIdx.x + 64 * k];

#pragma unroll
        for (int n = 0; n < 9; ++n) {
          Results[n] += Value * ShrMatB[k + 64 * n];
        }
      }

#pragma unroll
      for (int n = 0; n < 9; ++n) {
        GlobMatC[threadIdx.x + 64 * n] = Results[n] + GlobMatC[threadIdx.x + 64 * n];
      }
    }
  }
}
void sgemm_NT_NT_m56_64_n9_64_k49_nsp_01c641c(const float * __restrict__ A, int ExtraOffsetA, const float * B, int ExtraOffsetB, float ** C, int ExtraOffsetC, unsigned NumElements) {
  dim3 Block(64, 1, 1);
  dim3 Grid((NumElements + 1 - 1) / 1, 1, 1);
  kernel_sgemm_NT_NT_m56_64_n9_64_k49_nsp_01c641c<<<Grid,Block>>>(A, ExtraOffsetA, B, ExtraOffsetB, C, ExtraOffsetC, NumElements);
  CHECK_ERR;
}
