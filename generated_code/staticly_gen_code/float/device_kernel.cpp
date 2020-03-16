#include "device.h"
#include <cassert>
#include <cstring>
#include <cstdlib>
#include "device_subroutine.h"
#include "device_kernel.h"
namespace device_gen_code {
  constexpr unsigned long const kernel::volume::NonZeroFlops;
  constexpr unsigned long const kernel::volume::HardwareFlops;
  void kernel::volume::execute() {
    assert(I != nullptr);
    assert(Q != nullptr);
    assert(kDivM(2) != nullptr);
    assert(kDivM(0) != nullptr);
    assert(kDivM(1) != nullptr);
    assert(star(2) != nullptr);
    assert(star(0) != nullptr);
    assert(star(1) != nullptr);
    assert(kernel::volume::num_elements != 0);
    device::DeviceInstance& device = device::DeviceInstance::getInstance();
    float *_tmp0, *_tmp2, *_tmp4;
    float *d_buffer0 = (float*)device.api->getStackMemory(360 * num_elements * sizeof(float));
    _tmp0 = d_buffer0;
    {
    unsigned offset_I = I_offset;
    unsigned offset_star = star_offset(0);
    unsigned offset__tmp0 = 360;
    
    
    device.gemm(device::ColMajor, device::NoTrans, device::NoTrans, 40, 9, 9, 1.0, I, 56, star(0), 9, 0.0, _tmp0, 40, offset_I, offset_star, offset__tmp0, num_elements);
    }
    {
    unsigned offset_kDivM = 0;
    unsigned offset__tmp0 = 360;
    unsigned offset_Q = Q_offset;
    
    
    device.gemm(device::ColMajor, device::NoTrans, device::NoTrans, 56, 9, 35, 1.0, kDivM(0), 56, _tmp0, 40, 1.0, Q, 56, offset_kDivM, offset__tmp0, offset_Q, num_elements);
    }
    _tmp2 = d_buffer0;
    {
    unsigned offset_I = I_offset;
    unsigned offset_star = star_offset(1);
    unsigned offset__tmp2 = 360;
    
    
    device.gemm(device::ColMajor, device::NoTrans, device::NoTrans, 40, 9, 9, 1.0, I, 56, star(1), 9, 0.0, _tmp2, 40, offset_I, offset_star, offset__tmp2, num_elements);
    }
    {
    unsigned offset_kDivM = 0;
    unsigned offset__tmp2 = 360;
    unsigned offset_Q = Q_offset;
    
    
    device.gemm(device::ColMajor, device::NoTrans, device::NoTrans, 56, 9, 35, 1.0, kDivM(1), 56, _tmp2, 40, 1.0, Q, 56, offset_kDivM, offset__tmp2, offset_Q, num_elements);
    }
    _tmp4 = d_buffer0;
    {
    unsigned offset_I = I_offset;
    unsigned offset_star = star_offset(2);
    unsigned offset__tmp4 = 360;
    
    
    device.gemm(device::ColMajor, device::NoTrans, device::NoTrans, 40, 9, 9, 1.0, I, 56, star(2), 9, 0.0, _tmp4, 40, offset_I, offset_star, offset__tmp4, num_elements);
    }
    {
    unsigned offset_kDivM = 0;
    unsigned offset__tmp4 = 360;
    unsigned offset_Q = Q_offset;
    
    
    device.gemm(device::ColMajor, device::NoTrans, device::NoTrans, 56, 9, 35, 1.0, kDivM(2), 56, _tmp4, 40, 1.0, Q, 56, offset_kDivM, offset__tmp4, offset_Q, num_elements);
    }
    device.api->popStackMemory();
  }
  constexpr unsigned long const kernel::rotateGodunovStateLocal::NonZeroFlops;
  constexpr unsigned long const kernel::rotateGodunovStateLocal::HardwareFlops;
  void kernel::rotateGodunovStateLocal::execute() {
    assert(QgodLocal != nullptr);
    assert(Tinv != nullptr);
    assert(godunovMatrix != nullptr);
    assert(kernel::rotateGodunovStateLocal::num_elements != 0);
    device::DeviceInstance& device = device::DeviceInstance::getInstance();
    {
    unsigned offset_Tinv = Tinv_offset;
    unsigned offset_QgodLocal = QgodLocal_offset;
    unsigned offset_godunovMatrix = godunovMatrix_offset;
    
    
    device.gemm(device::ColMajor, device::Trans, device::NoTrans, 9, 9, 9, 1.0, Tinv, 9, QgodLocal, 9, 0.0, godunovMatrix, 9, offset_Tinv, offset_QgodLocal, offset_godunovMatrix, num_elements);
    }
  }
  constexpr unsigned long const kernel::rotateGodunovStateNeighbor::NonZeroFlops;
  constexpr unsigned long const kernel::rotateGodunovStateNeighbor::HardwareFlops;
  void kernel::rotateGodunovStateNeighbor::execute() {
    assert(QgodNeighbor != nullptr);
    assert(Tinv != nullptr);
    assert(godunovMatrix != nullptr);
    assert(kernel::rotateGodunovStateNeighbor::num_elements != 0);
    device::DeviceInstance& device = device::DeviceInstance::getInstance();
    {
    unsigned offset_Tinv = Tinv_offset;
    unsigned offset_QgodNeighbor = QgodNeighbor_offset;
    unsigned offset_godunovMatrix = godunovMatrix_offset;
    
    
    device.gemm(device::ColMajor, device::Trans, device::NoTrans, 9, 9, 9, 1.0, Tinv, 9, QgodNeighbor, 9, 0.0, godunovMatrix, 9, offset_Tinv, offset_QgodNeighbor, offset_godunovMatrix, num_elements);
    }
  }
  constexpr unsigned long const kernel::rotateFluxMatrix::NonZeroFlops;
  constexpr unsigned long const kernel::rotateFluxMatrix::HardwareFlops;
  void kernel::rotateFluxMatrix::execute() {
    assert(!std::isnan(fluxScale));
    assert(T != nullptr);
    assert(fluxSolver != nullptr);
    assert(star(0) != nullptr);
    assert(kernel::rotateFluxMatrix::num_elements != 0);
    device::DeviceInstance& device = device::DeviceInstance::getInstance();
    {
    unsigned offset_star = star_offset(0);
    unsigned offset_T = T_offset;
    unsigned offset_fluxSolver = fluxSolver_offset;
    
    
    device.gemm(device::ColMajor, device::NoTrans, device::Trans, 9, 9, 9, fluxScale, star(0), 9, T, 9, 0.0, fluxSolver, 9, offset_star, offset_T, offset_fluxSolver, num_elements);
    }
  }
  constexpr unsigned long const kernel::plConvertToNodal::NonZeroFlops;
  constexpr unsigned long const kernel::plConvertToNodal::HardwareFlops;
  void kernel::plConvertToNodal::execute() {
    assert(QStress != nullptr);
    assert(QStressNodal != nullptr);
    assert(v != nullptr);
    assert(kernel::plConvertToNodal::num_elements != 0);
    device::DeviceInstance& device = device::DeviceInstance::getInstance();
    {
    unsigned offset_v = 0;
    unsigned offset_QStress = QStress_offset;
    unsigned offset_QStressNodal = QStressNodal_offset;
    
    
    device.gemm(device::ColMajor, device::NoTrans, device::NoTrans, 56, 6, 56, 1.0, v, 56, QStress, 56, 0.0, QStressNodal, 56, offset_v, offset_QStress, offset_QStressNodal, num_elements);
    }
  }
  constexpr unsigned long const kernel::plConvertToModal::NonZeroFlops;
  constexpr unsigned long const kernel::plConvertToModal::HardwareFlops;
  void kernel::plConvertToModal::execute() {
    assert(QStress != nullptr);
    assert(QStressNodal != nullptr);
    assert(vInv != nullptr);
    assert(kernel::plConvertToModal::num_elements != 0);
    device::DeviceInstance& device = device::DeviceInstance::getInstance();
    {
    unsigned offset_vInv = 0;
    unsigned offset_QStressNodal = QStressNodal_offset;
    unsigned offset_QStress = QStress_offset;
    
    
    device.gemm(device::ColMajor, device::NoTrans, device::NoTrans, 56, 6, 56, 1.0, vInv, 56, QStressNodal, 56, 1.0, QStress, 56, offset_vInv, offset_QStressNodal, offset_QStress, num_elements);
    }
  }
  constexpr unsigned long const kernel::localFlux::NonZeroFlops[];
  constexpr unsigned long const kernel::localFlux::HardwareFlops[];
  constexpr kernel::localFlux::member_function_ptr kernel::localFlux::ExecutePtrs[];
  void kernel::localFlux::execute0() {
    assert(AplusT != nullptr);
    assert(I != nullptr);
    assert(Q != nullptr);
    assert(fMrT(0) != nullptr);
    assert(rDivM(0) != nullptr);
    assert(kernel::localFlux::num_elements != 0);
    device::DeviceInstance& device = device::DeviceInstance::getInstance();
    float *_tmp0, *_tmp1;
    float *d_buffer0 = (float*)device.api->getStackMemory(216 * num_elements * sizeof(float));
    float *d_buffer1 = (float*)device.api->getStackMemory(216 * num_elements * sizeof(float));
    _tmp0 = d_buffer0;
    {
    unsigned offset_fMrT = 0;
    unsigned offset_I = I_offset;
    unsigned offset__tmp0 = 216;
    
    
    device.gemm(device::ColMajor, device::NoTrans, device::NoTrans, 24, 9, 56, 1.0, fMrT(0), 24, I, 56, 0.0, _tmp0, 24, offset_fMrT, offset_I, offset__tmp0, num_elements);
    }
    _tmp1 = d_buffer1;
    {
    unsigned offset__tmp0 = 216;
    unsigned offset_AplusT = AplusT_offset;
    unsigned offset__tmp1 = 216;
    
    
    device.gemm(device::ColMajor, device::NoTrans, device::NoTrans, 24, 9, 9, 1.0, _tmp0, 24, AplusT, 9, 0.0, _tmp1, 24, offset__tmp0, offset_AplusT, offset__tmp1, num_elements);
    }
    {
    unsigned offset_rDivM = 0;
    unsigned offset__tmp1 = 216;
    unsigned offset_Q = Q_offset;
    
    
    device.gemm(device::ColMajor, device::NoTrans, device::NoTrans, 56, 9, 21, 1.0, rDivM(0), 56, _tmp1, 24, 1.0, Q, 56, offset_rDivM, offset__tmp1, offset_Q, num_elements);
    }
    device.api->popStackMemory();
    device.api->popStackMemory();
  }
  void kernel::localFlux::execute1() {
    assert(AplusT != nullptr);
    assert(I != nullptr);
    assert(Q != nullptr);
    assert(fMrT(1) != nullptr);
    assert(rDivM(1) != nullptr);
    assert(kernel::localFlux::num_elements != 0);
    device::DeviceInstance& device = device::DeviceInstance::getInstance();
    float *_tmp0, *_tmp1;
    float *d_buffer0 = (float*)device.api->getStackMemory(216 * num_elements * sizeof(float));
    float *d_buffer1 = (float*)device.api->getStackMemory(216 * num_elements * sizeof(float));
    _tmp0 = d_buffer0;
    {
    unsigned offset_fMrT = 0;
    unsigned offset_I = I_offset;
    unsigned offset__tmp0 = 216;
    
    
    device.gemm(device::ColMajor, device::NoTrans, device::NoTrans, 24, 9, 56, 1.0, fMrT(1), 24, I, 56, 0.0, _tmp0, 24, offset_fMrT, offset_I, offset__tmp0, num_elements);
    }
    _tmp1 = d_buffer1;
    {
    unsigned offset__tmp0 = 216;
    unsigned offset_AplusT = AplusT_offset;
    unsigned offset__tmp1 = 216;
    
    
    device.gemm(device::ColMajor, device::NoTrans, device::NoTrans, 24, 9, 9, 1.0, _tmp0, 24, AplusT, 9, 0.0, _tmp1, 24, offset__tmp0, offset_AplusT, offset__tmp1, num_elements);
    }
    {
    unsigned offset_rDivM = 0;
    unsigned offset__tmp1 = 216;
    unsigned offset_Q = Q_offset;
    
    
    device.gemm(device::ColMajor, device::NoTrans, device::NoTrans, 56, 9, 21, 1.0, rDivM(1), 56, _tmp1, 24, 1.0, Q, 56, offset_rDivM, offset__tmp1, offset_Q, num_elements);
    }
    device.api->popStackMemory();
    device.api->popStackMemory();
  }
  void kernel::localFlux::execute2() {
    assert(AplusT != nullptr);
    assert(I != nullptr);
    assert(Q != nullptr);
    assert(fMrT(2) != nullptr);
    assert(rDivM(2) != nullptr);
    assert(kernel::localFlux::num_elements != 0);
    device::DeviceInstance& device = device::DeviceInstance::getInstance();
    float *_tmp0, *_tmp1;
    float *d_buffer0 = (float*)device.api->getStackMemory(216 * num_elements * sizeof(float));
    float *d_buffer1 = (float*)device.api->getStackMemory(216 * num_elements * sizeof(float));
    _tmp0 = d_buffer0;
    {
    unsigned offset_fMrT = 0;
    unsigned offset_I = I_offset;
    unsigned offset__tmp0 = 216;
    
    
    device.gemm(device::ColMajor, device::NoTrans, device::NoTrans, 24, 9, 56, 1.0, fMrT(2), 24, I, 56, 0.0, _tmp0, 24, offset_fMrT, offset_I, offset__tmp0, num_elements);
    }
    _tmp1 = d_buffer1;
    {
    unsigned offset__tmp0 = 216;
    unsigned offset_AplusT = AplusT_offset;
    unsigned offset__tmp1 = 216;
    
    
    device.gemm(device::ColMajor, device::NoTrans, device::NoTrans, 24, 9, 9, 1.0, _tmp0, 24, AplusT, 9, 0.0, _tmp1, 24, offset__tmp0, offset_AplusT, offset__tmp1, num_elements);
    }
    {
    unsigned offset_rDivM = 0;
    unsigned offset__tmp1 = 216;
    unsigned offset_Q = Q_offset;
    
    
    device.gemm(device::ColMajor, device::NoTrans, device::NoTrans, 56, 9, 21, 1.0, rDivM(2), 56, _tmp1, 24, 1.0, Q, 56, offset_rDivM, offset__tmp1, offset_Q, num_elements);
    }
    device.api->popStackMemory();
    device.api->popStackMemory();
  }
  void kernel::localFlux::execute3() {
    assert(AplusT != nullptr);
    assert(I != nullptr);
    assert(Q != nullptr);
    assert(fMrT(3) != nullptr);
    assert(rDivM(3) != nullptr);
    assert(kernel::localFlux::num_elements != 0);
    device::DeviceInstance& device = device::DeviceInstance::getInstance();
    float *_tmp0, *_tmp1;
    float *d_buffer0 = (float*)device.api->getStackMemory(216 * num_elements * sizeof(float));
    float *d_buffer1 = (float*)device.api->getStackMemory(216 * num_elements * sizeof(float));
    _tmp0 = d_buffer0;
    {
    unsigned offset_fMrT = 0;
    unsigned offset_I = I_offset;
    unsigned offset__tmp0 = 216;
    
    
    device.gemm(device::ColMajor, device::NoTrans, device::NoTrans, 24, 9, 56, 1.0, fMrT(3), 24, I, 56, 0.0, _tmp0, 24, offset_fMrT, offset_I, offset__tmp0, num_elements);
    }
    _tmp1 = d_buffer1;
    {
    unsigned offset__tmp0 = 216;
    unsigned offset_AplusT = AplusT_offset;
    unsigned offset__tmp1 = 216;
    
    
    device.gemm(device::ColMajor, device::NoTrans, device::NoTrans, 24, 9, 9, 1.0, _tmp0, 24, AplusT, 9, 0.0, _tmp1, 24, offset__tmp0, offset_AplusT, offset__tmp1, num_elements);
    }
    {
    unsigned offset_rDivM = 0;
    unsigned offset__tmp1 = 216;
    unsigned offset_Q = Q_offset;
    
    
    device.gemm(device::ColMajor, device::NoTrans, device::NoTrans, 56, 9, 21, 1.0, rDivM(3), 56, _tmp1, 24, 1.0, Q, 56, offset_rDivM, offset__tmp1, offset_Q, num_elements);
    }
    device.api->popStackMemory();
    device.api->popStackMemory();
  }
  constexpr unsigned long const kernel::neighboringFlux::NonZeroFlops[];
  constexpr unsigned long const kernel::neighboringFlux::HardwareFlops[];
  constexpr kernel::neighboringFlux::member_function_ptr kernel::neighboringFlux::ExecutePtrs[];
  void kernel::neighboringFlux::execute0() {
    assert(AminusT != nullptr);
    assert(I != nullptr);
    assert(Q != nullptr);
    assert(fP(0) != nullptr);
    assert(rDivM(0) != nullptr);
    assert(rT(0) != nullptr);
    assert(kernel::neighboringFlux::num_elements != 0);
    device::DeviceInstance& device = device::DeviceInstance::getInstance();
    float *_tmp0, *_tmp1, *_tmp2;
    float *d_buffer0 = (float*)device.api->getStackMemory(216 * num_elements * sizeof(float));
    float *d_buffer1 = (float*)device.api->getStackMemory(216 * num_elements * sizeof(float));
    _tmp0 = d_buffer0;
    {
    unsigned offset_rT = 0;
    unsigned offset_I = I_offset;
    unsigned offset__tmp0 = 216;
    
    
    device.gemm(device::ColMajor, device::NoTrans, device::NoTrans, 24, 9, 56, 1.0, rT(0), 24, I, 56, 0.0, _tmp0, 24, offset_rT, offset_I, offset__tmp0, num_elements);
    }
    _tmp1 = d_buffer1;
    {
    unsigned offset_fP = 0;
    unsigned offset__tmp0 = 216;
    unsigned offset__tmp1 = 216;
    
    
    device.gemm(device::ColMajor, device::NoTrans, device::NoTrans, 24, 9, 21, 1.0, fP(0), 24, _tmp0, 24, 0.0, _tmp1, 24, offset_fP, offset__tmp0, offset__tmp1, num_elements);
    }
    _tmp2 = d_buffer0;
    {
    unsigned offset__tmp1 = 216;
    unsigned offset_AminusT = AminusT_offset;
    unsigned offset__tmp2 = 216;
    
    
    device.gemm(device::ColMajor, device::NoTrans, device::NoTrans, 24, 9, 9, 1.0, _tmp1, 24, AminusT, 9, 0.0, _tmp2, 24, offset__tmp1, offset_AminusT, offset__tmp2, num_elements);
    }
    {
    unsigned offset_rDivM = 0;
    unsigned offset__tmp2 = 216;
    unsigned offset_Q = Q_offset;
    
    
    device.gemm(device::ColMajor, device::NoTrans, device::NoTrans, 56, 9, 21, 1.0, rDivM(0), 56, _tmp2, 24, 1.0, Q, 56, offset_rDivM, offset__tmp2, offset_Q, num_elements);
    }
    device.api->popStackMemory();
    device.api->popStackMemory();
  }
  void kernel::neighboringFlux::execute1() {
    assert(AminusT != nullptr);
    assert(I != nullptr);
    assert(Q != nullptr);
    assert(fP(1) != nullptr);
    assert(rDivM(0) != nullptr);
    assert(rT(0) != nullptr);
    assert(kernel::neighboringFlux::num_elements != 0);
    device::DeviceInstance& device = device::DeviceInstance::getInstance();
    float *_tmp0, *_tmp1, *_tmp2;
    float *d_buffer0 = (float*)device.api->getStackMemory(216 * num_elements * sizeof(float));
    float *d_buffer1 = (float*)device.api->getStackMemory(216 * num_elements * sizeof(float));
    _tmp0 = d_buffer0;
    {
    unsigned offset_rT = 0;
    unsigned offset_I = I_offset;
    unsigned offset__tmp0 = 216;
    
    
    device.gemm(device::ColMajor, device::NoTrans, device::NoTrans, 24, 9, 56, 1.0, rT(0), 24, I, 56, 0.0, _tmp0, 24, offset_rT, offset_I, offset__tmp0, num_elements);
    }
    _tmp1 = d_buffer1;
    {
    unsigned offset_fP = 0;
    unsigned offset__tmp0 = 216;
    unsigned offset__tmp1 = 216;
    
    
    device.gemm(device::ColMajor, device::NoTrans, device::NoTrans, 24, 9, 21, 1.0, fP(1), 24, _tmp0, 24, 0.0, _tmp1, 24, offset_fP, offset__tmp0, offset__tmp1, num_elements);
    }
    _tmp2 = d_buffer0;
    {
    unsigned offset__tmp1 = 216;
    unsigned offset_AminusT = AminusT_offset;
    unsigned offset__tmp2 = 216;
    
    
    device.gemm(device::ColMajor, device::NoTrans, device::NoTrans, 24, 9, 9, 1.0, _tmp1, 24, AminusT, 9, 0.0, _tmp2, 24, offset__tmp1, offset_AminusT, offset__tmp2, num_elements);
    }
    {
    unsigned offset_rDivM = 0;
    unsigned offset__tmp2 = 216;
    unsigned offset_Q = Q_offset;
    
    
    device.gemm(device::ColMajor, device::NoTrans, device::NoTrans, 56, 9, 21, 1.0, rDivM(0), 56, _tmp2, 24, 1.0, Q, 56, offset_rDivM, offset__tmp2, offset_Q, num_elements);
    }
    device.api->popStackMemory();
    device.api->popStackMemory();
  }
  void kernel::neighboringFlux::execute2() {
    assert(AminusT != nullptr);
    assert(I != nullptr);
    assert(Q != nullptr);
    assert(fP(2) != nullptr);
    assert(rDivM(0) != nullptr);
    assert(rT(0) != nullptr);
    assert(kernel::neighboringFlux::num_elements != 0);
    device::DeviceInstance& device = device::DeviceInstance::getInstance();
    float *_tmp0, *_tmp1, *_tmp2;
    float *d_buffer0 = (float*)device.api->getStackMemory(216 * num_elements * sizeof(float));
    float *d_buffer1 = (float*)device.api->getStackMemory(216 * num_elements * sizeof(float));
    _tmp0 = d_buffer0;
    {
    unsigned offset_rT = 0;
    unsigned offset_I = I_offset;
    unsigned offset__tmp0 = 216;
    
    
    device.gemm(device::ColMajor, device::NoTrans, device::NoTrans, 24, 9, 56, 1.0, rT(0), 24, I, 56, 0.0, _tmp0, 24, offset_rT, offset_I, offset__tmp0, num_elements);
    }
    _tmp1 = d_buffer1;
    {
    unsigned offset_fP = 0;
    unsigned offset__tmp0 = 216;
    unsigned offset__tmp1 = 216;
    
    
    device.gemm(device::ColMajor, device::NoTrans, device::NoTrans, 24, 9, 21, 1.0, fP(2), 24, _tmp0, 24, 0.0, _tmp1, 24, offset_fP, offset__tmp0, offset__tmp1, num_elements);
    }
    _tmp2 = d_buffer0;
    {
    unsigned offset__tmp1 = 216;
    unsigned offset_AminusT = AminusT_offset;
    unsigned offset__tmp2 = 216;
    
    
    device.gemm(device::ColMajor, device::NoTrans, device::NoTrans, 24, 9, 9, 1.0, _tmp1, 24, AminusT, 9, 0.0, _tmp2, 24, offset__tmp1, offset_AminusT, offset__tmp2, num_elements);
    }
    {
    unsigned offset_rDivM = 0;
    unsigned offset__tmp2 = 216;
    unsigned offset_Q = Q_offset;
    
    
    device.gemm(device::ColMajor, device::NoTrans, device::NoTrans, 56, 9, 21, 1.0, rDivM(0), 56, _tmp2, 24, 1.0, Q, 56, offset_rDivM, offset__tmp2, offset_Q, num_elements);
    }
    device.api->popStackMemory();
    device.api->popStackMemory();
  }
  void kernel::neighboringFlux::execute3() {
    assert(AminusT != nullptr);
    assert(I != nullptr);
    assert(Q != nullptr);
    assert(fP(0) != nullptr);
    assert(rDivM(0) != nullptr);
    assert(rT(1) != nullptr);
    assert(kernel::neighboringFlux::num_elements != 0);
    device::DeviceInstance& device = device::DeviceInstance::getInstance();
    float *_tmp0, *_tmp1, *_tmp2;
    float *d_buffer0 = (float*)device.api->getStackMemory(216 * num_elements * sizeof(float));
    float *d_buffer1 = (float*)device.api->getStackMemory(216 * num_elements * sizeof(float));
    _tmp0 = d_buffer0;
    {
    unsigned offset_rT = 0;
    unsigned offset_I = I_offset;
    unsigned offset__tmp0 = 216;
    
    
    device.gemm(device::ColMajor, device::NoTrans, device::NoTrans, 24, 9, 56, 1.0, rT(1), 24, I, 56, 0.0, _tmp0, 24, offset_rT, offset_I, offset__tmp0, num_elements);
    }
    _tmp1 = d_buffer1;
    {
    unsigned offset_fP = 0;
    unsigned offset__tmp0 = 216;
    unsigned offset__tmp1 = 216;
    
    
    device.gemm(device::ColMajor, device::NoTrans, device::NoTrans, 24, 9, 21, 1.0, fP(0), 24, _tmp0, 24, 0.0, _tmp1, 24, offset_fP, offset__tmp0, offset__tmp1, num_elements);
    }
    _tmp2 = d_buffer0;
    {
    unsigned offset__tmp1 = 216;
    unsigned offset_AminusT = AminusT_offset;
    unsigned offset__tmp2 = 216;
    
    
    device.gemm(device::ColMajor, device::NoTrans, device::NoTrans, 24, 9, 9, 1.0, _tmp1, 24, AminusT, 9, 0.0, _tmp2, 24, offset__tmp1, offset_AminusT, offset__tmp2, num_elements);
    }
    {
    unsigned offset_rDivM = 0;
    unsigned offset__tmp2 = 216;
    unsigned offset_Q = Q_offset;
    
    
    device.gemm(device::ColMajor, device::NoTrans, device::NoTrans, 56, 9, 21, 1.0, rDivM(0), 56, _tmp2, 24, 1.0, Q, 56, offset_rDivM, offset__tmp2, offset_Q, num_elements);
    }
    device.api->popStackMemory();
    device.api->popStackMemory();
  }
  void kernel::neighboringFlux::execute4() {
    assert(AminusT != nullptr);
    assert(I != nullptr);
    assert(Q != nullptr);
    assert(fP(1) != nullptr);
    assert(rDivM(0) != nullptr);
    assert(rT(1) != nullptr);
    assert(kernel::neighboringFlux::num_elements != 0);
    device::DeviceInstance& device = device::DeviceInstance::getInstance();
    float *_tmp0, *_tmp1, *_tmp2;
    float *d_buffer0 = (float*)device.api->getStackMemory(216 * num_elements * sizeof(float));
    float *d_buffer1 = (float*)device.api->getStackMemory(216 * num_elements * sizeof(float));
    _tmp0 = d_buffer0;
    {
    unsigned offset_rT = 0;
    unsigned offset_I = I_offset;
    unsigned offset__tmp0 = 216;
    
    
    device.gemm(device::ColMajor, device::NoTrans, device::NoTrans, 24, 9, 56, 1.0, rT(1), 24, I, 56, 0.0, _tmp0, 24, offset_rT, offset_I, offset__tmp0, num_elements);
    }
    _tmp1 = d_buffer1;
    {
    unsigned offset_fP = 0;
    unsigned offset__tmp0 = 216;
    unsigned offset__tmp1 = 216;
    
    
    device.gemm(device::ColMajor, device::NoTrans, device::NoTrans, 24, 9, 21, 1.0, fP(1), 24, _tmp0, 24, 0.0, _tmp1, 24, offset_fP, offset__tmp0, offset__tmp1, num_elements);
    }
    _tmp2 = d_buffer0;
    {
    unsigned offset__tmp1 = 216;
    unsigned offset_AminusT = AminusT_offset;
    unsigned offset__tmp2 = 216;
    
    
    device.gemm(device::ColMajor, device::NoTrans, device::NoTrans, 24, 9, 9, 1.0, _tmp1, 24, AminusT, 9, 0.0, _tmp2, 24, offset__tmp1, offset_AminusT, offset__tmp2, num_elements);
    }
    {
    unsigned offset_rDivM = 0;
    unsigned offset__tmp2 = 216;
    unsigned offset_Q = Q_offset;
    
    
    device.gemm(device::ColMajor, device::NoTrans, device::NoTrans, 56, 9, 21, 1.0, rDivM(0), 56, _tmp2, 24, 1.0, Q, 56, offset_rDivM, offset__tmp2, offset_Q, num_elements);
    }
    device.api->popStackMemory();
    device.api->popStackMemory();
  }
  void kernel::neighboringFlux::execute5() {
    assert(AminusT != nullptr);
    assert(I != nullptr);
    assert(Q != nullptr);
    assert(fP(2) != nullptr);
    assert(rDivM(0) != nullptr);
    assert(rT(1) != nullptr);
    assert(kernel::neighboringFlux::num_elements != 0);
    device::DeviceInstance& device = device::DeviceInstance::getInstance();
    float *_tmp0, *_tmp1, *_tmp2;
    float *d_buffer0 = (float*)device.api->getStackMemory(216 * num_elements * sizeof(float));
    float *d_buffer1 = (float*)device.api->getStackMemory(216 * num_elements * sizeof(float));
    _tmp0 = d_buffer0;
    {
    unsigned offset_rT = 0;
    unsigned offset_I = I_offset;
    unsigned offset__tmp0 = 216;
    
    
    device.gemm(device::ColMajor, device::NoTrans, device::NoTrans, 24, 9, 56, 1.0, rT(1), 24, I, 56, 0.0, _tmp0, 24, offset_rT, offset_I, offset__tmp0, num_elements);
    }
    _tmp1 = d_buffer1;
    {
    unsigned offset_fP = 0;
    unsigned offset__tmp0 = 216;
    unsigned offset__tmp1 = 216;
    
    
    device.gemm(device::ColMajor, device::NoTrans, device::NoTrans, 24, 9, 21, 1.0, fP(2), 24, _tmp0, 24, 0.0, _tmp1, 24, offset_fP, offset__tmp0, offset__tmp1, num_elements);
    }
    _tmp2 = d_buffer0;
    {
    unsigned offset__tmp1 = 216;
    unsigned offset_AminusT = AminusT_offset;
    unsigned offset__tmp2 = 216;
    
    
    device.gemm(device::ColMajor, device::NoTrans, device::NoTrans, 24, 9, 9, 1.0, _tmp1, 24, AminusT, 9, 0.0, _tmp2, 24, offset__tmp1, offset_AminusT, offset__tmp2, num_elements);
    }
    {
    unsigned offset_rDivM = 0;
    unsigned offset__tmp2 = 216;
    unsigned offset_Q = Q_offset;
    
    
    device.gemm(device::ColMajor, device::NoTrans, device::NoTrans, 56, 9, 21, 1.0, rDivM(0), 56, _tmp2, 24, 1.0, Q, 56, offset_rDivM, offset__tmp2, offset_Q, num_elements);
    }
    device.api->popStackMemory();
    device.api->popStackMemory();
  }
  void kernel::neighboringFlux::execute6() {
    assert(AminusT != nullptr);
    assert(I != nullptr);
    assert(Q != nullptr);
    assert(fP(0) != nullptr);
    assert(rDivM(0) != nullptr);
    assert(rT(2) != nullptr);
    assert(kernel::neighboringFlux::num_elements != 0);
    device::DeviceInstance& device = device::DeviceInstance::getInstance();
    float *_tmp0, *_tmp1, *_tmp2;
    float *d_buffer0 = (float*)device.api->getStackMemory(216 * num_elements * sizeof(float));
    float *d_buffer1 = (float*)device.api->getStackMemory(216 * num_elements * sizeof(float));
    _tmp0 = d_buffer0;
    {
    unsigned offset_rT = 0;
    unsigned offset_I = I_offset;
    unsigned offset__tmp0 = 216;
    
    
    device.gemm(device::ColMajor, device::NoTrans, device::NoTrans, 24, 9, 56, 1.0, rT(2), 24, I, 56, 0.0, _tmp0, 24, offset_rT, offset_I, offset__tmp0, num_elements);
    }
    _tmp1 = d_buffer1;
    {
    unsigned offset_fP = 0;
    unsigned offset__tmp0 = 216;
    unsigned offset__tmp1 = 216;
    
    
    device.gemm(device::ColMajor, device::NoTrans, device::NoTrans, 24, 9, 21, 1.0, fP(0), 24, _tmp0, 24, 0.0, _tmp1, 24, offset_fP, offset__tmp0, offset__tmp1, num_elements);
    }
    _tmp2 = d_buffer0;
    {
    unsigned offset__tmp1 = 216;
    unsigned offset_AminusT = AminusT_offset;
    unsigned offset__tmp2 = 216;
    
    
    device.gemm(device::ColMajor, device::NoTrans, device::NoTrans, 24, 9, 9, 1.0, _tmp1, 24, AminusT, 9, 0.0, _tmp2, 24, offset__tmp1, offset_AminusT, offset__tmp2, num_elements);
    }
    {
    unsigned offset_rDivM = 0;
    unsigned offset__tmp2 = 216;
    unsigned offset_Q = Q_offset;
    
    
    device.gemm(device::ColMajor, device::NoTrans, device::NoTrans, 56, 9, 21, 1.0, rDivM(0), 56, _tmp2, 24, 1.0, Q, 56, offset_rDivM, offset__tmp2, offset_Q, num_elements);
    }
    device.api->popStackMemory();
    device.api->popStackMemory();
  }
  void kernel::neighboringFlux::execute7() {
    assert(AminusT != nullptr);
    assert(I != nullptr);
    assert(Q != nullptr);
    assert(fP(1) != nullptr);
    assert(rDivM(0) != nullptr);
    assert(rT(2) != nullptr);
    assert(kernel::neighboringFlux::num_elements != 0);
    device::DeviceInstance& device = device::DeviceInstance::getInstance();
    float *_tmp0, *_tmp1, *_tmp2;
    float *d_buffer0 = (float*)device.api->getStackMemory(216 * num_elements * sizeof(float));
    float *d_buffer1 = (float*)device.api->getStackMemory(216 * num_elements * sizeof(float));
    _tmp0 = d_buffer0;
    {
    unsigned offset_rT = 0;
    unsigned offset_I = I_offset;
    unsigned offset__tmp0 = 216;
    
    
    device.gemm(device::ColMajor, device::NoTrans, device::NoTrans, 24, 9, 56, 1.0, rT(2), 24, I, 56, 0.0, _tmp0, 24, offset_rT, offset_I, offset__tmp0, num_elements);
    }
    _tmp1 = d_buffer1;
    {
    unsigned offset_fP = 0;
    unsigned offset__tmp0 = 216;
    unsigned offset__tmp1 = 216;
    
    
    device.gemm(device::ColMajor, device::NoTrans, device::NoTrans, 24, 9, 21, 1.0, fP(1), 24, _tmp0, 24, 0.0, _tmp1, 24, offset_fP, offset__tmp0, offset__tmp1, num_elements);
    }
    _tmp2 = d_buffer0;
    {
    unsigned offset__tmp1 = 216;
    unsigned offset_AminusT = AminusT_offset;
    unsigned offset__tmp2 = 216;
    
    
    device.gemm(device::ColMajor, device::NoTrans, device::NoTrans, 24, 9, 9, 1.0, _tmp1, 24, AminusT, 9, 0.0, _tmp2, 24, offset__tmp1, offset_AminusT, offset__tmp2, num_elements);
    }
    {
    unsigned offset_rDivM = 0;
    unsigned offset__tmp2 = 216;
    unsigned offset_Q = Q_offset;
    
    
    device.gemm(device::ColMajor, device::NoTrans, device::NoTrans, 56, 9, 21, 1.0, rDivM(0), 56, _tmp2, 24, 1.0, Q, 56, offset_rDivM, offset__tmp2, offset_Q, num_elements);
    }
    device.api->popStackMemory();
    device.api->popStackMemory();
  }
  void kernel::neighboringFlux::execute8() {
    assert(AminusT != nullptr);
    assert(I != nullptr);
    assert(Q != nullptr);
    assert(fP(2) != nullptr);
    assert(rDivM(0) != nullptr);
    assert(rT(2) != nullptr);
    assert(kernel::neighboringFlux::num_elements != 0);
    device::DeviceInstance& device = device::DeviceInstance::getInstance();
    float *_tmp0, *_tmp1, *_tmp2;
    float *d_buffer0 = (float*)device.api->getStackMemory(216 * num_elements * sizeof(float));
    float *d_buffer1 = (float*)device.api->getStackMemory(216 * num_elements * sizeof(float));
    _tmp0 = d_buffer0;
    {
    unsigned offset_rT = 0;
    unsigned offset_I = I_offset;
    unsigned offset__tmp0 = 216;
    
    
    device.gemm(device::ColMajor, device::NoTrans, device::NoTrans, 24, 9, 56, 1.0, rT(2), 24, I, 56, 0.0, _tmp0, 24, offset_rT, offset_I, offset__tmp0, num_elements);
    }
    _tmp1 = d_buffer1;
    {
    unsigned offset_fP = 0;
    unsigned offset__tmp0 = 216;
    unsigned offset__tmp1 = 216;
    
    
    device.gemm(device::ColMajor, device::NoTrans, device::NoTrans, 24, 9, 21, 1.0, fP(2), 24, _tmp0, 24, 0.0, _tmp1, 24, offset_fP, offset__tmp0, offset__tmp1, num_elements);
    }
    _tmp2 = d_buffer0;
    {
    unsigned offset__tmp1 = 216;
    unsigned offset_AminusT = AminusT_offset;
    unsigned offset__tmp2 = 216;
    
    
    device.gemm(device::ColMajor, device::NoTrans, device::NoTrans, 24, 9, 9, 1.0, _tmp1, 24, AminusT, 9, 0.0, _tmp2, 24, offset__tmp1, offset_AminusT, offset__tmp2, num_elements);
    }
    {
    unsigned offset_rDivM = 0;
    unsigned offset__tmp2 = 216;
    unsigned offset_Q = Q_offset;
    
    
    device.gemm(device::ColMajor, device::NoTrans, device::NoTrans, 56, 9, 21, 1.0, rDivM(0), 56, _tmp2, 24, 1.0, Q, 56, offset_rDivM, offset__tmp2, offset_Q, num_elements);
    }
    device.api->popStackMemory();
    device.api->popStackMemory();
  }
  void kernel::neighboringFlux::execute9() {
    assert(AminusT != nullptr);
    assert(I != nullptr);
    assert(Q != nullptr);
    assert(fP(0) != nullptr);
    assert(rDivM(0) != nullptr);
    assert(rT(3) != nullptr);
    assert(kernel::neighboringFlux::num_elements != 0);
    device::DeviceInstance& device = device::DeviceInstance::getInstance();
    float *_tmp0, *_tmp1, *_tmp2;
    float *d_buffer0 = (float*)device.api->getStackMemory(216 * num_elements * sizeof(float));
    float *d_buffer1 = (float*)device.api->getStackMemory(216 * num_elements * sizeof(float));
    _tmp0 = d_buffer0;
    {
    unsigned offset_rT = 0;
    unsigned offset_I = I_offset;
    unsigned offset__tmp0 = 216;
    
    
    device.gemm(device::ColMajor, device::NoTrans, device::NoTrans, 24, 9, 56, 1.0, rT(3), 24, I, 56, 0.0, _tmp0, 24, offset_rT, offset_I, offset__tmp0, num_elements);
    }
    _tmp1 = d_buffer1;
    {
    unsigned offset_fP = 0;
    unsigned offset__tmp0 = 216;
    unsigned offset__tmp1 = 216;
    
    
    device.gemm(device::ColMajor, device::NoTrans, device::NoTrans, 24, 9, 21, 1.0, fP(0), 24, _tmp0, 24, 0.0, _tmp1, 24, offset_fP, offset__tmp0, offset__tmp1, num_elements);
    }
    _tmp2 = d_buffer0;
    {
    unsigned offset__tmp1 = 216;
    unsigned offset_AminusT = AminusT_offset;
    unsigned offset__tmp2 = 216;
    
    
    device.gemm(device::ColMajor, device::NoTrans, device::NoTrans, 24, 9, 9, 1.0, _tmp1, 24, AminusT, 9, 0.0, _tmp2, 24, offset__tmp1, offset_AminusT, offset__tmp2, num_elements);
    }
    {
    unsigned offset_rDivM = 0;
    unsigned offset__tmp2 = 216;
    unsigned offset_Q = Q_offset;
    
    
    device.gemm(device::ColMajor, device::NoTrans, device::NoTrans, 56, 9, 21, 1.0, rDivM(0), 56, _tmp2, 24, 1.0, Q, 56, offset_rDivM, offset__tmp2, offset_Q, num_elements);
    }
    device.api->popStackMemory();
    device.api->popStackMemory();
  }
  void kernel::neighboringFlux::execute10() {
    assert(AminusT != nullptr);
    assert(I != nullptr);
    assert(Q != nullptr);
    assert(fP(1) != nullptr);
    assert(rDivM(0) != nullptr);
    assert(rT(3) != nullptr);
    assert(kernel::neighboringFlux::num_elements != 0);
    device::DeviceInstance& device = device::DeviceInstance::getInstance();
    float *_tmp0, *_tmp1, *_tmp2;
    float *d_buffer0 = (float*)device.api->getStackMemory(216 * num_elements * sizeof(float));
    float *d_buffer1 = (float*)device.api->getStackMemory(216 * num_elements * sizeof(float));
    _tmp0 = d_buffer0;
    {
    unsigned offset_rT = 0;
    unsigned offset_I = I_offset;
    unsigned offset__tmp0 = 216;
    
    
    device.gemm(device::ColMajor, device::NoTrans, device::NoTrans, 24, 9, 56, 1.0, rT(3), 24, I, 56, 0.0, _tmp0, 24, offset_rT, offset_I, offset__tmp0, num_elements);
    }
    _tmp1 = d_buffer1;
    {
    unsigned offset_fP = 0;
    unsigned offset__tmp0 = 216;
    unsigned offset__tmp1 = 216;
    
    
    device.gemm(device::ColMajor, device::NoTrans, device::NoTrans, 24, 9, 21, 1.0, fP(1), 24, _tmp0, 24, 0.0, _tmp1, 24, offset_fP, offset__tmp0, offset__tmp1, num_elements);
    }
    _tmp2 = d_buffer0;
    {
    unsigned offset__tmp1 = 216;
    unsigned offset_AminusT = AminusT_offset;
    unsigned offset__tmp2 = 216;
    
    
    device.gemm(device::ColMajor, device::NoTrans, device::NoTrans, 24, 9, 9, 1.0, _tmp1, 24, AminusT, 9, 0.0, _tmp2, 24, offset__tmp1, offset_AminusT, offset__tmp2, num_elements);
    }
    {
    unsigned offset_rDivM = 0;
    unsigned offset__tmp2 = 216;
    unsigned offset_Q = Q_offset;
    
    
    device.gemm(device::ColMajor, device::NoTrans, device::NoTrans, 56, 9, 21, 1.0, rDivM(0), 56, _tmp2, 24, 1.0, Q, 56, offset_rDivM, offset__tmp2, offset_Q, num_elements);
    }
    device.api->popStackMemory();
    device.api->popStackMemory();
  }
  void kernel::neighboringFlux::execute11() {
    assert(AminusT != nullptr);
    assert(I != nullptr);
    assert(Q != nullptr);
    assert(fP(2) != nullptr);
    assert(rDivM(0) != nullptr);
    assert(rT(3) != nullptr);
    assert(kernel::neighboringFlux::num_elements != 0);
    device::DeviceInstance& device = device::DeviceInstance::getInstance();
    float *_tmp0, *_tmp1, *_tmp2;
    float *d_buffer0 = (float*)device.api->getStackMemory(216 * num_elements * sizeof(float));
    float *d_buffer1 = (float*)device.api->getStackMemory(216 * num_elements * sizeof(float));
    _tmp0 = d_buffer0;
    {
    unsigned offset_rT = 0;
    unsigned offset_I = I_offset;
    unsigned offset__tmp0 = 216;
    
    
    device.gemm(device::ColMajor, device::NoTrans, device::NoTrans, 24, 9, 56, 1.0, rT(3), 24, I, 56, 0.0, _tmp0, 24, offset_rT, offset_I, offset__tmp0, num_elements);
    }
    _tmp1 = d_buffer1;
    {
    unsigned offset_fP = 0;
    unsigned offset__tmp0 = 216;
    unsigned offset__tmp1 = 216;
    
    
    device.gemm(device::ColMajor, device::NoTrans, device::NoTrans, 24, 9, 21, 1.0, fP(2), 24, _tmp0, 24, 0.0, _tmp1, 24, offset_fP, offset__tmp0, offset__tmp1, num_elements);
    }
    _tmp2 = d_buffer0;
    {
    unsigned offset__tmp1 = 216;
    unsigned offset_AminusT = AminusT_offset;
    unsigned offset__tmp2 = 216;
    
    
    device.gemm(device::ColMajor, device::NoTrans, device::NoTrans, 24, 9, 9, 1.0, _tmp1, 24, AminusT, 9, 0.0, _tmp2, 24, offset__tmp1, offset_AminusT, offset__tmp2, num_elements);
    }
    {
    unsigned offset_rDivM = 0;
    unsigned offset__tmp2 = 216;
    unsigned offset_Q = Q_offset;
    
    
    device.gemm(device::ColMajor, device::NoTrans, device::NoTrans, 56, 9, 21, 1.0, rDivM(0), 56, _tmp2, 24, 1.0, Q, 56, offset_rDivM, offset__tmp2, offset_Q, num_elements);
    }
    device.api->popStackMemory();
    device.api->popStackMemory();
  }
  void kernel::neighboringFlux::execute12() {
    assert(AminusT != nullptr);
    assert(I != nullptr);
    assert(Q != nullptr);
    assert(fP(0) != nullptr);
    assert(rDivM(1) != nullptr);
    assert(rT(0) != nullptr);
    assert(kernel::neighboringFlux::num_elements != 0);
    device::DeviceInstance& device = device::DeviceInstance::getInstance();
    float *_tmp0, *_tmp1, *_tmp2;
    float *d_buffer0 = (float*)device.api->getStackMemory(216 * num_elements * sizeof(float));
    float *d_buffer1 = (float*)device.api->getStackMemory(216 * num_elements * sizeof(float));
    _tmp0 = d_buffer0;
    {
    unsigned offset_rT = 0;
    unsigned offset_I = I_offset;
    unsigned offset__tmp0 = 216;
    
    
    device.gemm(device::ColMajor, device::NoTrans, device::NoTrans, 24, 9, 56, 1.0, rT(0), 24, I, 56, 0.0, _tmp0, 24, offset_rT, offset_I, offset__tmp0, num_elements);
    }
    _tmp1 = d_buffer1;
    {
    unsigned offset_fP = 0;
    unsigned offset__tmp0 = 216;
    unsigned offset__tmp1 = 216;
    
    
    device.gemm(device::ColMajor, device::NoTrans, device::NoTrans, 24, 9, 21, 1.0, fP(0), 24, _tmp0, 24, 0.0, _tmp1, 24, offset_fP, offset__tmp0, offset__tmp1, num_elements);
    }
    _tmp2 = d_buffer0;
    {
    unsigned offset__tmp1 = 216;
    unsigned offset_AminusT = AminusT_offset;
    unsigned offset__tmp2 = 216;
    
    
    device.gemm(device::ColMajor, device::NoTrans, device::NoTrans, 24, 9, 9, 1.0, _tmp1, 24, AminusT, 9, 0.0, _tmp2, 24, offset__tmp1, offset_AminusT, offset__tmp2, num_elements);
    }
    {
    unsigned offset_rDivM = 0;
    unsigned offset__tmp2 = 216;
    unsigned offset_Q = Q_offset;
    
    
    device.gemm(device::ColMajor, device::NoTrans, device::NoTrans, 56, 9, 21, 1.0, rDivM(1), 56, _tmp2, 24, 1.0, Q, 56, offset_rDivM, offset__tmp2, offset_Q, num_elements);
    }
    device.api->popStackMemory();
    device.api->popStackMemory();
  }
  void kernel::neighboringFlux::execute13() {
    assert(AminusT != nullptr);
    assert(I != nullptr);
    assert(Q != nullptr);
    assert(fP(1) != nullptr);
    assert(rDivM(1) != nullptr);
    assert(rT(0) != nullptr);
    assert(kernel::neighboringFlux::num_elements != 0);
    device::DeviceInstance& device = device::DeviceInstance::getInstance();
    float *_tmp0, *_tmp1, *_tmp2;
    float *d_buffer0 = (float*)device.api->getStackMemory(216 * num_elements * sizeof(float));
    float *d_buffer1 = (float*)device.api->getStackMemory(216 * num_elements * sizeof(float));
    _tmp0 = d_buffer0;
    {
    unsigned offset_rT = 0;
    unsigned offset_I = I_offset;
    unsigned offset__tmp0 = 216;
    
    
    device.gemm(device::ColMajor, device::NoTrans, device::NoTrans, 24, 9, 56, 1.0, rT(0), 24, I, 56, 0.0, _tmp0, 24, offset_rT, offset_I, offset__tmp0, num_elements);
    }
    _tmp1 = d_buffer1;
    {
    unsigned offset_fP = 0;
    unsigned offset__tmp0 = 216;
    unsigned offset__tmp1 = 216;
    
    
    device.gemm(device::ColMajor, device::NoTrans, device::NoTrans, 24, 9, 21, 1.0, fP(1), 24, _tmp0, 24, 0.0, _tmp1, 24, offset_fP, offset__tmp0, offset__tmp1, num_elements);
    }
    _tmp2 = d_buffer0;
    {
    unsigned offset__tmp1 = 216;
    unsigned offset_AminusT = AminusT_offset;
    unsigned offset__tmp2 = 216;
    
    
    device.gemm(device::ColMajor, device::NoTrans, device::NoTrans, 24, 9, 9, 1.0, _tmp1, 24, AminusT, 9, 0.0, _tmp2, 24, offset__tmp1, offset_AminusT, offset__tmp2, num_elements);
    }
    {
    unsigned offset_rDivM = 0;
    unsigned offset__tmp2 = 216;
    unsigned offset_Q = Q_offset;
    
    
    device.gemm(device::ColMajor, device::NoTrans, device::NoTrans, 56, 9, 21, 1.0, rDivM(1), 56, _tmp2, 24, 1.0, Q, 56, offset_rDivM, offset__tmp2, offset_Q, num_elements);
    }
    device.api->popStackMemory();
    device.api->popStackMemory();
  }
  void kernel::neighboringFlux::execute14() {
    assert(AminusT != nullptr);
    assert(I != nullptr);
    assert(Q != nullptr);
    assert(fP(2) != nullptr);
    assert(rDivM(1) != nullptr);
    assert(rT(0) != nullptr);
    assert(kernel::neighboringFlux::num_elements != 0);
    device::DeviceInstance& device = device::DeviceInstance::getInstance();
    float *_tmp0, *_tmp1, *_tmp2;
    float *d_buffer0 = (float*)device.api->getStackMemory(216 * num_elements * sizeof(float));
    float *d_buffer1 = (float*)device.api->getStackMemory(216 * num_elements * sizeof(float));
    _tmp0 = d_buffer0;
    {
    unsigned offset_rT = 0;
    unsigned offset_I = I_offset;
    unsigned offset__tmp0 = 216;
    
    
    device.gemm(device::ColMajor, device::NoTrans, device::NoTrans, 24, 9, 56, 1.0, rT(0), 24, I, 56, 0.0, _tmp0, 24, offset_rT, offset_I, offset__tmp0, num_elements);
    }
    _tmp1 = d_buffer1;
    {
    unsigned offset_fP = 0;
    unsigned offset__tmp0 = 216;
    unsigned offset__tmp1 = 216;
    
    
    device.gemm(device::ColMajor, device::NoTrans, device::NoTrans, 24, 9, 21, 1.0, fP(2), 24, _tmp0, 24, 0.0, _tmp1, 24, offset_fP, offset__tmp0, offset__tmp1, num_elements);
    }
    _tmp2 = d_buffer0;
    {
    unsigned offset__tmp1 = 216;
    unsigned offset_AminusT = AminusT_offset;
    unsigned offset__tmp2 = 216;
    
    
    device.gemm(device::ColMajor, device::NoTrans, device::NoTrans, 24, 9, 9, 1.0, _tmp1, 24, AminusT, 9, 0.0, _tmp2, 24, offset__tmp1, offset_AminusT, offset__tmp2, num_elements);
    }
    {
    unsigned offset_rDivM = 0;
    unsigned offset__tmp2 = 216;
    unsigned offset_Q = Q_offset;
    
    
    device.gemm(device::ColMajor, device::NoTrans, device::NoTrans, 56, 9, 21, 1.0, rDivM(1), 56, _tmp2, 24, 1.0, Q, 56, offset_rDivM, offset__tmp2, offset_Q, num_elements);
    }
    device.api->popStackMemory();
    device.api->popStackMemory();
  }
  void kernel::neighboringFlux::execute15() {
    assert(AminusT != nullptr);
    assert(I != nullptr);
    assert(Q != nullptr);
    assert(fP(0) != nullptr);
    assert(rDivM(1) != nullptr);
    assert(rT(1) != nullptr);
    assert(kernel::neighboringFlux::num_elements != 0);
    device::DeviceInstance& device = device::DeviceInstance::getInstance();
    float *_tmp0, *_tmp1, *_tmp2;
    float *d_buffer0 = (float*)device.api->getStackMemory(216 * num_elements * sizeof(float));
    float *d_buffer1 = (float*)device.api->getStackMemory(216 * num_elements * sizeof(float));
    _tmp0 = d_buffer0;
    {
    unsigned offset_rT = 0;
    unsigned offset_I = I_offset;
    unsigned offset__tmp0 = 216;
    
    
    device.gemm(device::ColMajor, device::NoTrans, device::NoTrans, 24, 9, 56, 1.0, rT(1), 24, I, 56, 0.0, _tmp0, 24, offset_rT, offset_I, offset__tmp0, num_elements);
    }
    _tmp1 = d_buffer1;
    {
    unsigned offset_fP = 0;
    unsigned offset__tmp0 = 216;
    unsigned offset__tmp1 = 216;
    
    
    device.gemm(device::ColMajor, device::NoTrans, device::NoTrans, 24, 9, 21, 1.0, fP(0), 24, _tmp0, 24, 0.0, _tmp1, 24, offset_fP, offset__tmp0, offset__tmp1, num_elements);
    }
    _tmp2 = d_buffer0;
    {
    unsigned offset__tmp1 = 216;
    unsigned offset_AminusT = AminusT_offset;
    unsigned offset__tmp2 = 216;
    
    
    device.gemm(device::ColMajor, device::NoTrans, device::NoTrans, 24, 9, 9, 1.0, _tmp1, 24, AminusT, 9, 0.0, _tmp2, 24, offset__tmp1, offset_AminusT, offset__tmp2, num_elements);
    }
    {
    unsigned offset_rDivM = 0;
    unsigned offset__tmp2 = 216;
    unsigned offset_Q = Q_offset;
    
    
    device.gemm(device::ColMajor, device::NoTrans, device::NoTrans, 56, 9, 21, 1.0, rDivM(1), 56, _tmp2, 24, 1.0, Q, 56, offset_rDivM, offset__tmp2, offset_Q, num_elements);
    }
    device.api->popStackMemory();
    device.api->popStackMemory();
  }
  void kernel::neighboringFlux::execute16() {
    assert(AminusT != nullptr);
    assert(I != nullptr);
    assert(Q != nullptr);
    assert(fP(1) != nullptr);
    assert(rDivM(1) != nullptr);
    assert(rT(1) != nullptr);
    assert(kernel::neighboringFlux::num_elements != 0);
    device::DeviceInstance& device = device::DeviceInstance::getInstance();
    float *_tmp0, *_tmp1, *_tmp2;
    float *d_buffer0 = (float*)device.api->getStackMemory(216 * num_elements * sizeof(float));
    float *d_buffer1 = (float*)device.api->getStackMemory(216 * num_elements * sizeof(float));
    _tmp0 = d_buffer0;
    {
    unsigned offset_rT = 0;
    unsigned offset_I = I_offset;
    unsigned offset__tmp0 = 216;
    
    
    device.gemm(device::ColMajor, device::NoTrans, device::NoTrans, 24, 9, 56, 1.0, rT(1), 24, I, 56, 0.0, _tmp0, 24, offset_rT, offset_I, offset__tmp0, num_elements);
    }
    _tmp1 = d_buffer1;
    {
    unsigned offset_fP = 0;
    unsigned offset__tmp0 = 216;
    unsigned offset__tmp1 = 216;
    
    
    device.gemm(device::ColMajor, device::NoTrans, device::NoTrans, 24, 9, 21, 1.0, fP(1), 24, _tmp0, 24, 0.0, _tmp1, 24, offset_fP, offset__tmp0, offset__tmp1, num_elements);
    }
    _tmp2 = d_buffer0;
    {
    unsigned offset__tmp1 = 216;
    unsigned offset_AminusT = AminusT_offset;
    unsigned offset__tmp2 = 216;
    
    
    device.gemm(device::ColMajor, device::NoTrans, device::NoTrans, 24, 9, 9, 1.0, _tmp1, 24, AminusT, 9, 0.0, _tmp2, 24, offset__tmp1, offset_AminusT, offset__tmp2, num_elements);
    }
    {
    unsigned offset_rDivM = 0;
    unsigned offset__tmp2 = 216;
    unsigned offset_Q = Q_offset;
    
    
    device.gemm(device::ColMajor, device::NoTrans, device::NoTrans, 56, 9, 21, 1.0, rDivM(1), 56, _tmp2, 24, 1.0, Q, 56, offset_rDivM, offset__tmp2, offset_Q, num_elements);
    }
    device.api->popStackMemory();
    device.api->popStackMemory();
  }
  void kernel::neighboringFlux::execute17() {
    assert(AminusT != nullptr);
    assert(I != nullptr);
    assert(Q != nullptr);
    assert(fP(2) != nullptr);
    assert(rDivM(1) != nullptr);
    assert(rT(1) != nullptr);
    assert(kernel::neighboringFlux::num_elements != 0);
    device::DeviceInstance& device = device::DeviceInstance::getInstance();
    float *_tmp0, *_tmp1, *_tmp2;
    float *d_buffer0 = (float*)device.api->getStackMemory(216 * num_elements * sizeof(float));
    float *d_buffer1 = (float*)device.api->getStackMemory(216 * num_elements * sizeof(float));
    _tmp0 = d_buffer0;
    {
    unsigned offset_rT = 0;
    unsigned offset_I = I_offset;
    unsigned offset__tmp0 = 216;
    
    
    device.gemm(device::ColMajor, device::NoTrans, device::NoTrans, 24, 9, 56, 1.0, rT(1), 24, I, 56, 0.0, _tmp0, 24, offset_rT, offset_I, offset__tmp0, num_elements);
    }
    _tmp1 = d_buffer1;
    {
    unsigned offset_fP = 0;
    unsigned offset__tmp0 = 216;
    unsigned offset__tmp1 = 216;
    
    
    device.gemm(device::ColMajor, device::NoTrans, device::NoTrans, 24, 9, 21, 1.0, fP(2), 24, _tmp0, 24, 0.0, _tmp1, 24, offset_fP, offset__tmp0, offset__tmp1, num_elements);
    }
    _tmp2 = d_buffer0;
    {
    unsigned offset__tmp1 = 216;
    unsigned offset_AminusT = AminusT_offset;
    unsigned offset__tmp2 = 216;
    
    
    device.gemm(device::ColMajor, device::NoTrans, device::NoTrans, 24, 9, 9, 1.0, _tmp1, 24, AminusT, 9, 0.0, _tmp2, 24, offset__tmp1, offset_AminusT, offset__tmp2, num_elements);
    }
    {
    unsigned offset_rDivM = 0;
    unsigned offset__tmp2 = 216;
    unsigned offset_Q = Q_offset;
    
    
    device.gemm(device::ColMajor, device::NoTrans, device::NoTrans, 56, 9, 21, 1.0, rDivM(1), 56, _tmp2, 24, 1.0, Q, 56, offset_rDivM, offset__tmp2, offset_Q, num_elements);
    }
    device.api->popStackMemory();
    device.api->popStackMemory();
  }
  void kernel::neighboringFlux::execute18() {
    assert(AminusT != nullptr);
    assert(I != nullptr);
    assert(Q != nullptr);
    assert(fP(0) != nullptr);
    assert(rDivM(1) != nullptr);
    assert(rT(2) != nullptr);
    assert(kernel::neighboringFlux::num_elements != 0);
    device::DeviceInstance& device = device::DeviceInstance::getInstance();
    float *_tmp0, *_tmp1, *_tmp2;
    float *d_buffer0 = (float*)device.api->getStackMemory(216 * num_elements * sizeof(float));
    float *d_buffer1 = (float*)device.api->getStackMemory(216 * num_elements * sizeof(float));
    _tmp0 = d_buffer0;
    {
    unsigned offset_rT = 0;
    unsigned offset_I = I_offset;
    unsigned offset__tmp0 = 216;
    
    
    device.gemm(device::ColMajor, device::NoTrans, device::NoTrans, 24, 9, 56, 1.0, rT(2), 24, I, 56, 0.0, _tmp0, 24, offset_rT, offset_I, offset__tmp0, num_elements);
    }
    _tmp1 = d_buffer1;
    {
    unsigned offset_fP = 0;
    unsigned offset__tmp0 = 216;
    unsigned offset__tmp1 = 216;
    
    
    device.gemm(device::ColMajor, device::NoTrans, device::NoTrans, 24, 9, 21, 1.0, fP(0), 24, _tmp0, 24, 0.0, _tmp1, 24, offset_fP, offset__tmp0, offset__tmp1, num_elements);
    }
    _tmp2 = d_buffer0;
    {
    unsigned offset__tmp1 = 216;
    unsigned offset_AminusT = AminusT_offset;
    unsigned offset__tmp2 = 216;
    
    
    device.gemm(device::ColMajor, device::NoTrans, device::NoTrans, 24, 9, 9, 1.0, _tmp1, 24, AminusT, 9, 0.0, _tmp2, 24, offset__tmp1, offset_AminusT, offset__tmp2, num_elements);
    }
    {
    unsigned offset_rDivM = 0;
    unsigned offset__tmp2 = 216;
    unsigned offset_Q = Q_offset;
    
    
    device.gemm(device::ColMajor, device::NoTrans, device::NoTrans, 56, 9, 21, 1.0, rDivM(1), 56, _tmp2, 24, 1.0, Q, 56, offset_rDivM, offset__tmp2, offset_Q, num_elements);
    }
    device.api->popStackMemory();
    device.api->popStackMemory();
  }
  void kernel::neighboringFlux::execute19() {
    assert(AminusT != nullptr);
    assert(I != nullptr);
    assert(Q != nullptr);
    assert(fP(1) != nullptr);
    assert(rDivM(1) != nullptr);
    assert(rT(2) != nullptr);
    assert(kernel::neighboringFlux::num_elements != 0);
    device::DeviceInstance& device = device::DeviceInstance::getInstance();
    float *_tmp0, *_tmp1, *_tmp2;
    float *d_buffer0 = (float*)device.api->getStackMemory(216 * num_elements * sizeof(float));
    float *d_buffer1 = (float*)device.api->getStackMemory(216 * num_elements * sizeof(float));
    _tmp0 = d_buffer0;
    {
    unsigned offset_rT = 0;
    unsigned offset_I = I_offset;
    unsigned offset__tmp0 = 216;
    
    
    device.gemm(device::ColMajor, device::NoTrans, device::NoTrans, 24, 9, 56, 1.0, rT(2), 24, I, 56, 0.0, _tmp0, 24, offset_rT, offset_I, offset__tmp0, num_elements);
    }
    _tmp1 = d_buffer1;
    {
    unsigned offset_fP = 0;
    unsigned offset__tmp0 = 216;
    unsigned offset__tmp1 = 216;
    
    
    device.gemm(device::ColMajor, device::NoTrans, device::NoTrans, 24, 9, 21, 1.0, fP(1), 24, _tmp0, 24, 0.0, _tmp1, 24, offset_fP, offset__tmp0, offset__tmp1, num_elements);
    }
    _tmp2 = d_buffer0;
    {
    unsigned offset__tmp1 = 216;
    unsigned offset_AminusT = AminusT_offset;
    unsigned offset__tmp2 = 216;
    
    
    device.gemm(device::ColMajor, device::NoTrans, device::NoTrans, 24, 9, 9, 1.0, _tmp1, 24, AminusT, 9, 0.0, _tmp2, 24, offset__tmp1, offset_AminusT, offset__tmp2, num_elements);
    }
    {
    unsigned offset_rDivM = 0;
    unsigned offset__tmp2 = 216;
    unsigned offset_Q = Q_offset;
    
    
    device.gemm(device::ColMajor, device::NoTrans, device::NoTrans, 56, 9, 21, 1.0, rDivM(1), 56, _tmp2, 24, 1.0, Q, 56, offset_rDivM, offset__tmp2, offset_Q, num_elements);
    }
    device.api->popStackMemory();
    device.api->popStackMemory();
  }
  void kernel::neighboringFlux::execute20() {
    assert(AminusT != nullptr);
    assert(I != nullptr);
    assert(Q != nullptr);
    assert(fP(2) != nullptr);
    assert(rDivM(1) != nullptr);
    assert(rT(2) != nullptr);
    assert(kernel::neighboringFlux::num_elements != 0);
    device::DeviceInstance& device = device::DeviceInstance::getInstance();
    float *_tmp0, *_tmp1, *_tmp2;
    float *d_buffer0 = (float*)device.api->getStackMemory(216 * num_elements * sizeof(float));
    float *d_buffer1 = (float*)device.api->getStackMemory(216 * num_elements * sizeof(float));
    _tmp0 = d_buffer0;
    {
    unsigned offset_rT = 0;
    unsigned offset_I = I_offset;
    unsigned offset__tmp0 = 216;
    
    
    device.gemm(device::ColMajor, device::NoTrans, device::NoTrans, 24, 9, 56, 1.0, rT(2), 24, I, 56, 0.0, _tmp0, 24, offset_rT, offset_I, offset__tmp0, num_elements);
    }
    _tmp1 = d_buffer1;
    {
    unsigned offset_fP = 0;
    unsigned offset__tmp0 = 216;
    unsigned offset__tmp1 = 216;
    
    
    device.gemm(device::ColMajor, device::NoTrans, device::NoTrans, 24, 9, 21, 1.0, fP(2), 24, _tmp0, 24, 0.0, _tmp1, 24, offset_fP, offset__tmp0, offset__tmp1, num_elements);
    }
    _tmp2 = d_buffer0;
    {
    unsigned offset__tmp1 = 216;
    unsigned offset_AminusT = AminusT_offset;
    unsigned offset__tmp2 = 216;
    
    
    device.gemm(device::ColMajor, device::NoTrans, device::NoTrans, 24, 9, 9, 1.0, _tmp1, 24, AminusT, 9, 0.0, _tmp2, 24, offset__tmp1, offset_AminusT, offset__tmp2, num_elements);
    }
    {
    unsigned offset_rDivM = 0;
    unsigned offset__tmp2 = 216;
    unsigned offset_Q = Q_offset;
    
    
    device.gemm(device::ColMajor, device::NoTrans, device::NoTrans, 56, 9, 21, 1.0, rDivM(1), 56, _tmp2, 24, 1.0, Q, 56, offset_rDivM, offset__tmp2, offset_Q, num_elements);
    }
    device.api->popStackMemory();
    device.api->popStackMemory();
  }
  void kernel::neighboringFlux::execute21() {
    assert(AminusT != nullptr);
    assert(I != nullptr);
    assert(Q != nullptr);
    assert(fP(0) != nullptr);
    assert(rDivM(1) != nullptr);
    assert(rT(3) != nullptr);
    assert(kernel::neighboringFlux::num_elements != 0);
    device::DeviceInstance& device = device::DeviceInstance::getInstance();
    float *_tmp0, *_tmp1, *_tmp2;
    float *d_buffer0 = (float*)device.api->getStackMemory(216 * num_elements * sizeof(float));
    float *d_buffer1 = (float*)device.api->getStackMemory(216 * num_elements * sizeof(float));
    _tmp0 = d_buffer0;
    {
    unsigned offset_rT = 0;
    unsigned offset_I = I_offset;
    unsigned offset__tmp0 = 216;
    
    
    device.gemm(device::ColMajor, device::NoTrans, device::NoTrans, 24, 9, 56, 1.0, rT(3), 24, I, 56, 0.0, _tmp0, 24, offset_rT, offset_I, offset__tmp0, num_elements);
    }
    _tmp1 = d_buffer1;
    {
    unsigned offset_fP = 0;
    unsigned offset__tmp0 = 216;
    unsigned offset__tmp1 = 216;
    
    
    device.gemm(device::ColMajor, device::NoTrans, device::NoTrans, 24, 9, 21, 1.0, fP(0), 24, _tmp0, 24, 0.0, _tmp1, 24, offset_fP, offset__tmp0, offset__tmp1, num_elements);
    }
    _tmp2 = d_buffer0;
    {
    unsigned offset__tmp1 = 216;
    unsigned offset_AminusT = AminusT_offset;
    unsigned offset__tmp2 = 216;
    
    
    device.gemm(device::ColMajor, device::NoTrans, device::NoTrans, 24, 9, 9, 1.0, _tmp1, 24, AminusT, 9, 0.0, _tmp2, 24, offset__tmp1, offset_AminusT, offset__tmp2, num_elements);
    }
    {
    unsigned offset_rDivM = 0;
    unsigned offset__tmp2 = 216;
    unsigned offset_Q = Q_offset;
    
    
    device.gemm(device::ColMajor, device::NoTrans, device::NoTrans, 56, 9, 21, 1.0, rDivM(1), 56, _tmp2, 24, 1.0, Q, 56, offset_rDivM, offset__tmp2, offset_Q, num_elements);
    }
    device.api->popStackMemory();
    device.api->popStackMemory();
  }
  void kernel::neighboringFlux::execute22() {
    assert(AminusT != nullptr);
    assert(I != nullptr);
    assert(Q != nullptr);
    assert(fP(1) != nullptr);
    assert(rDivM(1) != nullptr);
    assert(rT(3) != nullptr);
    assert(kernel::neighboringFlux::num_elements != 0);
    device::DeviceInstance& device = device::DeviceInstance::getInstance();
    float *_tmp0, *_tmp1, *_tmp2;
    float *d_buffer0 = (float*)device.api->getStackMemory(216 * num_elements * sizeof(float));
    float *d_buffer1 = (float*)device.api->getStackMemory(216 * num_elements * sizeof(float));
    _tmp0 = d_buffer0;
    {
    unsigned offset_rT = 0;
    unsigned offset_I = I_offset;
    unsigned offset__tmp0 = 216;
    
    
    device.gemm(device::ColMajor, device::NoTrans, device::NoTrans, 24, 9, 56, 1.0, rT(3), 24, I, 56, 0.0, _tmp0, 24, offset_rT, offset_I, offset__tmp0, num_elements);
    }
    _tmp1 = d_buffer1;
    {
    unsigned offset_fP = 0;
    unsigned offset__tmp0 = 216;
    unsigned offset__tmp1 = 216;
    
    
    device.gemm(device::ColMajor, device::NoTrans, device::NoTrans, 24, 9, 21, 1.0, fP(1), 24, _tmp0, 24, 0.0, _tmp1, 24, offset_fP, offset__tmp0, offset__tmp1, num_elements);
    }
    _tmp2 = d_buffer0;
    {
    unsigned offset__tmp1 = 216;
    unsigned offset_AminusT = AminusT_offset;
    unsigned offset__tmp2 = 216;
    
    
    device.gemm(device::ColMajor, device::NoTrans, device::NoTrans, 24, 9, 9, 1.0, _tmp1, 24, AminusT, 9, 0.0, _tmp2, 24, offset__tmp1, offset_AminusT, offset__tmp2, num_elements);
    }
    {
    unsigned offset_rDivM = 0;
    unsigned offset__tmp2 = 216;
    unsigned offset_Q = Q_offset;
    
    
    device.gemm(device::ColMajor, device::NoTrans, device::NoTrans, 56, 9, 21, 1.0, rDivM(1), 56, _tmp2, 24, 1.0, Q, 56, offset_rDivM, offset__tmp2, offset_Q, num_elements);
    }
    device.api->popStackMemory();
    device.api->popStackMemory();
  }
  void kernel::neighboringFlux::execute23() {
    assert(AminusT != nullptr);
    assert(I != nullptr);
    assert(Q != nullptr);
    assert(fP(2) != nullptr);
    assert(rDivM(1) != nullptr);
    assert(rT(3) != nullptr);
    assert(kernel::neighboringFlux::num_elements != 0);
    device::DeviceInstance& device = device::DeviceInstance::getInstance();
    float *_tmp0, *_tmp1, *_tmp2;
    float *d_buffer0 = (float*)device.api->getStackMemory(216 * num_elements * sizeof(float));
    float *d_buffer1 = (float*)device.api->getStackMemory(216 * num_elements * sizeof(float));
    _tmp0 = d_buffer0;
    {
    unsigned offset_rT = 0;
    unsigned offset_I = I_offset;
    unsigned offset__tmp0 = 216;
    
    
    device.gemm(device::ColMajor, device::NoTrans, device::NoTrans, 24, 9, 56, 1.0, rT(3), 24, I, 56, 0.0, _tmp0, 24, offset_rT, offset_I, offset__tmp0, num_elements);
    }
    _tmp1 = d_buffer1;
    {
    unsigned offset_fP = 0;
    unsigned offset__tmp0 = 216;
    unsigned offset__tmp1 = 216;
    
    
    device.gemm(device::ColMajor, device::NoTrans, device::NoTrans, 24, 9, 21, 1.0, fP(2), 24, _tmp0, 24, 0.0, _tmp1, 24, offset_fP, offset__tmp0, offset__tmp1, num_elements);
    }
    _tmp2 = d_buffer0;
    {
    unsigned offset__tmp1 = 216;
    unsigned offset_AminusT = AminusT_offset;
    unsigned offset__tmp2 = 216;
    
    
    device.gemm(device::ColMajor, device::NoTrans, device::NoTrans, 24, 9, 9, 1.0, _tmp1, 24, AminusT, 9, 0.0, _tmp2, 24, offset__tmp1, offset_AminusT, offset__tmp2, num_elements);
    }
    {
    unsigned offset_rDivM = 0;
    unsigned offset__tmp2 = 216;
    unsigned offset_Q = Q_offset;
    
    
    device.gemm(device::ColMajor, device::NoTrans, device::NoTrans, 56, 9, 21, 1.0, rDivM(1), 56, _tmp2, 24, 1.0, Q, 56, offset_rDivM, offset__tmp2, offset_Q, num_elements);
    }
    device.api->popStackMemory();
    device.api->popStackMemory();
  }
  void kernel::neighboringFlux::execute24() {
    assert(AminusT != nullptr);
    assert(I != nullptr);
    assert(Q != nullptr);
    assert(fP(0) != nullptr);
    assert(rDivM(2) != nullptr);
    assert(rT(0) != nullptr);
    assert(kernel::neighboringFlux::num_elements != 0);
    device::DeviceInstance& device = device::DeviceInstance::getInstance();
    float *_tmp0, *_tmp1, *_tmp2;
    float *d_buffer0 = (float*)device.api->getStackMemory(216 * num_elements * sizeof(float));
    float *d_buffer1 = (float*)device.api->getStackMemory(216 * num_elements * sizeof(float));
    _tmp0 = d_buffer0;
    {
    unsigned offset_rT = 0;
    unsigned offset_I = I_offset;
    unsigned offset__tmp0 = 216;
    
    
    device.gemm(device::ColMajor, device::NoTrans, device::NoTrans, 24, 9, 56, 1.0, rT(0), 24, I, 56, 0.0, _tmp0, 24, offset_rT, offset_I, offset__tmp0, num_elements);
    }
    _tmp1 = d_buffer1;
    {
    unsigned offset_fP = 0;
    unsigned offset__tmp0 = 216;
    unsigned offset__tmp1 = 216;
    
    
    device.gemm(device::ColMajor, device::NoTrans, device::NoTrans, 24, 9, 21, 1.0, fP(0), 24, _tmp0, 24, 0.0, _tmp1, 24, offset_fP, offset__tmp0, offset__tmp1, num_elements);
    }
    _tmp2 = d_buffer0;
    {
    unsigned offset__tmp1 = 216;
    unsigned offset_AminusT = AminusT_offset;
    unsigned offset__tmp2 = 216;
    
    
    device.gemm(device::ColMajor, device::NoTrans, device::NoTrans, 24, 9, 9, 1.0, _tmp1, 24, AminusT, 9, 0.0, _tmp2, 24, offset__tmp1, offset_AminusT, offset__tmp2, num_elements);
    }
    {
    unsigned offset_rDivM = 0;
    unsigned offset__tmp2 = 216;
    unsigned offset_Q = Q_offset;
    
    
    device.gemm(device::ColMajor, device::NoTrans, device::NoTrans, 56, 9, 21, 1.0, rDivM(2), 56, _tmp2, 24, 1.0, Q, 56, offset_rDivM, offset__tmp2, offset_Q, num_elements);
    }
    device.api->popStackMemory();
    device.api->popStackMemory();
  }
  void kernel::neighboringFlux::execute25() {
    assert(AminusT != nullptr);
    assert(I != nullptr);
    assert(Q != nullptr);
    assert(fP(1) != nullptr);
    assert(rDivM(2) != nullptr);
    assert(rT(0) != nullptr);
    assert(kernel::neighboringFlux::num_elements != 0);
    device::DeviceInstance& device = device::DeviceInstance::getInstance();
    float *_tmp0, *_tmp1, *_tmp2;
    float *d_buffer0 = (float*)device.api->getStackMemory(216 * num_elements * sizeof(float));
    float *d_buffer1 = (float*)device.api->getStackMemory(216 * num_elements * sizeof(float));
    _tmp0 = d_buffer0;
    {
    unsigned offset_rT = 0;
    unsigned offset_I = I_offset;
    unsigned offset__tmp0 = 216;
    
    
    device.gemm(device::ColMajor, device::NoTrans, device::NoTrans, 24, 9, 56, 1.0, rT(0), 24, I, 56, 0.0, _tmp0, 24, offset_rT, offset_I, offset__tmp0, num_elements);
    }
    _tmp1 = d_buffer1;
    {
    unsigned offset_fP = 0;
    unsigned offset__tmp0 = 216;
    unsigned offset__tmp1 = 216;
    
    
    device.gemm(device::ColMajor, device::NoTrans, device::NoTrans, 24, 9, 21, 1.0, fP(1), 24, _tmp0, 24, 0.0, _tmp1, 24, offset_fP, offset__tmp0, offset__tmp1, num_elements);
    }
    _tmp2 = d_buffer0;
    {
    unsigned offset__tmp1 = 216;
    unsigned offset_AminusT = AminusT_offset;
    unsigned offset__tmp2 = 216;
    
    
    device.gemm(device::ColMajor, device::NoTrans, device::NoTrans, 24, 9, 9, 1.0, _tmp1, 24, AminusT, 9, 0.0, _tmp2, 24, offset__tmp1, offset_AminusT, offset__tmp2, num_elements);
    }
    {
    unsigned offset_rDivM = 0;
    unsigned offset__tmp2 = 216;
    unsigned offset_Q = Q_offset;
    
    
    device.gemm(device::ColMajor, device::NoTrans, device::NoTrans, 56, 9, 21, 1.0, rDivM(2), 56, _tmp2, 24, 1.0, Q, 56, offset_rDivM, offset__tmp2, offset_Q, num_elements);
    }
    device.api->popStackMemory();
    device.api->popStackMemory();
  }
  void kernel::neighboringFlux::execute26() {
    assert(AminusT != nullptr);
    assert(I != nullptr);
    assert(Q != nullptr);
    assert(fP(2) != nullptr);
    assert(rDivM(2) != nullptr);
    assert(rT(0) != nullptr);
    assert(kernel::neighboringFlux::num_elements != 0);
    device::DeviceInstance& device = device::DeviceInstance::getInstance();
    float *_tmp0, *_tmp1, *_tmp2;
    float *d_buffer0 = (float*)device.api->getStackMemory(216 * num_elements * sizeof(float));
    float *d_buffer1 = (float*)device.api->getStackMemory(216 * num_elements * sizeof(float));
    _tmp0 = d_buffer0;
    {
    unsigned offset_rT = 0;
    unsigned offset_I = I_offset;
    unsigned offset__tmp0 = 216;
    
    
    device.gemm(device::ColMajor, device::NoTrans, device::NoTrans, 24, 9, 56, 1.0, rT(0), 24, I, 56, 0.0, _tmp0, 24, offset_rT, offset_I, offset__tmp0, num_elements);
    }
    _tmp1 = d_buffer1;
    {
    unsigned offset_fP = 0;
    unsigned offset__tmp0 = 216;
    unsigned offset__tmp1 = 216;
    
    
    device.gemm(device::ColMajor, device::NoTrans, device::NoTrans, 24, 9, 21, 1.0, fP(2), 24, _tmp0, 24, 0.0, _tmp1, 24, offset_fP, offset__tmp0, offset__tmp1, num_elements);
    }
    _tmp2 = d_buffer0;
    {
    unsigned offset__tmp1 = 216;
    unsigned offset_AminusT = AminusT_offset;
    unsigned offset__tmp2 = 216;
    
    
    device.gemm(device::ColMajor, device::NoTrans, device::NoTrans, 24, 9, 9, 1.0, _tmp1, 24, AminusT, 9, 0.0, _tmp2, 24, offset__tmp1, offset_AminusT, offset__tmp2, num_elements);
    }
    {
    unsigned offset_rDivM = 0;
    unsigned offset__tmp2 = 216;
    unsigned offset_Q = Q_offset;
    
    
    device.gemm(device::ColMajor, device::NoTrans, device::NoTrans, 56, 9, 21, 1.0, rDivM(2), 56, _tmp2, 24, 1.0, Q, 56, offset_rDivM, offset__tmp2, offset_Q, num_elements);
    }
    device.api->popStackMemory();
    device.api->popStackMemory();
  }
  void kernel::neighboringFlux::execute27() {
    assert(AminusT != nullptr);
    assert(I != nullptr);
    assert(Q != nullptr);
    assert(fP(0) != nullptr);
    assert(rDivM(2) != nullptr);
    assert(rT(1) != nullptr);
    assert(kernel::neighboringFlux::num_elements != 0);
    device::DeviceInstance& device = device::DeviceInstance::getInstance();
    float *_tmp0, *_tmp1, *_tmp2;
    float *d_buffer0 = (float*)device.api->getStackMemory(216 * num_elements * sizeof(float));
    float *d_buffer1 = (float*)device.api->getStackMemory(216 * num_elements * sizeof(float));
    _tmp0 = d_buffer0;
    {
    unsigned offset_rT = 0;
    unsigned offset_I = I_offset;
    unsigned offset__tmp0 = 216;
    
    
    device.gemm(device::ColMajor, device::NoTrans, device::NoTrans, 24, 9, 56, 1.0, rT(1), 24, I, 56, 0.0, _tmp0, 24, offset_rT, offset_I, offset__tmp0, num_elements);
    }
    _tmp1 = d_buffer1;
    {
    unsigned offset_fP = 0;
    unsigned offset__tmp0 = 216;
    unsigned offset__tmp1 = 216;
    
    
    device.gemm(device::ColMajor, device::NoTrans, device::NoTrans, 24, 9, 21, 1.0, fP(0), 24, _tmp0, 24, 0.0, _tmp1, 24, offset_fP, offset__tmp0, offset__tmp1, num_elements);
    }
    _tmp2 = d_buffer0;
    {
    unsigned offset__tmp1 = 216;
    unsigned offset_AminusT = AminusT_offset;
    unsigned offset__tmp2 = 216;
    
    
    device.gemm(device::ColMajor, device::NoTrans, device::NoTrans, 24, 9, 9, 1.0, _tmp1, 24, AminusT, 9, 0.0, _tmp2, 24, offset__tmp1, offset_AminusT, offset__tmp2, num_elements);
    }
    {
    unsigned offset_rDivM = 0;
    unsigned offset__tmp2 = 216;
    unsigned offset_Q = Q_offset;
    
    
    device.gemm(device::ColMajor, device::NoTrans, device::NoTrans, 56, 9, 21, 1.0, rDivM(2), 56, _tmp2, 24, 1.0, Q, 56, offset_rDivM, offset__tmp2, offset_Q, num_elements);
    }
    device.api->popStackMemory();
    device.api->popStackMemory();
  }
  void kernel::neighboringFlux::execute28() {
    assert(AminusT != nullptr);
    assert(I != nullptr);
    assert(Q != nullptr);
    assert(fP(1) != nullptr);
    assert(rDivM(2) != nullptr);
    assert(rT(1) != nullptr);
    assert(kernel::neighboringFlux::num_elements != 0);
    device::DeviceInstance& device = device::DeviceInstance::getInstance();
    float *_tmp0, *_tmp1, *_tmp2;
    float *d_buffer0 = (float*)device.api->getStackMemory(216 * num_elements * sizeof(float));
    float *d_buffer1 = (float*)device.api->getStackMemory(216 * num_elements * sizeof(float));
    _tmp0 = d_buffer0;
    {
    unsigned offset_rT = 0;
    unsigned offset_I = I_offset;
    unsigned offset__tmp0 = 216;
    
    
    device.gemm(device::ColMajor, device::NoTrans, device::NoTrans, 24, 9, 56, 1.0, rT(1), 24, I, 56, 0.0, _tmp0, 24, offset_rT, offset_I, offset__tmp0, num_elements);
    }
    _tmp1 = d_buffer1;
    {
    unsigned offset_fP = 0;
    unsigned offset__tmp0 = 216;
    unsigned offset__tmp1 = 216;
    
    
    device.gemm(device::ColMajor, device::NoTrans, device::NoTrans, 24, 9, 21, 1.0, fP(1), 24, _tmp0, 24, 0.0, _tmp1, 24, offset_fP, offset__tmp0, offset__tmp1, num_elements);
    }
    _tmp2 = d_buffer0;
    {
    unsigned offset__tmp1 = 216;
    unsigned offset_AminusT = AminusT_offset;
    unsigned offset__tmp2 = 216;
    
    
    device.gemm(device::ColMajor, device::NoTrans, device::NoTrans, 24, 9, 9, 1.0, _tmp1, 24, AminusT, 9, 0.0, _tmp2, 24, offset__tmp1, offset_AminusT, offset__tmp2, num_elements);
    }
    {
    unsigned offset_rDivM = 0;
    unsigned offset__tmp2 = 216;
    unsigned offset_Q = Q_offset;
    
    
    device.gemm(device::ColMajor, device::NoTrans, device::NoTrans, 56, 9, 21, 1.0, rDivM(2), 56, _tmp2, 24, 1.0, Q, 56, offset_rDivM, offset__tmp2, offset_Q, num_elements);
    }
    device.api->popStackMemory();
    device.api->popStackMemory();
  }
  void kernel::neighboringFlux::execute29() {
    assert(AminusT != nullptr);
    assert(I != nullptr);
    assert(Q != nullptr);
    assert(fP(2) != nullptr);
    assert(rDivM(2) != nullptr);
    assert(rT(1) != nullptr);
    assert(kernel::neighboringFlux::num_elements != 0);
    device::DeviceInstance& device = device::DeviceInstance::getInstance();
    float *_tmp0, *_tmp1, *_tmp2;
    float *d_buffer0 = (float*)device.api->getStackMemory(216 * num_elements * sizeof(float));
    float *d_buffer1 = (float*)device.api->getStackMemory(216 * num_elements * sizeof(float));
    _tmp0 = d_buffer0;
    {
    unsigned offset_rT = 0;
    unsigned offset_I = I_offset;
    unsigned offset__tmp0 = 216;
    
    
    device.gemm(device::ColMajor, device::NoTrans, device::NoTrans, 24, 9, 56, 1.0, rT(1), 24, I, 56, 0.0, _tmp0, 24, offset_rT, offset_I, offset__tmp0, num_elements);
    }
    _tmp1 = d_buffer1;
    {
    unsigned offset_fP = 0;
    unsigned offset__tmp0 = 216;
    unsigned offset__tmp1 = 216;
    
    
    device.gemm(device::ColMajor, device::NoTrans, device::NoTrans, 24, 9, 21, 1.0, fP(2), 24, _tmp0, 24, 0.0, _tmp1, 24, offset_fP, offset__tmp0, offset__tmp1, num_elements);
    }
    _tmp2 = d_buffer0;
    {
    unsigned offset__tmp1 = 216;
    unsigned offset_AminusT = AminusT_offset;
    unsigned offset__tmp2 = 216;
    
    
    device.gemm(device::ColMajor, device::NoTrans, device::NoTrans, 24, 9, 9, 1.0, _tmp1, 24, AminusT, 9, 0.0, _tmp2, 24, offset__tmp1, offset_AminusT, offset__tmp2, num_elements);
    }
    {
    unsigned offset_rDivM = 0;
    unsigned offset__tmp2 = 216;
    unsigned offset_Q = Q_offset;
    
    
    device.gemm(device::ColMajor, device::NoTrans, device::NoTrans, 56, 9, 21, 1.0, rDivM(2), 56, _tmp2, 24, 1.0, Q, 56, offset_rDivM, offset__tmp2, offset_Q, num_elements);
    }
    device.api->popStackMemory();
    device.api->popStackMemory();
  }
  void kernel::neighboringFlux::execute30() {
    assert(AminusT != nullptr);
    assert(I != nullptr);
    assert(Q != nullptr);
    assert(fP(0) != nullptr);
    assert(rDivM(2) != nullptr);
    assert(rT(2) != nullptr);
    assert(kernel::neighboringFlux::num_elements != 0);
    device::DeviceInstance& device = device::DeviceInstance::getInstance();
    float *_tmp0, *_tmp1, *_tmp2;
    float *d_buffer0 = (float*)device.api->getStackMemory(216 * num_elements * sizeof(float));
    float *d_buffer1 = (float*)device.api->getStackMemory(216 * num_elements * sizeof(float));
    _tmp0 = d_buffer0;
    {
    unsigned offset_rT = 0;
    unsigned offset_I = I_offset;
    unsigned offset__tmp0 = 216;
    
    
    device.gemm(device::ColMajor, device::NoTrans, device::NoTrans, 24, 9, 56, 1.0, rT(2), 24, I, 56, 0.0, _tmp0, 24, offset_rT, offset_I, offset__tmp0, num_elements);
    }
    _tmp1 = d_buffer1;
    {
    unsigned offset_fP = 0;
    unsigned offset__tmp0 = 216;
    unsigned offset__tmp1 = 216;
    
    
    device.gemm(device::ColMajor, device::NoTrans, device::NoTrans, 24, 9, 21, 1.0, fP(0), 24, _tmp0, 24, 0.0, _tmp1, 24, offset_fP, offset__tmp0, offset__tmp1, num_elements);
    }
    _tmp2 = d_buffer0;
    {
    unsigned offset__tmp1 = 216;
    unsigned offset_AminusT = AminusT_offset;
    unsigned offset__tmp2 = 216;
    
    
    device.gemm(device::ColMajor, device::NoTrans, device::NoTrans, 24, 9, 9, 1.0, _tmp1, 24, AminusT, 9, 0.0, _tmp2, 24, offset__tmp1, offset_AminusT, offset__tmp2, num_elements);
    }
    {
    unsigned offset_rDivM = 0;
    unsigned offset__tmp2 = 216;
    unsigned offset_Q = Q_offset;
    
    
    device.gemm(device::ColMajor, device::NoTrans, device::NoTrans, 56, 9, 21, 1.0, rDivM(2), 56, _tmp2, 24, 1.0, Q, 56, offset_rDivM, offset__tmp2, offset_Q, num_elements);
    }
    device.api->popStackMemory();
    device.api->popStackMemory();
  }
  void kernel::neighboringFlux::execute31() {
    assert(AminusT != nullptr);
    assert(I != nullptr);
    assert(Q != nullptr);
    assert(fP(1) != nullptr);
    assert(rDivM(2) != nullptr);
    assert(rT(2) != nullptr);
    assert(kernel::neighboringFlux::num_elements != 0);
    device::DeviceInstance& device = device::DeviceInstance::getInstance();
    float *_tmp0, *_tmp1, *_tmp2;
    float *d_buffer0 = (float*)device.api->getStackMemory(216 * num_elements * sizeof(float));
    float *d_buffer1 = (float*)device.api->getStackMemory(216 * num_elements * sizeof(float));
    _tmp0 = d_buffer0;
    {
    unsigned offset_rT = 0;
    unsigned offset_I = I_offset;
    unsigned offset__tmp0 = 216;
    
    
    device.gemm(device::ColMajor, device::NoTrans, device::NoTrans, 24, 9, 56, 1.0, rT(2), 24, I, 56, 0.0, _tmp0, 24, offset_rT, offset_I, offset__tmp0, num_elements);
    }
    _tmp1 = d_buffer1;
    {
    unsigned offset_fP = 0;
    unsigned offset__tmp0 = 216;
    unsigned offset__tmp1 = 216;
    
    
    device.gemm(device::ColMajor, device::NoTrans, device::NoTrans, 24, 9, 21, 1.0, fP(1), 24, _tmp0, 24, 0.0, _tmp1, 24, offset_fP, offset__tmp0, offset__tmp1, num_elements);
    }
    _tmp2 = d_buffer0;
    {
    unsigned offset__tmp1 = 216;
    unsigned offset_AminusT = AminusT_offset;
    unsigned offset__tmp2 = 216;
    
    
    device.gemm(device::ColMajor, device::NoTrans, device::NoTrans, 24, 9, 9, 1.0, _tmp1, 24, AminusT, 9, 0.0, _tmp2, 24, offset__tmp1, offset_AminusT, offset__tmp2, num_elements);
    }
    {
    unsigned offset_rDivM = 0;
    unsigned offset__tmp2 = 216;
    unsigned offset_Q = Q_offset;
    
    
    device.gemm(device::ColMajor, device::NoTrans, device::NoTrans, 56, 9, 21, 1.0, rDivM(2), 56, _tmp2, 24, 1.0, Q, 56, offset_rDivM, offset__tmp2, offset_Q, num_elements);
    }
    device.api->popStackMemory();
    device.api->popStackMemory();
  }
  void kernel::neighboringFlux::execute32() {
    assert(AminusT != nullptr);
    assert(I != nullptr);
    assert(Q != nullptr);
    assert(fP(2) != nullptr);
    assert(rDivM(2) != nullptr);
    assert(rT(2) != nullptr);
    assert(kernel::neighboringFlux::num_elements != 0);
    device::DeviceInstance& device = device::DeviceInstance::getInstance();
    float *_tmp0, *_tmp1, *_tmp2;
    float *d_buffer0 = (float*)device.api->getStackMemory(216 * num_elements * sizeof(float));
    float *d_buffer1 = (float*)device.api->getStackMemory(216 * num_elements * sizeof(float));
    _tmp0 = d_buffer0;
    {
    unsigned offset_rT = 0;
    unsigned offset_I = I_offset;
    unsigned offset__tmp0 = 216;
    
    
    device.gemm(device::ColMajor, device::NoTrans, device::NoTrans, 24, 9, 56, 1.0, rT(2), 24, I, 56, 0.0, _tmp0, 24, offset_rT, offset_I, offset__tmp0, num_elements);
    }
    _tmp1 = d_buffer1;
    {
    unsigned offset_fP = 0;
    unsigned offset__tmp0 = 216;
    unsigned offset__tmp1 = 216;
    
    
    device.gemm(device::ColMajor, device::NoTrans, device::NoTrans, 24, 9, 21, 1.0, fP(2), 24, _tmp0, 24, 0.0, _tmp1, 24, offset_fP, offset__tmp0, offset__tmp1, num_elements);
    }
    _tmp2 = d_buffer0;
    {
    unsigned offset__tmp1 = 216;
    unsigned offset_AminusT = AminusT_offset;
    unsigned offset__tmp2 = 216;
    
    
    device.gemm(device::ColMajor, device::NoTrans, device::NoTrans, 24, 9, 9, 1.0, _tmp1, 24, AminusT, 9, 0.0, _tmp2, 24, offset__tmp1, offset_AminusT, offset__tmp2, num_elements);
    }
    {
    unsigned offset_rDivM = 0;
    unsigned offset__tmp2 = 216;
    unsigned offset_Q = Q_offset;
    
    
    device.gemm(device::ColMajor, device::NoTrans, device::NoTrans, 56, 9, 21, 1.0, rDivM(2), 56, _tmp2, 24, 1.0, Q, 56, offset_rDivM, offset__tmp2, offset_Q, num_elements);
    }
    device.api->popStackMemory();
    device.api->popStackMemory();
  }
  void kernel::neighboringFlux::execute33() {
    assert(AminusT != nullptr);
    assert(I != nullptr);
    assert(Q != nullptr);
    assert(fP(0) != nullptr);
    assert(rDivM(2) != nullptr);
    assert(rT(3) != nullptr);
    assert(kernel::neighboringFlux::num_elements != 0);
    device::DeviceInstance& device = device::DeviceInstance::getInstance();
    float *_tmp0, *_tmp1, *_tmp2;
    float *d_buffer0 = (float*)device.api->getStackMemory(216 * num_elements * sizeof(float));
    float *d_buffer1 = (float*)device.api->getStackMemory(216 * num_elements * sizeof(float));
    _tmp0 = d_buffer0;
    {
    unsigned offset_rT = 0;
    unsigned offset_I = I_offset;
    unsigned offset__tmp0 = 216;
    
    
    device.gemm(device::ColMajor, device::NoTrans, device::NoTrans, 24, 9, 56, 1.0, rT(3), 24, I, 56, 0.0, _tmp0, 24, offset_rT, offset_I, offset__tmp0, num_elements);
    }
    _tmp1 = d_buffer1;
    {
    unsigned offset_fP = 0;
    unsigned offset__tmp0 = 216;
    unsigned offset__tmp1 = 216;
    
    
    device.gemm(device::ColMajor, device::NoTrans, device::NoTrans, 24, 9, 21, 1.0, fP(0), 24, _tmp0, 24, 0.0, _tmp1, 24, offset_fP, offset__tmp0, offset__tmp1, num_elements);
    }
    _tmp2 = d_buffer0;
    {
    unsigned offset__tmp1 = 216;
    unsigned offset_AminusT = AminusT_offset;
    unsigned offset__tmp2 = 216;
    
    
    device.gemm(device::ColMajor, device::NoTrans, device::NoTrans, 24, 9, 9, 1.0, _tmp1, 24, AminusT, 9, 0.0, _tmp2, 24, offset__tmp1, offset_AminusT, offset__tmp2, num_elements);
    }
    {
    unsigned offset_rDivM = 0;
    unsigned offset__tmp2 = 216;
    unsigned offset_Q = Q_offset;
    
    
    device.gemm(device::ColMajor, device::NoTrans, device::NoTrans, 56, 9, 21, 1.0, rDivM(2), 56, _tmp2, 24, 1.0, Q, 56, offset_rDivM, offset__tmp2, offset_Q, num_elements);
    }
    device.api->popStackMemory();
    device.api->popStackMemory();
  }
  void kernel::neighboringFlux::execute34() {
    assert(AminusT != nullptr);
    assert(I != nullptr);
    assert(Q != nullptr);
    assert(fP(1) != nullptr);
    assert(rDivM(2) != nullptr);
    assert(rT(3) != nullptr);
    assert(kernel::neighboringFlux::num_elements != 0);
    device::DeviceInstance& device = device::DeviceInstance::getInstance();
    float *_tmp0, *_tmp1, *_tmp2;
    float *d_buffer0 = (float*)device.api->getStackMemory(216 * num_elements * sizeof(float));
    float *d_buffer1 = (float*)device.api->getStackMemory(216 * num_elements * sizeof(float));
    _tmp0 = d_buffer0;
    {
    unsigned offset_rT = 0;
    unsigned offset_I = I_offset;
    unsigned offset__tmp0 = 216;
    
    
    device.gemm(device::ColMajor, device::NoTrans, device::NoTrans, 24, 9, 56, 1.0, rT(3), 24, I, 56, 0.0, _tmp0, 24, offset_rT, offset_I, offset__tmp0, num_elements);
    }
    _tmp1 = d_buffer1;
    {
    unsigned offset_fP = 0;
    unsigned offset__tmp0 = 216;
    unsigned offset__tmp1 = 216;
    
    
    device.gemm(device::ColMajor, device::NoTrans, device::NoTrans, 24, 9, 21, 1.0, fP(1), 24, _tmp0, 24, 0.0, _tmp1, 24, offset_fP, offset__tmp0, offset__tmp1, num_elements);
    }
    _tmp2 = d_buffer0;
    {
    unsigned offset__tmp1 = 216;
    unsigned offset_AminusT = AminusT_offset;
    unsigned offset__tmp2 = 216;
    
    
    device.gemm(device::ColMajor, device::NoTrans, device::NoTrans, 24, 9, 9, 1.0, _tmp1, 24, AminusT, 9, 0.0, _tmp2, 24, offset__tmp1, offset_AminusT, offset__tmp2, num_elements);
    }
    {
    unsigned offset_rDivM = 0;
    unsigned offset__tmp2 = 216;
    unsigned offset_Q = Q_offset;
    
    
    device.gemm(device::ColMajor, device::NoTrans, device::NoTrans, 56, 9, 21, 1.0, rDivM(2), 56, _tmp2, 24, 1.0, Q, 56, offset_rDivM, offset__tmp2, offset_Q, num_elements);
    }
    device.api->popStackMemory();
    device.api->popStackMemory();
  }
  void kernel::neighboringFlux::execute35() {
    assert(AminusT != nullptr);
    assert(I != nullptr);
    assert(Q != nullptr);
    assert(fP(2) != nullptr);
    assert(rDivM(2) != nullptr);
    assert(rT(3) != nullptr);
    assert(kernel::neighboringFlux::num_elements != 0);
    device::DeviceInstance& device = device::DeviceInstance::getInstance();
    float *_tmp0, *_tmp1, *_tmp2;
    float *d_buffer0 = (float*)device.api->getStackMemory(216 * num_elements * sizeof(float));
    float *d_buffer1 = (float*)device.api->getStackMemory(216 * num_elements * sizeof(float));
    _tmp0 = d_buffer0;
    {
    unsigned offset_rT = 0;
    unsigned offset_I = I_offset;
    unsigned offset__tmp0 = 216;
    
    
    device.gemm(device::ColMajor, device::NoTrans, device::NoTrans, 24, 9, 56, 1.0, rT(3), 24, I, 56, 0.0, _tmp0, 24, offset_rT, offset_I, offset__tmp0, num_elements);
    }
    _tmp1 = d_buffer1;
    {
    unsigned offset_fP = 0;
    unsigned offset__tmp0 = 216;
    unsigned offset__tmp1 = 216;
    
    
    device.gemm(device::ColMajor, device::NoTrans, device::NoTrans, 24, 9, 21, 1.0, fP(2), 24, _tmp0, 24, 0.0, _tmp1, 24, offset_fP, offset__tmp0, offset__tmp1, num_elements);
    }
    _tmp2 = d_buffer0;
    {
    unsigned offset__tmp1 = 216;
    unsigned offset_AminusT = AminusT_offset;
    unsigned offset__tmp2 = 216;
    
    
    device.gemm(device::ColMajor, device::NoTrans, device::NoTrans, 24, 9, 9, 1.0, _tmp1, 24, AminusT, 9, 0.0, _tmp2, 24, offset__tmp1, offset_AminusT, offset__tmp2, num_elements);
    }
    {
    unsigned offset_rDivM = 0;
    unsigned offset__tmp2 = 216;
    unsigned offset_Q = Q_offset;
    
    
    device.gemm(device::ColMajor, device::NoTrans, device::NoTrans, 56, 9, 21, 1.0, rDivM(2), 56, _tmp2, 24, 1.0, Q, 56, offset_rDivM, offset__tmp2, offset_Q, num_elements);
    }
    device.api->popStackMemory();
    device.api->popStackMemory();
  }
  void kernel::neighboringFlux::execute36() {
    assert(AminusT != nullptr);
    assert(I != nullptr);
    assert(Q != nullptr);
    assert(fP(0) != nullptr);
    assert(rDivM(3) != nullptr);
    assert(rT(0) != nullptr);
    assert(kernel::neighboringFlux::num_elements != 0);
    device::DeviceInstance& device = device::DeviceInstance::getInstance();
    float *_tmp0, *_tmp1, *_tmp2;
    float *d_buffer0 = (float*)device.api->getStackMemory(216 * num_elements * sizeof(float));
    float *d_buffer1 = (float*)device.api->getStackMemory(216 * num_elements * sizeof(float));
    _tmp0 = d_buffer0;
    {
    unsigned offset_rT = 0;
    unsigned offset_I = I_offset;
    unsigned offset__tmp0 = 216;
    
    
    device.gemm(device::ColMajor, device::NoTrans, device::NoTrans, 24, 9, 56, 1.0, rT(0), 24, I, 56, 0.0, _tmp0, 24, offset_rT, offset_I, offset__tmp0, num_elements);
    }
    _tmp1 = d_buffer1;
    {
    unsigned offset_fP = 0;
    unsigned offset__tmp0 = 216;
    unsigned offset__tmp1 = 216;
    
    
    device.gemm(device::ColMajor, device::NoTrans, device::NoTrans, 24, 9, 21, 1.0, fP(0), 24, _tmp0, 24, 0.0, _tmp1, 24, offset_fP, offset__tmp0, offset__tmp1, num_elements);
    }
    _tmp2 = d_buffer0;
    {
    unsigned offset__tmp1 = 216;
    unsigned offset_AminusT = AminusT_offset;
    unsigned offset__tmp2 = 216;
    
    
    device.gemm(device::ColMajor, device::NoTrans, device::NoTrans, 24, 9, 9, 1.0, _tmp1, 24, AminusT, 9, 0.0, _tmp2, 24, offset__tmp1, offset_AminusT, offset__tmp2, num_elements);
    }
    {
    unsigned offset_rDivM = 0;
    unsigned offset__tmp2 = 216;
    unsigned offset_Q = Q_offset;
    
    
    device.gemm(device::ColMajor, device::NoTrans, device::NoTrans, 56, 9, 21, 1.0, rDivM(3), 56, _tmp2, 24, 1.0, Q, 56, offset_rDivM, offset__tmp2, offset_Q, num_elements);
    }
    device.api->popStackMemory();
    device.api->popStackMemory();
  }
  void kernel::neighboringFlux::execute37() {
    assert(AminusT != nullptr);
    assert(I != nullptr);
    assert(Q != nullptr);
    assert(fP(1) != nullptr);
    assert(rDivM(3) != nullptr);
    assert(rT(0) != nullptr);
    assert(kernel::neighboringFlux::num_elements != 0);
    device::DeviceInstance& device = device::DeviceInstance::getInstance();
    float *_tmp0, *_tmp1, *_tmp2;
    float *d_buffer0 = (float*)device.api->getStackMemory(216 * num_elements * sizeof(float));
    float *d_buffer1 = (float*)device.api->getStackMemory(216 * num_elements * sizeof(float));
    _tmp0 = d_buffer0;
    {
    unsigned offset_rT = 0;
    unsigned offset_I = I_offset;
    unsigned offset__tmp0 = 216;
    
    
    device.gemm(device::ColMajor, device::NoTrans, device::NoTrans, 24, 9, 56, 1.0, rT(0), 24, I, 56, 0.0, _tmp0, 24, offset_rT, offset_I, offset__tmp0, num_elements);
    }
    _tmp1 = d_buffer1;
    {
    unsigned offset_fP = 0;
    unsigned offset__tmp0 = 216;
    unsigned offset__tmp1 = 216;
    
    
    device.gemm(device::ColMajor, device::NoTrans, device::NoTrans, 24, 9, 21, 1.0, fP(1), 24, _tmp0, 24, 0.0, _tmp1, 24, offset_fP, offset__tmp0, offset__tmp1, num_elements);
    }
    _tmp2 = d_buffer0;
    {
    unsigned offset__tmp1 = 216;
    unsigned offset_AminusT = AminusT_offset;
    unsigned offset__tmp2 = 216;
    
    
    device.gemm(device::ColMajor, device::NoTrans, device::NoTrans, 24, 9, 9, 1.0, _tmp1, 24, AminusT, 9, 0.0, _tmp2, 24, offset__tmp1, offset_AminusT, offset__tmp2, num_elements);
    }
    {
    unsigned offset_rDivM = 0;
    unsigned offset__tmp2 = 216;
    unsigned offset_Q = Q_offset;
    
    
    device.gemm(device::ColMajor, device::NoTrans, device::NoTrans, 56, 9, 21, 1.0, rDivM(3), 56, _tmp2, 24, 1.0, Q, 56, offset_rDivM, offset__tmp2, offset_Q, num_elements);
    }
    device.api->popStackMemory();
    device.api->popStackMemory();
  }
  void kernel::neighboringFlux::execute38() {
    assert(AminusT != nullptr);
    assert(I != nullptr);
    assert(Q != nullptr);
    assert(fP(2) != nullptr);
    assert(rDivM(3) != nullptr);
    assert(rT(0) != nullptr);
    assert(kernel::neighboringFlux::num_elements != 0);
    device::DeviceInstance& device = device::DeviceInstance::getInstance();
    float *_tmp0, *_tmp1, *_tmp2;
    float *d_buffer0 = (float*)device.api->getStackMemory(216 * num_elements * sizeof(float));
    float *d_buffer1 = (float*)device.api->getStackMemory(216 * num_elements * sizeof(float));
    _tmp0 = d_buffer0;
    {
    unsigned offset_rT = 0;
    unsigned offset_I = I_offset;
    unsigned offset__tmp0 = 216;
    
    
    device.gemm(device::ColMajor, device::NoTrans, device::NoTrans, 24, 9, 56, 1.0, rT(0), 24, I, 56, 0.0, _tmp0, 24, offset_rT, offset_I, offset__tmp0, num_elements);
    }
    _tmp1 = d_buffer1;
    {
    unsigned offset_fP = 0;
    unsigned offset__tmp0 = 216;
    unsigned offset__tmp1 = 216;
    
    
    device.gemm(device::ColMajor, device::NoTrans, device::NoTrans, 24, 9, 21, 1.0, fP(2), 24, _tmp0, 24, 0.0, _tmp1, 24, offset_fP, offset__tmp0, offset__tmp1, num_elements);
    }
    _tmp2 = d_buffer0;
    {
    unsigned offset__tmp1 = 216;
    unsigned offset_AminusT = AminusT_offset;
    unsigned offset__tmp2 = 216;
    
    
    device.gemm(device::ColMajor, device::NoTrans, device::NoTrans, 24, 9, 9, 1.0, _tmp1, 24, AminusT, 9, 0.0, _tmp2, 24, offset__tmp1, offset_AminusT, offset__tmp2, num_elements);
    }
    {
    unsigned offset_rDivM = 0;
    unsigned offset__tmp2 = 216;
    unsigned offset_Q = Q_offset;
    
    
    device.gemm(device::ColMajor, device::NoTrans, device::NoTrans, 56, 9, 21, 1.0, rDivM(3), 56, _tmp2, 24, 1.0, Q, 56, offset_rDivM, offset__tmp2, offset_Q, num_elements);
    }
    device.api->popStackMemory();
    device.api->popStackMemory();
  }
  void kernel::neighboringFlux::execute39() {
    assert(AminusT != nullptr);
    assert(I != nullptr);
    assert(Q != nullptr);
    assert(fP(0) != nullptr);
    assert(rDivM(3) != nullptr);
    assert(rT(1) != nullptr);
    assert(kernel::neighboringFlux::num_elements != 0);
    device::DeviceInstance& device = device::DeviceInstance::getInstance();
    float *_tmp0, *_tmp1, *_tmp2;
    float *d_buffer0 = (float*)device.api->getStackMemory(216 * num_elements * sizeof(float));
    float *d_buffer1 = (float*)device.api->getStackMemory(216 * num_elements * sizeof(float));
    _tmp0 = d_buffer0;
    {
    unsigned offset_rT = 0;
    unsigned offset_I = I_offset;
    unsigned offset__tmp0 = 216;
    
    
    device.gemm(device::ColMajor, device::NoTrans, device::NoTrans, 24, 9, 56, 1.0, rT(1), 24, I, 56, 0.0, _tmp0, 24, offset_rT, offset_I, offset__tmp0, num_elements);
    }
    _tmp1 = d_buffer1;
    {
    unsigned offset_fP = 0;
    unsigned offset__tmp0 = 216;
    unsigned offset__tmp1 = 216;
    
    
    device.gemm(device::ColMajor, device::NoTrans, device::NoTrans, 24, 9, 21, 1.0, fP(0), 24, _tmp0, 24, 0.0, _tmp1, 24, offset_fP, offset__tmp0, offset__tmp1, num_elements);
    }
    _tmp2 = d_buffer0;
    {
    unsigned offset__tmp1 = 216;
    unsigned offset_AminusT = AminusT_offset;
    unsigned offset__tmp2 = 216;
    
    
    device.gemm(device::ColMajor, device::NoTrans, device::NoTrans, 24, 9, 9, 1.0, _tmp1, 24, AminusT, 9, 0.0, _tmp2, 24, offset__tmp1, offset_AminusT, offset__tmp2, num_elements);
    }
    {
    unsigned offset_rDivM = 0;
    unsigned offset__tmp2 = 216;
    unsigned offset_Q = Q_offset;
    
    
    device.gemm(device::ColMajor, device::NoTrans, device::NoTrans, 56, 9, 21, 1.0, rDivM(3), 56, _tmp2, 24, 1.0, Q, 56, offset_rDivM, offset__tmp2, offset_Q, num_elements);
    }
    device.api->popStackMemory();
    device.api->popStackMemory();
  }
  void kernel::neighboringFlux::execute40() {
    assert(AminusT != nullptr);
    assert(I != nullptr);
    assert(Q != nullptr);
    assert(fP(1) != nullptr);
    assert(rDivM(3) != nullptr);
    assert(rT(1) != nullptr);
    assert(kernel::neighboringFlux::num_elements != 0);
    device::DeviceInstance& device = device::DeviceInstance::getInstance();
    float *_tmp0, *_tmp1, *_tmp2;
    float *d_buffer0 = (float*)device.api->getStackMemory(216 * num_elements * sizeof(float));
    float *d_buffer1 = (float*)device.api->getStackMemory(216 * num_elements * sizeof(float));
    _tmp0 = d_buffer0;
    {
    unsigned offset_rT = 0;
    unsigned offset_I = I_offset;
    unsigned offset__tmp0 = 216;
    
    
    device.gemm(device::ColMajor, device::NoTrans, device::NoTrans, 24, 9, 56, 1.0, rT(1), 24, I, 56, 0.0, _tmp0, 24, offset_rT, offset_I, offset__tmp0, num_elements);
    }
    _tmp1 = d_buffer1;
    {
    unsigned offset_fP = 0;
    unsigned offset__tmp0 = 216;
    unsigned offset__tmp1 = 216;
    
    
    device.gemm(device::ColMajor, device::NoTrans, device::NoTrans, 24, 9, 21, 1.0, fP(1), 24, _tmp0, 24, 0.0, _tmp1, 24, offset_fP, offset__tmp0, offset__tmp1, num_elements);
    }
    _tmp2 = d_buffer0;
    {
    unsigned offset__tmp1 = 216;
    unsigned offset_AminusT = AminusT_offset;
    unsigned offset__tmp2 = 216;
    
    
    device.gemm(device::ColMajor, device::NoTrans, device::NoTrans, 24, 9, 9, 1.0, _tmp1, 24, AminusT, 9, 0.0, _tmp2, 24, offset__tmp1, offset_AminusT, offset__tmp2, num_elements);
    }
    {
    unsigned offset_rDivM = 0;
    unsigned offset__tmp2 = 216;
    unsigned offset_Q = Q_offset;
    
    
    device.gemm(device::ColMajor, device::NoTrans, device::NoTrans, 56, 9, 21, 1.0, rDivM(3), 56, _tmp2, 24, 1.0, Q, 56, offset_rDivM, offset__tmp2, offset_Q, num_elements);
    }
    device.api->popStackMemory();
    device.api->popStackMemory();
  }
  void kernel::neighboringFlux::execute41() {
    assert(AminusT != nullptr);
    assert(I != nullptr);
    assert(Q != nullptr);
    assert(fP(2) != nullptr);
    assert(rDivM(3) != nullptr);
    assert(rT(1) != nullptr);
    assert(kernel::neighboringFlux::num_elements != 0);
    device::DeviceInstance& device = device::DeviceInstance::getInstance();
    float *_tmp0, *_tmp1, *_tmp2;
    float *d_buffer0 = (float*)device.api->getStackMemory(216 * num_elements * sizeof(float));
    float *d_buffer1 = (float*)device.api->getStackMemory(216 * num_elements * sizeof(float));
    _tmp0 = d_buffer0;
    {
    unsigned offset_rT = 0;
    unsigned offset_I = I_offset;
    unsigned offset__tmp0 = 216;
    
    
    device.gemm(device::ColMajor, device::NoTrans, device::NoTrans, 24, 9, 56, 1.0, rT(1), 24, I, 56, 0.0, _tmp0, 24, offset_rT, offset_I, offset__tmp0, num_elements);
    }
    _tmp1 = d_buffer1;
    {
    unsigned offset_fP = 0;
    unsigned offset__tmp0 = 216;
    unsigned offset__tmp1 = 216;
    
    
    device.gemm(device::ColMajor, device::NoTrans, device::NoTrans, 24, 9, 21, 1.0, fP(2), 24, _tmp0, 24, 0.0, _tmp1, 24, offset_fP, offset__tmp0, offset__tmp1, num_elements);
    }
    _tmp2 = d_buffer0;
    {
    unsigned offset__tmp1 = 216;
    unsigned offset_AminusT = AminusT_offset;
    unsigned offset__tmp2 = 216;
    
    
    device.gemm(device::ColMajor, device::NoTrans, device::NoTrans, 24, 9, 9, 1.0, _tmp1, 24, AminusT, 9, 0.0, _tmp2, 24, offset__tmp1, offset_AminusT, offset__tmp2, num_elements);
    }
    {
    unsigned offset_rDivM = 0;
    unsigned offset__tmp2 = 216;
    unsigned offset_Q = Q_offset;
    
    
    device.gemm(device::ColMajor, device::NoTrans, device::NoTrans, 56, 9, 21, 1.0, rDivM(3), 56, _tmp2, 24, 1.0, Q, 56, offset_rDivM, offset__tmp2, offset_Q, num_elements);
    }
    device.api->popStackMemory();
    device.api->popStackMemory();
  }
  void kernel::neighboringFlux::execute42() {
    assert(AminusT != nullptr);
    assert(I != nullptr);
    assert(Q != nullptr);
    assert(fP(0) != nullptr);
    assert(rDivM(3) != nullptr);
    assert(rT(2) != nullptr);
    assert(kernel::neighboringFlux::num_elements != 0);
    device::DeviceInstance& device = device::DeviceInstance::getInstance();
    float *_tmp0, *_tmp1, *_tmp2;
    float *d_buffer0 = (float*)device.api->getStackMemory(216 * num_elements * sizeof(float));
    float *d_buffer1 = (float*)device.api->getStackMemory(216 * num_elements * sizeof(float));
    _tmp0 = d_buffer0;
    {
    unsigned offset_rT = 0;
    unsigned offset_I = I_offset;
    unsigned offset__tmp0 = 216;
    
    
    device.gemm(device::ColMajor, device::NoTrans, device::NoTrans, 24, 9, 56, 1.0, rT(2), 24, I, 56, 0.0, _tmp0, 24, offset_rT, offset_I, offset__tmp0, num_elements);
    }
    _tmp1 = d_buffer1;
    {
    unsigned offset_fP = 0;
    unsigned offset__tmp0 = 216;
    unsigned offset__tmp1 = 216;
    
    
    device.gemm(device::ColMajor, device::NoTrans, device::NoTrans, 24, 9, 21, 1.0, fP(0), 24, _tmp0, 24, 0.0, _tmp1, 24, offset_fP, offset__tmp0, offset__tmp1, num_elements);
    }
    _tmp2 = d_buffer0;
    {
    unsigned offset__tmp1 = 216;
    unsigned offset_AminusT = AminusT_offset;
    unsigned offset__tmp2 = 216;
    
    
    device.gemm(device::ColMajor, device::NoTrans, device::NoTrans, 24, 9, 9, 1.0, _tmp1, 24, AminusT, 9, 0.0, _tmp2, 24, offset__tmp1, offset_AminusT, offset__tmp2, num_elements);
    }
    {
    unsigned offset_rDivM = 0;
    unsigned offset__tmp2 = 216;
    unsigned offset_Q = Q_offset;
    
    
    device.gemm(device::ColMajor, device::NoTrans, device::NoTrans, 56, 9, 21, 1.0, rDivM(3), 56, _tmp2, 24, 1.0, Q, 56, offset_rDivM, offset__tmp2, offset_Q, num_elements);
    }
    device.api->popStackMemory();
    device.api->popStackMemory();
  }
  void kernel::neighboringFlux::execute43() {
    assert(AminusT != nullptr);
    assert(I != nullptr);
    assert(Q != nullptr);
    assert(fP(1) != nullptr);
    assert(rDivM(3) != nullptr);
    assert(rT(2) != nullptr);
    assert(kernel::neighboringFlux::num_elements != 0);
    device::DeviceInstance& device = device::DeviceInstance::getInstance();
    float *_tmp0, *_tmp1, *_tmp2;
    float *d_buffer0 = (float*)device.api->getStackMemory(216 * num_elements * sizeof(float));
    float *d_buffer1 = (float*)device.api->getStackMemory(216 * num_elements * sizeof(float));
    _tmp0 = d_buffer0;
    {
    unsigned offset_rT = 0;
    unsigned offset_I = I_offset;
    unsigned offset__tmp0 = 216;
    
    
    device.gemm(device::ColMajor, device::NoTrans, device::NoTrans, 24, 9, 56, 1.0, rT(2), 24, I, 56, 0.0, _tmp0, 24, offset_rT, offset_I, offset__tmp0, num_elements);
    }
    _tmp1 = d_buffer1;
    {
    unsigned offset_fP = 0;
    unsigned offset__tmp0 = 216;
    unsigned offset__tmp1 = 216;
    
    
    device.gemm(device::ColMajor, device::NoTrans, device::NoTrans, 24, 9, 21, 1.0, fP(1), 24, _tmp0, 24, 0.0, _tmp1, 24, offset_fP, offset__tmp0, offset__tmp1, num_elements);
    }
    _tmp2 = d_buffer0;
    {
    unsigned offset__tmp1 = 216;
    unsigned offset_AminusT = AminusT_offset;
    unsigned offset__tmp2 = 216;
    
    
    device.gemm(device::ColMajor, device::NoTrans, device::NoTrans, 24, 9, 9, 1.0, _tmp1, 24, AminusT, 9, 0.0, _tmp2, 24, offset__tmp1, offset_AminusT, offset__tmp2, num_elements);
    }
    {
    unsigned offset_rDivM = 0;
    unsigned offset__tmp2 = 216;
    unsigned offset_Q = Q_offset;
    
    
    device.gemm(device::ColMajor, device::NoTrans, device::NoTrans, 56, 9, 21, 1.0, rDivM(3), 56, _tmp2, 24, 1.0, Q, 56, offset_rDivM, offset__tmp2, offset_Q, num_elements);
    }
    device.api->popStackMemory();
    device.api->popStackMemory();
  }
  void kernel::neighboringFlux::execute44() {
    assert(AminusT != nullptr);
    assert(I != nullptr);
    assert(Q != nullptr);
    assert(fP(2) != nullptr);
    assert(rDivM(3) != nullptr);
    assert(rT(2) != nullptr);
    assert(kernel::neighboringFlux::num_elements != 0);
    device::DeviceInstance& device = device::DeviceInstance::getInstance();
    float *_tmp0, *_tmp1, *_tmp2;
    float *d_buffer0 = (float*)device.api->getStackMemory(216 * num_elements * sizeof(float));
    float *d_buffer1 = (float*)device.api->getStackMemory(216 * num_elements * sizeof(float));
    _tmp0 = d_buffer0;
    {
    unsigned offset_rT = 0;
    unsigned offset_I = I_offset;
    unsigned offset__tmp0 = 216;
    
    
    device.gemm(device::ColMajor, device::NoTrans, device::NoTrans, 24, 9, 56, 1.0, rT(2), 24, I, 56, 0.0, _tmp0, 24, offset_rT, offset_I, offset__tmp0, num_elements);
    }
    _tmp1 = d_buffer1;
    {
    unsigned offset_fP = 0;
    unsigned offset__tmp0 = 216;
    unsigned offset__tmp1 = 216;
    
    
    device.gemm(device::ColMajor, device::NoTrans, device::NoTrans, 24, 9, 21, 1.0, fP(2), 24, _tmp0, 24, 0.0, _tmp1, 24, offset_fP, offset__tmp0, offset__tmp1, num_elements);
    }
    _tmp2 = d_buffer0;
    {
    unsigned offset__tmp1 = 216;
    unsigned offset_AminusT = AminusT_offset;
    unsigned offset__tmp2 = 216;
    
    
    device.gemm(device::ColMajor, device::NoTrans, device::NoTrans, 24, 9, 9, 1.0, _tmp1, 24, AminusT, 9, 0.0, _tmp2, 24, offset__tmp1, offset_AminusT, offset__tmp2, num_elements);
    }
    {
    unsigned offset_rDivM = 0;
    unsigned offset__tmp2 = 216;
    unsigned offset_Q = Q_offset;
    
    
    device.gemm(device::ColMajor, device::NoTrans, device::NoTrans, 56, 9, 21, 1.0, rDivM(3), 56, _tmp2, 24, 1.0, Q, 56, offset_rDivM, offset__tmp2, offset_Q, num_elements);
    }
    device.api->popStackMemory();
    device.api->popStackMemory();
  }
  void kernel::neighboringFlux::execute45() {
    assert(AminusT != nullptr);
    assert(I != nullptr);
    assert(Q != nullptr);
    assert(fP(0) != nullptr);
    assert(rDivM(3) != nullptr);
    assert(rT(3) != nullptr);
    assert(kernel::neighboringFlux::num_elements != 0);
    device::DeviceInstance& device = device::DeviceInstance::getInstance();
    float *_tmp0, *_tmp1, *_tmp2;
    float *d_buffer0 = (float*)device.api->getStackMemory(216 * num_elements * sizeof(float));
    float *d_buffer1 = (float*)device.api->getStackMemory(216 * num_elements * sizeof(float));
    _tmp0 = d_buffer0;
    {
    unsigned offset_rT = 0;
    unsigned offset_I = I_offset;
    unsigned offset__tmp0 = 216;
    
    
    device.gemm(device::ColMajor, device::NoTrans, device::NoTrans, 24, 9, 56, 1.0, rT(3), 24, I, 56, 0.0, _tmp0, 24, offset_rT, offset_I, offset__tmp0, num_elements);
    }
    _tmp1 = d_buffer1;
    {
    unsigned offset_fP = 0;
    unsigned offset__tmp0 = 216;
    unsigned offset__tmp1 = 216;
    
    
    device.gemm(device::ColMajor, device::NoTrans, device::NoTrans, 24, 9, 21, 1.0, fP(0), 24, _tmp0, 24, 0.0, _tmp1, 24, offset_fP, offset__tmp0, offset__tmp1, num_elements);
    }
    _tmp2 = d_buffer0;
    {
    unsigned offset__tmp1 = 216;
    unsigned offset_AminusT = AminusT_offset;
    unsigned offset__tmp2 = 216;
    
    
    device.gemm(device::ColMajor, device::NoTrans, device::NoTrans, 24, 9, 9, 1.0, _tmp1, 24, AminusT, 9, 0.0, _tmp2, 24, offset__tmp1, offset_AminusT, offset__tmp2, num_elements);
    }
    {
    unsigned offset_rDivM = 0;
    unsigned offset__tmp2 = 216;
    unsigned offset_Q = Q_offset;
    
    
    device.gemm(device::ColMajor, device::NoTrans, device::NoTrans, 56, 9, 21, 1.0, rDivM(3), 56, _tmp2, 24, 1.0, Q, 56, offset_rDivM, offset__tmp2, offset_Q, num_elements);
    }
    device.api->popStackMemory();
    device.api->popStackMemory();
  }
  void kernel::neighboringFlux::execute46() {
    assert(AminusT != nullptr);
    assert(I != nullptr);
    assert(Q != nullptr);
    assert(fP(1) != nullptr);
    assert(rDivM(3) != nullptr);
    assert(rT(3) != nullptr);
    assert(kernel::neighboringFlux::num_elements != 0);
    device::DeviceInstance& device = device::DeviceInstance::getInstance();
    float *_tmp0, *_tmp1, *_tmp2;
    float *d_buffer0 = (float*)device.api->getStackMemory(216 * num_elements * sizeof(float));
    float *d_buffer1 = (float*)device.api->getStackMemory(216 * num_elements * sizeof(float));
    _tmp0 = d_buffer0;
    {
    unsigned offset_rT = 0;
    unsigned offset_I = I_offset;
    unsigned offset__tmp0 = 216;
    
    
    device.gemm(device::ColMajor, device::NoTrans, device::NoTrans, 24, 9, 56, 1.0, rT(3), 24, I, 56, 0.0, _tmp0, 24, offset_rT, offset_I, offset__tmp0, num_elements);
    }
    _tmp1 = d_buffer1;
    {
    unsigned offset_fP = 0;
    unsigned offset__tmp0 = 216;
    unsigned offset__tmp1 = 216;
    
    
    device.gemm(device::ColMajor, device::NoTrans, device::NoTrans, 24, 9, 21, 1.0, fP(1), 24, _tmp0, 24, 0.0, _tmp1, 24, offset_fP, offset__tmp0, offset__tmp1, num_elements);
    }
    _tmp2 = d_buffer0;
    {
    unsigned offset__tmp1 = 216;
    unsigned offset_AminusT = AminusT_offset;
    unsigned offset__tmp2 = 216;
    
    
    device.gemm(device::ColMajor, device::NoTrans, device::NoTrans, 24, 9, 9, 1.0, _tmp1, 24, AminusT, 9, 0.0, _tmp2, 24, offset__tmp1, offset_AminusT, offset__tmp2, num_elements);
    }
    {
    unsigned offset_rDivM = 0;
    unsigned offset__tmp2 = 216;
    unsigned offset_Q = Q_offset;
    
    
    device.gemm(device::ColMajor, device::NoTrans, device::NoTrans, 56, 9, 21, 1.0, rDivM(3), 56, _tmp2, 24, 1.0, Q, 56, offset_rDivM, offset__tmp2, offset_Q, num_elements);
    }
    device.api->popStackMemory();
    device.api->popStackMemory();
  }
  void kernel::neighboringFlux::execute47() {
    assert(AminusT != nullptr);
    assert(I != nullptr);
    assert(Q != nullptr);
    assert(fP(2) != nullptr);
    assert(rDivM(3) != nullptr);
    assert(rT(3) != nullptr);
    assert(kernel::neighboringFlux::num_elements != 0);
    device::DeviceInstance& device = device::DeviceInstance::getInstance();
    float *_tmp0, *_tmp1, *_tmp2;
    float *d_buffer0 = (float*)device.api->getStackMemory(216 * num_elements * sizeof(float));
    float *d_buffer1 = (float*)device.api->getStackMemory(216 * num_elements * sizeof(float));
    _tmp0 = d_buffer0;
    {
    unsigned offset_rT = 0;
    unsigned offset_I = I_offset;
    unsigned offset__tmp0 = 216;
    
    
    device.gemm(device::ColMajor, device::NoTrans, device::NoTrans, 24, 9, 56, 1.0, rT(3), 24, I, 56, 0.0, _tmp0, 24, offset_rT, offset_I, offset__tmp0, num_elements);
    }
    _tmp1 = d_buffer1;
    {
    unsigned offset_fP = 0;
    unsigned offset__tmp0 = 216;
    unsigned offset__tmp1 = 216;
    
    
    device.gemm(device::ColMajor, device::NoTrans, device::NoTrans, 24, 9, 21, 1.0, fP(2), 24, _tmp0, 24, 0.0, _tmp1, 24, offset_fP, offset__tmp0, offset__tmp1, num_elements);
    }
    _tmp2 = d_buffer0;
    {
    unsigned offset__tmp1 = 216;
    unsigned offset_AminusT = AminusT_offset;
    unsigned offset__tmp2 = 216;
    
    
    device.gemm(device::ColMajor, device::NoTrans, device::NoTrans, 24, 9, 9, 1.0, _tmp1, 24, AminusT, 9, 0.0, _tmp2, 24, offset__tmp1, offset_AminusT, offset__tmp2, num_elements);
    }
    {
    unsigned offset_rDivM = 0;
    unsigned offset__tmp2 = 216;
    unsigned offset_Q = Q_offset;
    
    
    device.gemm(device::ColMajor, device::NoTrans, device::NoTrans, 56, 9, 21, 1.0, rDivM(3), 56, _tmp2, 24, 1.0, Q, 56, offset_rDivM, offset__tmp2, offset_Q, num_elements);
    }
    device.api->popStackMemory();
    device.api->popStackMemory();
  }
  constexpr unsigned long const kernel::derivativeTaylorExpansion::NonZeroFlops[];
  constexpr unsigned long const kernel::derivativeTaylorExpansion::HardwareFlops[];
  constexpr kernel::derivativeTaylorExpansion::member_function_ptr kernel::derivativeTaylorExpansion::ExecutePtrs[];
  void kernel::derivativeTaylorExpansion::execute0() {
    assert(!std::isnan(power));
    assert(I != nullptr);
    assert(dQ(0) != nullptr);
    assert(kernel::derivativeTaylorExpansion::num_elements != 0);
    device::DeviceInstance& device = device::DeviceInstance::getInstance();
    device.copyAddScale(56, 9, power, dQ(0) + 0, 56, 0.0, I + 0, 56, dQ_offset(0), I_offset, num_elements);
  }
  void kernel::derivativeTaylorExpansion::execute1() {
    assert(!std::isnan(power));
    assert(I != nullptr);
    assert(dQ(1) != nullptr);
    assert(kernel::derivativeTaylorExpansion::num_elements != 0);
    device::DeviceInstance& device = device::DeviceInstance::getInstance();
    device.copyAddScale(35, 9, power, dQ(1) + 0, 40, 1.0, I + 0, 56, dQ_offset(1), I_offset, num_elements);
  }
  void kernel::derivativeTaylorExpansion::execute2() {
    assert(!std::isnan(power));
    assert(I != nullptr);
    assert(dQ(2) != nullptr);
    assert(kernel::derivativeTaylorExpansion::num_elements != 0);
    device::DeviceInstance& device = device::DeviceInstance::getInstance();
    device.copyAddScale(20, 9, power, dQ(2) + 0, 24, 1.0, I + 0, 56, dQ_offset(2), I_offset, num_elements);
  }
  void kernel::derivativeTaylorExpansion::execute3() {
    assert(!std::isnan(power));
    assert(I != nullptr);
    assert(dQ(3) != nullptr);
    assert(kernel::derivativeTaylorExpansion::num_elements != 0);
    device::DeviceInstance& device = device::DeviceInstance::getInstance();
    device.copyAddScale(10, 9, power, dQ(3) + 0, 16, 1.0, I + 0, 56, dQ_offset(3), I_offset, num_elements);
  }
  void kernel::derivativeTaylorExpansion::execute4() {
    assert(!std::isnan(power));
    assert(I != nullptr);
    assert(dQ(4) != nullptr);
    assert(kernel::derivativeTaylorExpansion::num_elements != 0);
    device::DeviceInstance& device = device::DeviceInstance::getInstance();
    device.copyAddScale(4, 9, power, dQ(4) + 0, 8, 1.0, I + 0, 56, dQ_offset(4), I_offset, num_elements);
  }
  void kernel::derivativeTaylorExpansion::execute5() {
    assert(!std::isnan(power));
    assert(I != nullptr);
    assert(dQ(5) != nullptr);
    assert(kernel::derivativeTaylorExpansion::num_elements != 0);
    device::DeviceInstance& device = device::DeviceInstance::getInstance();
    device.copyAddScale(1, 9, power, dQ(5) + 0, 8, 1.0, I + 0, 56, dQ_offset(5), I_offset, num_elements);
  }
  constexpr unsigned long const kernel::derivative::NonZeroFlops[];
  constexpr unsigned long const kernel::derivative::HardwareFlops[];
  constexpr kernel::derivative::member_function_ptr kernel::derivative::ExecutePtrs[];
  void kernel::derivative::execute1() {
    assert(dQ(0) != nullptr);
    assert(dQ(1) != nullptr);
    assert(kDivMT(2) != nullptr);
    assert(kDivMT(0) != nullptr);
    assert(kDivMT(1) != nullptr);
    assert(star(2) != nullptr);
    assert(star(0) != nullptr);
    assert(star(1) != nullptr);
    assert(kernel::derivative::num_elements != 0);
    device::DeviceInstance& device = device::DeviceInstance::getInstance();
    float *_tmp0, *_tmp2, *_tmp4;
    float *d_buffer0 = (float*)device.api->getStackMemory(360 * num_elements * sizeof(float));
    _tmp0 = d_buffer0;
    {
    unsigned offset_kDivMT = 0;
    unsigned offset_dQ = dQ_offset(0);
    unsigned offset__tmp0 = 360;
    
    
    device.gemm(device::ColMajor, device::NoTrans, device::NoTrans, 40, 9, 53, 1.0, kDivMT(0), 40, dQ(0), 56, 0.0, _tmp0, 40, offset_kDivMT, offset_dQ + 1, offset__tmp0, num_elements);
    }
    {
    unsigned offset__tmp0 = 360;
    unsigned offset_star = star_offset(0);
    unsigned offset_dQ = dQ_offset(1);
    
    
    device.gemm(device::ColMajor, device::NoTrans, device::NoTrans, 40, 9, 9, 1.0, _tmp0, 40, star(0), 9, 0.0, dQ(1), 40, offset__tmp0, offset_star, offset_dQ, num_elements);
    }
    _tmp2 = d_buffer0;
    {
    unsigned offset_kDivMT = 0;
    unsigned offset_dQ = dQ_offset(0);
    unsigned offset__tmp2 = 360;
    
    
    device.gemm(device::ColMajor, device::NoTrans, device::NoTrans, 40, 9, 54, 1.0, kDivMT(1), 40, dQ(0), 56, 0.0, _tmp2, 40, offset_kDivMT, offset_dQ + 1, offset__tmp2, num_elements);
    }
    {
    unsigned offset__tmp2 = 360;
    unsigned offset_star = star_offset(1);
    unsigned offset_dQ = dQ_offset(1);
    
    
    device.gemm(device::ColMajor, device::NoTrans, device::NoTrans, 40, 9, 9, 1.0, _tmp2, 40, star(1), 9, 1.0, dQ(1), 40, offset__tmp2, offset_star, offset_dQ, num_elements);
    }
    _tmp4 = d_buffer0;
    {
    unsigned offset_kDivMT = 0;
    unsigned offset_dQ = dQ_offset(0);
    unsigned offset__tmp4 = 360;
    
    
    device.gemm(device::ColMajor, device::NoTrans, device::NoTrans, 40, 9, 55, 1.0, kDivMT(2), 40, dQ(0), 56, 0.0, _tmp4, 40, offset_kDivMT, offset_dQ + 1, offset__tmp4, num_elements);
    }
    {
    unsigned offset__tmp4 = 360;
    unsigned offset_star = star_offset(2);
    unsigned offset_dQ = dQ_offset(1);
    
    
    device.gemm(device::ColMajor, device::NoTrans, device::NoTrans, 40, 9, 9, 1.0, _tmp4, 40, star(2), 9, 1.0, dQ(1), 40, offset__tmp4, offset_star, offset_dQ, num_elements);
    }
    device.api->popStackMemory();
  }
  void kernel::derivative::execute2() {
    assert(dQ(2) != nullptr);
    assert(dQ(1) != nullptr);
    assert(kDivMT(2) != nullptr);
    assert(kDivMT(0) != nullptr);
    assert(kDivMT(1) != nullptr);
    assert(star(2) != nullptr);
    assert(star(0) != nullptr);
    assert(star(1) != nullptr);
    assert(kernel::derivative::num_elements != 0);
    device::DeviceInstance& device = device::DeviceInstance::getInstance();
    float *_tmp0, *_tmp2, *_tmp4;
    float *d_buffer0 = (float*)device.api->getStackMemory(216 * num_elements * sizeof(float));
    _tmp0 = d_buffer0;
    {
    unsigned offset_kDivMT = 0;
    unsigned offset_dQ = dQ_offset(1);
    unsigned offset__tmp0 = 216;
    
    
    device.gemm(device::ColMajor, device::NoTrans, device::NoTrans, 24, 9, 32, 1.0, kDivMT(0), 40, dQ(1), 40, 0.0, _tmp0, 24, offset_kDivMT, offset_dQ + 1, offset__tmp0, num_elements);
    }
    {
    unsigned offset__tmp0 = 216;
    unsigned offset_star = star_offset(0);
    unsigned offset_dQ = dQ_offset(2);
    
    
    device.gemm(device::ColMajor, device::NoTrans, device::NoTrans, 24, 9, 9, 1.0, _tmp0, 24, star(0), 9, 0.0, dQ(2), 24, offset__tmp0, offset_star, offset_dQ, num_elements);
    }
    _tmp2 = d_buffer0;
    {
    unsigned offset_kDivMT = 0;
    unsigned offset_dQ = dQ_offset(1);
    unsigned offset__tmp2 = 216;
    
    
    device.gemm(device::ColMajor, device::NoTrans, device::NoTrans, 24, 9, 33, 1.0, kDivMT(1), 40, dQ(1), 40, 0.0, _tmp2, 24, offset_kDivMT, offset_dQ + 1, offset__tmp2, num_elements);
    }
    {
    unsigned offset__tmp2 = 216;
    unsigned offset_star = star_offset(1);
    unsigned offset_dQ = dQ_offset(2);
    
    
    device.gemm(device::ColMajor, device::NoTrans, device::NoTrans, 24, 9, 9, 1.0, _tmp2, 24, star(1), 9, 1.0, dQ(2), 24, offset__tmp2, offset_star, offset_dQ, num_elements);
    }
    _tmp4 = d_buffer0;
    {
    unsigned offset_kDivMT = 0;
    unsigned offset_dQ = dQ_offset(1);
    unsigned offset__tmp4 = 216;
    
    
    device.gemm(device::ColMajor, device::NoTrans, device::NoTrans, 24, 9, 34, 1.0, kDivMT(2), 40, dQ(1), 40, 0.0, _tmp4, 24, offset_kDivMT, offset_dQ + 1, offset__tmp4, num_elements);
    }
    {
    unsigned offset__tmp4 = 216;
    unsigned offset_star = star_offset(2);
    unsigned offset_dQ = dQ_offset(2);
    
    
    device.gemm(device::ColMajor, device::NoTrans, device::NoTrans, 24, 9, 9, 1.0, _tmp4, 24, star(2), 9, 1.0, dQ(2), 24, offset__tmp4, offset_star, offset_dQ, num_elements);
    }
    device.api->popStackMemory();
  }
  void kernel::derivative::execute3() {
    assert(dQ(2) != nullptr);
    assert(dQ(3) != nullptr);
    assert(kDivMT(2) != nullptr);
    assert(kDivMT(0) != nullptr);
    assert(kDivMT(1) != nullptr);
    assert(star(2) != nullptr);
    assert(star(0) != nullptr);
    assert(star(1) != nullptr);
    assert(kernel::derivative::num_elements != 0);
    device::DeviceInstance& device = device::DeviceInstance::getInstance();
    float *_tmp0, *_tmp2, *_tmp4;
    float *d_buffer0 = (float*)device.api->getStackMemory(144 * num_elements * sizeof(float));
    _tmp0 = d_buffer0;
    {
    unsigned offset_kDivMT = 0;
    unsigned offset_dQ = dQ_offset(2);
    unsigned offset__tmp0 = 144;
    
    
    device.gemm(device::ColMajor, device::NoTrans, device::NoTrans, 16, 9, 17, 1.0, kDivMT(0), 40, dQ(2), 24, 0.0, _tmp0, 16, offset_kDivMT, offset_dQ + 1, offset__tmp0, num_elements);
    }
    {
    unsigned offset__tmp0 = 144;
    unsigned offset_star = star_offset(0);
    unsigned offset_dQ = dQ_offset(3);
    
    
    device.gemm(device::ColMajor, device::NoTrans, device::NoTrans, 16, 9, 9, 1.0, _tmp0, 16, star(0), 9, 0.0, dQ(3), 16, offset__tmp0, offset_star, offset_dQ, num_elements);
    }
    _tmp2 = d_buffer0;
    {
    unsigned offset_kDivMT = 0;
    unsigned offset_dQ = dQ_offset(2);
    unsigned offset__tmp2 = 144;
    
    
    device.gemm(device::ColMajor, device::NoTrans, device::NoTrans, 16, 9, 18, 1.0, kDivMT(1), 40, dQ(2), 24, 0.0, _tmp2, 16, offset_kDivMT, offset_dQ + 1, offset__tmp2, num_elements);
    }
    {
    unsigned offset__tmp2 = 144;
    unsigned offset_star = star_offset(1);
    unsigned offset_dQ = dQ_offset(3);
    
    
    device.gemm(device::ColMajor, device::NoTrans, device::NoTrans, 16, 9, 9, 1.0, _tmp2, 16, star(1), 9, 1.0, dQ(3), 16, offset__tmp2, offset_star, offset_dQ, num_elements);
    }
    _tmp4 = d_buffer0;
    {
    unsigned offset_kDivMT = 0;
    unsigned offset_dQ = dQ_offset(2);
    unsigned offset__tmp4 = 144;
    
    
    device.gemm(device::ColMajor, device::NoTrans, device::NoTrans, 16, 9, 19, 1.0, kDivMT(2), 40, dQ(2), 24, 0.0, _tmp4, 16, offset_kDivMT, offset_dQ + 1, offset__tmp4, num_elements);
    }
    {
    unsigned offset__tmp4 = 144;
    unsigned offset_star = star_offset(2);
    unsigned offset_dQ = dQ_offset(3);
    
    
    device.gemm(device::ColMajor, device::NoTrans, device::NoTrans, 16, 9, 9, 1.0, _tmp4, 16, star(2), 9, 1.0, dQ(3), 16, offset__tmp4, offset_star, offset_dQ, num_elements);
    }
    device.api->popStackMemory();
  }
  void kernel::derivative::execute4() {
    assert(dQ(3) != nullptr);
    assert(dQ(4) != nullptr);
    assert(kDivMT(2) != nullptr);
    assert(kDivMT(0) != nullptr);
    assert(kDivMT(1) != nullptr);
    assert(star(2) != nullptr);
    assert(star(0) != nullptr);
    assert(star(1) != nullptr);
    assert(kernel::derivative::num_elements != 0);
    device::DeviceInstance& device = device::DeviceInstance::getInstance();
    float *_tmp0, *_tmp2, *_tmp4;
    float *d_buffer0 = (float*)device.api->getStackMemory(72 * num_elements * sizeof(float));
    _tmp0 = d_buffer0;
    {
    unsigned offset_kDivMT = 0;
    unsigned offset_dQ = dQ_offset(3);
    unsigned offset__tmp0 = 72;
    
    
    device.gemm(device::ColMajor, device::NoTrans, device::NoTrans, 8, 9, 7, 1.0, kDivMT(0), 40, dQ(3), 16, 0.0, _tmp0, 8, offset_kDivMT, offset_dQ + 1, offset__tmp0, num_elements);
    }
    {
    unsigned offset__tmp0 = 72;
    unsigned offset_star = star_offset(0);
    unsigned offset_dQ = dQ_offset(4);
    
    
    device.gemm(device::ColMajor, device::NoTrans, device::NoTrans, 8, 9, 9, 1.0, _tmp0, 8, star(0), 9, 0.0, dQ(4), 8, offset__tmp0, offset_star, offset_dQ, num_elements);
    }
    _tmp2 = d_buffer0;
    {
    unsigned offset_kDivMT = 0;
    unsigned offset_dQ = dQ_offset(3);
    unsigned offset__tmp2 = 72;
    
    
    device.gemm(device::ColMajor, device::NoTrans, device::NoTrans, 8, 9, 8, 1.0, kDivMT(1), 40, dQ(3), 16, 0.0, _tmp2, 8, offset_kDivMT, offset_dQ + 1, offset__tmp2, num_elements);
    }
    {
    unsigned offset__tmp2 = 72;
    unsigned offset_star = star_offset(1);
    unsigned offset_dQ = dQ_offset(4);
    
    
    device.gemm(device::ColMajor, device::NoTrans, device::NoTrans, 8, 9, 9, 1.0, _tmp2, 8, star(1), 9, 1.0, dQ(4), 8, offset__tmp2, offset_star, offset_dQ, num_elements);
    }
    _tmp4 = d_buffer0;
    {
    unsigned offset_kDivMT = 0;
    unsigned offset_dQ = dQ_offset(3);
    unsigned offset__tmp4 = 72;
    
    
    device.gemm(device::ColMajor, device::NoTrans, device::NoTrans, 8, 9, 9, 1.0, kDivMT(2), 40, dQ(3), 16, 0.0, _tmp4, 8, offset_kDivMT, offset_dQ + 1, offset__tmp4, num_elements);
    }
    {
    unsigned offset__tmp4 = 72;
    unsigned offset_star = star_offset(2);
    unsigned offset_dQ = dQ_offset(4);
    
    
    device.gemm(device::ColMajor, device::NoTrans, device::NoTrans, 8, 9, 9, 1.0, _tmp4, 8, star(2), 9, 1.0, dQ(4), 8, offset__tmp4, offset_star, offset_dQ, num_elements);
    }
    device.api->popStackMemory();
  }
  void kernel::derivative::execute5() {
    assert(dQ(5) != nullptr);
    assert(dQ(4) != nullptr);
    assert(kDivMT(2) != nullptr);
    assert(kDivMT(0) != nullptr);
    assert(kDivMT(1) != nullptr);
    assert(star(2) != nullptr);
    assert(star(0) != nullptr);
    assert(star(1) != nullptr);
    assert(kernel::derivative::num_elements != 0);
    device::DeviceInstance& device = device::DeviceInstance::getInstance();
    float *_tmp0, *_tmp2, *_tmp4;
    float *d_buffer0 = (float*)device.api->getStackMemory(72 * num_elements * sizeof(float));
    _tmp0 = d_buffer0;
    {
    unsigned offset_kDivMT = 0;
    unsigned offset_dQ = dQ_offset(4);
    unsigned offset__tmp0 = 72;
    
    
    device.gemm(device::ColMajor, device::NoTrans, device::NoTrans, 8, 9, 1, 1.0, kDivMT(0), 40, dQ(4), 8, 0.0, _tmp0, 8, offset_kDivMT, offset_dQ + 1, offset__tmp0, num_elements);
    }
    {
    unsigned offset__tmp0 = 72;
    unsigned offset_star = star_offset(0);
    unsigned offset_dQ = dQ_offset(5);
    
    
    device.gemm(device::ColMajor, device::NoTrans, device::NoTrans, 8, 9, 9, 1.0, _tmp0, 8, star(0), 9, 0.0, dQ(5), 8, offset__tmp0, offset_star, offset_dQ, num_elements);
    }
    _tmp2 = d_buffer0;
    {
    unsigned offset_kDivMT = 0;
    unsigned offset_dQ = dQ_offset(4);
    unsigned offset__tmp2 = 72;
    
    
    device.gemm(device::ColMajor, device::NoTrans, device::NoTrans, 8, 9, 2, 1.0, kDivMT(1), 40, dQ(4), 8, 0.0, _tmp2, 8, offset_kDivMT, offset_dQ + 1, offset__tmp2, num_elements);
    }
    {
    unsigned offset__tmp2 = 72;
    unsigned offset_star = star_offset(1);
    unsigned offset_dQ = dQ_offset(5);
    
    
    device.gemm(device::ColMajor, device::NoTrans, device::NoTrans, 8, 9, 9, 1.0, _tmp2, 8, star(1), 9, 1.0, dQ(5), 8, offset__tmp2, offset_star, offset_dQ, num_elements);
    }
    _tmp4 = d_buffer0;
    {
    unsigned offset_kDivMT = 0;
    unsigned offset_dQ = dQ_offset(4);
    unsigned offset__tmp4 = 72;
    
    
    device.gemm(device::ColMajor, device::NoTrans, device::NoTrans, 8, 9, 3, 1.0, kDivMT(2), 40, dQ(4), 8, 0.0, _tmp4, 8, offset_kDivMT, offset_dQ + 1, offset__tmp4, num_elements);
    }
    {
    unsigned offset__tmp4 = 72;
    unsigned offset_star = star_offset(2);
    unsigned offset_dQ = dQ_offset(5);
    
    
    device.gemm(device::ColMajor, device::NoTrans, device::NoTrans, 8, 9, 9, 1.0, _tmp4, 8, star(2), 9, 1.0, dQ(5), 8, offset__tmp4, offset_star, offset_dQ, num_elements);
    }
    device.api->popStackMemory();
  }
  constexpr unsigned long const kernel::godunovState::NonZeroFlops[];
  constexpr unsigned long const kernel::godunovState::HardwareFlops[];
  constexpr kernel::godunovState::member_function_ptr kernel::godunovState::ExecutePtrs[];
  void kernel::godunovState::execute0() {
    assert(Q != nullptr);
    assert(V3mTo2n(0,0) != nullptr);
    assert(godunovMatrix != nullptr);
    assert(godunovState != nullptr);
    assert(kernel::godunovState::num_elements != 0);
    device::DeviceInstance& device = device::DeviceInstance::getInstance();
    float *_tmp0;
    float *d_buffer0 = (float*)device.api->getStackMemory(504 * num_elements * sizeof(float));
    _tmp0 = d_buffer0;
    {
    unsigned offset_V3mTo2n = 0;
    unsigned offset_Q = Q_offset;
    unsigned offset__tmp0 = 504;
    
    
    device.gemm(device::ColMajor, device::NoTrans, device::NoTrans, 56, 9, 56, 1.0, V3mTo2n(0,0), 56, Q, 56, 0.0, _tmp0, 56, offset_V3mTo2n, offset_Q, offset__tmp0, num_elements);
    }
    {
    unsigned offset__tmp0 = 504;
    unsigned offset_godunovMatrix = godunovMatrix_offset;
    unsigned offset_godunovState = godunovState_offset;
    
    
    device.gemm(device::ColMajor, device::NoTrans, device::NoTrans, 56, 9, 9, 1.0, _tmp0, 56, godunovMatrix, 9, 0.0, godunovState, 56, offset__tmp0, offset_godunovMatrix, offset_godunovState, num_elements);
    }
    device.api->popStackMemory();
  }
  void kernel::godunovState::execute1() {
    assert(Q != nullptr);
    assert(V3mTo2n(1,0) != nullptr);
    assert(godunovMatrix != nullptr);
    assert(godunovState != nullptr);
    assert(kernel::godunovState::num_elements != 0);
    device::DeviceInstance& device = device::DeviceInstance::getInstance();
    float *_tmp0;
    float *d_buffer0 = (float*)device.api->getStackMemory(504 * num_elements * sizeof(float));
    _tmp0 = d_buffer0;
    {
    unsigned offset_V3mTo2n = 0;
    unsigned offset_Q = Q_offset;
    unsigned offset__tmp0 = 504;
    
    
    device.gemm(device::ColMajor, device::NoTrans, device::NoTrans, 56, 9, 56, 1.0, V3mTo2n(1,0), 56, Q, 56, 0.0, _tmp0, 56, offset_V3mTo2n, offset_Q, offset__tmp0, num_elements);
    }
    {
    unsigned offset__tmp0 = 504;
    unsigned offset_godunovMatrix = godunovMatrix_offset;
    unsigned offset_godunovState = godunovState_offset;
    
    
    device.gemm(device::ColMajor, device::NoTrans, device::NoTrans, 56, 9, 9, 1.0, _tmp0, 56, godunovMatrix, 9, 0.0, godunovState, 56, offset__tmp0, offset_godunovMatrix, offset_godunovState, num_elements);
    }
    device.api->popStackMemory();
  }
  void kernel::godunovState::execute2() {
    assert(Q != nullptr);
    assert(V3mTo2n(2,0) != nullptr);
    assert(godunovMatrix != nullptr);
    assert(godunovState != nullptr);
    assert(kernel::godunovState::num_elements != 0);
    device::DeviceInstance& device = device::DeviceInstance::getInstance();
    float *_tmp0;
    float *d_buffer0 = (float*)device.api->getStackMemory(504 * num_elements * sizeof(float));
    _tmp0 = d_buffer0;
    {
    unsigned offset_V3mTo2n = 0;
    unsigned offset_Q = Q_offset;
    unsigned offset__tmp0 = 504;
    
    
    device.gemm(device::ColMajor, device::NoTrans, device::NoTrans, 56, 9, 56, 1.0, V3mTo2n(2,0), 56, Q, 56, 0.0, _tmp0, 56, offset_V3mTo2n, offset_Q, offset__tmp0, num_elements);
    }
    {
    unsigned offset__tmp0 = 504;
    unsigned offset_godunovMatrix = godunovMatrix_offset;
    unsigned offset_godunovState = godunovState_offset;
    
    
    device.gemm(device::ColMajor, device::NoTrans, device::NoTrans, 56, 9, 9, 1.0, _tmp0, 56, godunovMatrix, 9, 0.0, godunovState, 56, offset__tmp0, offset_godunovMatrix, offset_godunovState, num_elements);
    }
    device.api->popStackMemory();
  }
  void kernel::godunovState::execute3() {
    assert(Q != nullptr);
    assert(V3mTo2n(3,0) != nullptr);
    assert(godunovMatrix != nullptr);
    assert(godunovState != nullptr);
    assert(kernel::godunovState::num_elements != 0);
    device::DeviceInstance& device = device::DeviceInstance::getInstance();
    float *_tmp0;
    float *d_buffer0 = (float*)device.api->getStackMemory(504 * num_elements * sizeof(float));
    _tmp0 = d_buffer0;
    {
    unsigned offset_V3mTo2n = 0;
    unsigned offset_Q = Q_offset;
    unsigned offset__tmp0 = 504;
    
    
    device.gemm(device::ColMajor, device::NoTrans, device::NoTrans, 56, 9, 56, 1.0, V3mTo2n(3,0), 56, Q, 56, 0.0, _tmp0, 56, offset_V3mTo2n, offset_Q, offset__tmp0, num_elements);
    }
    {
    unsigned offset__tmp0 = 504;
    unsigned offset_godunovMatrix = godunovMatrix_offset;
    unsigned offset_godunovState = godunovState_offset;
    
    
    device.gemm(device::ColMajor, device::NoTrans, device::NoTrans, 56, 9, 9, 1.0, _tmp0, 56, godunovMatrix, 9, 0.0, godunovState, 56, offset__tmp0, offset_godunovMatrix, offset_godunovState, num_elements);
    }
    device.api->popStackMemory();
  }
  void kernel::godunovState::execute4() {
    assert(Q != nullptr);
    assert(V3mTo2n(0,1) != nullptr);
    assert(godunovMatrix != nullptr);
    assert(godunovState != nullptr);
    assert(kernel::godunovState::num_elements != 0);
    device::DeviceInstance& device = device::DeviceInstance::getInstance();
    float *_tmp0;
    float *d_buffer0 = (float*)device.api->getStackMemory(504 * num_elements * sizeof(float));
    _tmp0 = d_buffer0;
    {
    unsigned offset_V3mTo2n = 0;
    unsigned offset_Q = Q_offset;
    unsigned offset__tmp0 = 504;
    
    
    device.gemm(device::ColMajor, device::NoTrans, device::NoTrans, 56, 9, 56, 1.0, V3mTo2n(0,1), 56, Q, 56, 0.0, _tmp0, 56, offset_V3mTo2n, offset_Q, offset__tmp0, num_elements);
    }
    {
    unsigned offset__tmp0 = 504;
    unsigned offset_godunovMatrix = godunovMatrix_offset;
    unsigned offset_godunovState = godunovState_offset;
    
    
    device.gemm(device::ColMajor, device::NoTrans, device::NoTrans, 56, 9, 9, 1.0, _tmp0, 56, godunovMatrix, 9, 1.0, godunovState, 56, offset__tmp0, offset_godunovMatrix, offset_godunovState, num_elements);
    }
    device.api->popStackMemory();
  }
  void kernel::godunovState::execute5() {
    assert(Q != nullptr);
    assert(V3mTo2n(1,1) != nullptr);
    assert(godunovMatrix != nullptr);
    assert(godunovState != nullptr);
    assert(kernel::godunovState::num_elements != 0);
    device::DeviceInstance& device = device::DeviceInstance::getInstance();
    float *_tmp0;
    float *d_buffer0 = (float*)device.api->getStackMemory(504 * num_elements * sizeof(float));
    _tmp0 = d_buffer0;
    {
    unsigned offset_V3mTo2n = 0;
    unsigned offset_Q = Q_offset;
    unsigned offset__tmp0 = 504;
    
    
    device.gemm(device::ColMajor, device::NoTrans, device::NoTrans, 56, 9, 56, 1.0, V3mTo2n(1,1), 56, Q, 56, 0.0, _tmp0, 56, offset_V3mTo2n, offset_Q, offset__tmp0, num_elements);
    }
    {
    unsigned offset__tmp0 = 504;
    unsigned offset_godunovMatrix = godunovMatrix_offset;
    unsigned offset_godunovState = godunovState_offset;
    
    
    device.gemm(device::ColMajor, device::NoTrans, device::NoTrans, 56, 9, 9, 1.0, _tmp0, 56, godunovMatrix, 9, 1.0, godunovState, 56, offset__tmp0, offset_godunovMatrix, offset_godunovState, num_elements);
    }
    device.api->popStackMemory();
  }
  void kernel::godunovState::execute6() {
    assert(Q != nullptr);
    assert(V3mTo2n(2,1) != nullptr);
    assert(godunovMatrix != nullptr);
    assert(godunovState != nullptr);
    assert(kernel::godunovState::num_elements != 0);
    device::DeviceInstance& device = device::DeviceInstance::getInstance();
    float *_tmp0;
    float *d_buffer0 = (float*)device.api->getStackMemory(504 * num_elements * sizeof(float));
    _tmp0 = d_buffer0;
    {
    unsigned offset_V3mTo2n = 0;
    unsigned offset_Q = Q_offset;
    unsigned offset__tmp0 = 504;
    
    
    device.gemm(device::ColMajor, device::NoTrans, device::NoTrans, 56, 9, 56, 1.0, V3mTo2n(2,1), 56, Q, 56, 0.0, _tmp0, 56, offset_V3mTo2n, offset_Q, offset__tmp0, num_elements);
    }
    {
    unsigned offset__tmp0 = 504;
    unsigned offset_godunovMatrix = godunovMatrix_offset;
    unsigned offset_godunovState = godunovState_offset;
    
    
    device.gemm(device::ColMajor, device::NoTrans, device::NoTrans, 56, 9, 9, 1.0, _tmp0, 56, godunovMatrix, 9, 1.0, godunovState, 56, offset__tmp0, offset_godunovMatrix, offset_godunovState, num_elements);
    }
    device.api->popStackMemory();
  }
  void kernel::godunovState::execute7() {
    assert(Q != nullptr);
    assert(V3mTo2n(3,1) != nullptr);
    assert(godunovMatrix != nullptr);
    assert(godunovState != nullptr);
    assert(kernel::godunovState::num_elements != 0);
    device::DeviceInstance& device = device::DeviceInstance::getInstance();
    float *_tmp0;
    float *d_buffer0 = (float*)device.api->getStackMemory(504 * num_elements * sizeof(float));
    _tmp0 = d_buffer0;
    {
    unsigned offset_V3mTo2n = 0;
    unsigned offset_Q = Q_offset;
    unsigned offset__tmp0 = 504;
    
    
    device.gemm(device::ColMajor, device::NoTrans, device::NoTrans, 56, 9, 56, 1.0, V3mTo2n(3,1), 56, Q, 56, 0.0, _tmp0, 56, offset_V3mTo2n, offset_Q, offset__tmp0, num_elements);
    }
    {
    unsigned offset__tmp0 = 504;
    unsigned offset_godunovMatrix = godunovMatrix_offset;
    unsigned offset_godunovState = godunovState_offset;
    
    
    device.gemm(device::ColMajor, device::NoTrans, device::NoTrans, 56, 9, 9, 1.0, _tmp0, 56, godunovMatrix, 9, 1.0, godunovState, 56, offset__tmp0, offset_godunovMatrix, offset_godunovState, num_elements);
    }
    device.api->popStackMemory();
  }
  void kernel::godunovState::execute8() {
    assert(Q != nullptr);
    assert(V3mTo2n(0,2) != nullptr);
    assert(godunovMatrix != nullptr);
    assert(godunovState != nullptr);
    assert(kernel::godunovState::num_elements != 0);
    device::DeviceInstance& device = device::DeviceInstance::getInstance();
    float *_tmp0;
    float *d_buffer0 = (float*)device.api->getStackMemory(504 * num_elements * sizeof(float));
    _tmp0 = d_buffer0;
    {
    unsigned offset_V3mTo2n = 0;
    unsigned offset_Q = Q_offset;
    unsigned offset__tmp0 = 504;
    
    
    device.gemm(device::ColMajor, device::NoTrans, device::NoTrans, 56, 9, 56, 1.0, V3mTo2n(0,2), 56, Q, 56, 0.0, _tmp0, 56, offset_V3mTo2n, offset_Q, offset__tmp0, num_elements);
    }
    {
    unsigned offset__tmp0 = 504;
    unsigned offset_godunovMatrix = godunovMatrix_offset;
    unsigned offset_godunovState = godunovState_offset;
    
    
    device.gemm(device::ColMajor, device::NoTrans, device::NoTrans, 56, 9, 9, 1.0, _tmp0, 56, godunovMatrix, 9, 1.0, godunovState, 56, offset__tmp0, offset_godunovMatrix, offset_godunovState, num_elements);
    }
    device.api->popStackMemory();
  }
  void kernel::godunovState::execute9() {
    assert(Q != nullptr);
    assert(V3mTo2n(1,2) != nullptr);
    assert(godunovMatrix != nullptr);
    assert(godunovState != nullptr);
    assert(kernel::godunovState::num_elements != 0);
    device::DeviceInstance& device = device::DeviceInstance::getInstance();
    float *_tmp0;
    float *d_buffer0 = (float*)device.api->getStackMemory(504 * num_elements * sizeof(float));
    _tmp0 = d_buffer0;
    {
    unsigned offset_V3mTo2n = 0;
    unsigned offset_Q = Q_offset;
    unsigned offset__tmp0 = 504;
    
    
    device.gemm(device::ColMajor, device::NoTrans, device::NoTrans, 56, 9, 56, 1.0, V3mTo2n(1,2), 56, Q, 56, 0.0, _tmp0, 56, offset_V3mTo2n, offset_Q, offset__tmp0, num_elements);
    }
    {
    unsigned offset__tmp0 = 504;
    unsigned offset_godunovMatrix = godunovMatrix_offset;
    unsigned offset_godunovState = godunovState_offset;
    
    
    device.gemm(device::ColMajor, device::NoTrans, device::NoTrans, 56, 9, 9, 1.0, _tmp0, 56, godunovMatrix, 9, 1.0, godunovState, 56, offset__tmp0, offset_godunovMatrix, offset_godunovState, num_elements);
    }
    device.api->popStackMemory();
  }
  void kernel::godunovState::execute10() {
    assert(Q != nullptr);
    assert(V3mTo2n(2,2) != nullptr);
    assert(godunovMatrix != nullptr);
    assert(godunovState != nullptr);
    assert(kernel::godunovState::num_elements != 0);
    device::DeviceInstance& device = device::DeviceInstance::getInstance();
    float *_tmp0;
    float *d_buffer0 = (float*)device.api->getStackMemory(504 * num_elements * sizeof(float));
    _tmp0 = d_buffer0;
    {
    unsigned offset_V3mTo2n = 0;
    unsigned offset_Q = Q_offset;
    unsigned offset__tmp0 = 504;
    
    
    device.gemm(device::ColMajor, device::NoTrans, device::NoTrans, 56, 9, 56, 1.0, V3mTo2n(2,2), 56, Q, 56, 0.0, _tmp0, 56, offset_V3mTo2n, offset_Q, offset__tmp0, num_elements);
    }
    {
    unsigned offset__tmp0 = 504;
    unsigned offset_godunovMatrix = godunovMatrix_offset;
    unsigned offset_godunovState = godunovState_offset;
    
    
    device.gemm(device::ColMajor, device::NoTrans, device::NoTrans, 56, 9, 9, 1.0, _tmp0, 56, godunovMatrix, 9, 1.0, godunovState, 56, offset__tmp0, offset_godunovMatrix, offset_godunovState, num_elements);
    }
    device.api->popStackMemory();
  }
  void kernel::godunovState::execute11() {
    assert(Q != nullptr);
    assert(V3mTo2n(3,2) != nullptr);
    assert(godunovMatrix != nullptr);
    assert(godunovState != nullptr);
    assert(kernel::godunovState::num_elements != 0);
    device::DeviceInstance& device = device::DeviceInstance::getInstance();
    float *_tmp0;
    float *d_buffer0 = (float*)device.api->getStackMemory(504 * num_elements * sizeof(float));
    _tmp0 = d_buffer0;
    {
    unsigned offset_V3mTo2n = 0;
    unsigned offset_Q = Q_offset;
    unsigned offset__tmp0 = 504;
    
    
    device.gemm(device::ColMajor, device::NoTrans, device::NoTrans, 56, 9, 56, 1.0, V3mTo2n(3,2), 56, Q, 56, 0.0, _tmp0, 56, offset_V3mTo2n, offset_Q, offset__tmp0, num_elements);
    }
    {
    unsigned offset__tmp0 = 504;
    unsigned offset_godunovMatrix = godunovMatrix_offset;
    unsigned offset_godunovState = godunovState_offset;
    
    
    device.gemm(device::ColMajor, device::NoTrans, device::NoTrans, 56, 9, 9, 1.0, _tmp0, 56, godunovMatrix, 9, 1.0, godunovState, 56, offset__tmp0, offset_godunovMatrix, offset_godunovState, num_elements);
    }
    device.api->popStackMemory();
  }
  void kernel::godunovState::execute12() {
    assert(Q != nullptr);
    assert(V3mTo2n(0,3) != nullptr);
    assert(godunovMatrix != nullptr);
    assert(godunovState != nullptr);
    assert(kernel::godunovState::num_elements != 0);
    device::DeviceInstance& device = device::DeviceInstance::getInstance();
    float *_tmp0;
    float *d_buffer0 = (float*)device.api->getStackMemory(504 * num_elements * sizeof(float));
    _tmp0 = d_buffer0;
    {
    unsigned offset_V3mTo2n = 0;
    unsigned offset_Q = Q_offset;
    unsigned offset__tmp0 = 504;
    
    
    device.gemm(device::ColMajor, device::NoTrans, device::NoTrans, 56, 9, 56, 1.0, V3mTo2n(0,3), 56, Q, 56, 0.0, _tmp0, 56, offset_V3mTo2n, offset_Q, offset__tmp0, num_elements);
    }
    {
    unsigned offset__tmp0 = 504;
    unsigned offset_godunovMatrix = godunovMatrix_offset;
    unsigned offset_godunovState = godunovState_offset;
    
    
    device.gemm(device::ColMajor, device::NoTrans, device::NoTrans, 56, 9, 9, 1.0, _tmp0, 56, godunovMatrix, 9, 1.0, godunovState, 56, offset__tmp0, offset_godunovMatrix, offset_godunovState, num_elements);
    }
    device.api->popStackMemory();
  }
  void kernel::godunovState::execute13() {
    assert(Q != nullptr);
    assert(V3mTo2n(1,3) != nullptr);
    assert(godunovMatrix != nullptr);
    assert(godunovState != nullptr);
    assert(kernel::godunovState::num_elements != 0);
    device::DeviceInstance& device = device::DeviceInstance::getInstance();
    float *_tmp0;
    float *d_buffer0 = (float*)device.api->getStackMemory(504 * num_elements * sizeof(float));
    _tmp0 = d_buffer0;
    {
    unsigned offset_V3mTo2n = 0;
    unsigned offset_Q = Q_offset;
    unsigned offset__tmp0 = 504;
    
    
    device.gemm(device::ColMajor, device::NoTrans, device::NoTrans, 56, 9, 56, 1.0, V3mTo2n(1,3), 56, Q, 56, 0.0, _tmp0, 56, offset_V3mTo2n, offset_Q, offset__tmp0, num_elements);
    }
    {
    unsigned offset__tmp0 = 504;
    unsigned offset_godunovMatrix = godunovMatrix_offset;
    unsigned offset_godunovState = godunovState_offset;
    
    
    device.gemm(device::ColMajor, device::NoTrans, device::NoTrans, 56, 9, 9, 1.0, _tmp0, 56, godunovMatrix, 9, 1.0, godunovState, 56, offset__tmp0, offset_godunovMatrix, offset_godunovState, num_elements);
    }
    device.api->popStackMemory();
  }
  void kernel::godunovState::execute14() {
    assert(Q != nullptr);
    assert(V3mTo2n(2,3) != nullptr);
    assert(godunovMatrix != nullptr);
    assert(godunovState != nullptr);
    assert(kernel::godunovState::num_elements != 0);
    device::DeviceInstance& device = device::DeviceInstance::getInstance();
    float *_tmp0;
    float *d_buffer0 = (float*)device.api->getStackMemory(504 * num_elements * sizeof(float));
    _tmp0 = d_buffer0;
    {
    unsigned offset_V3mTo2n = 0;
    unsigned offset_Q = Q_offset;
    unsigned offset__tmp0 = 504;
    
    
    device.gemm(device::ColMajor, device::NoTrans, device::NoTrans, 56, 9, 56, 1.0, V3mTo2n(2,3), 56, Q, 56, 0.0, _tmp0, 56, offset_V3mTo2n, offset_Q, offset__tmp0, num_elements);
    }
    {
    unsigned offset__tmp0 = 504;
    unsigned offset_godunovMatrix = godunovMatrix_offset;
    unsigned offset_godunovState = godunovState_offset;
    
    
    device.gemm(device::ColMajor, device::NoTrans, device::NoTrans, 56, 9, 9, 1.0, _tmp0, 56, godunovMatrix, 9, 1.0, godunovState, 56, offset__tmp0, offset_godunovMatrix, offset_godunovState, num_elements);
    }
    device.api->popStackMemory();
  }
  void kernel::godunovState::execute15() {
    assert(Q != nullptr);
    assert(V3mTo2n(3,3) != nullptr);
    assert(godunovMatrix != nullptr);
    assert(godunovState != nullptr);
    assert(kernel::godunovState::num_elements != 0);
    device::DeviceInstance& device = device::DeviceInstance::getInstance();
    float *_tmp0;
    float *d_buffer0 = (float*)device.api->getStackMemory(504 * num_elements * sizeof(float));
    _tmp0 = d_buffer0;
    {
    unsigned offset_V3mTo2n = 0;
    unsigned offset_Q = Q_offset;
    unsigned offset__tmp0 = 504;
    
    
    device.gemm(device::ColMajor, device::NoTrans, device::NoTrans, 56, 9, 56, 1.0, V3mTo2n(3,3), 56, Q, 56, 0.0, _tmp0, 56, offset_V3mTo2n, offset_Q, offset__tmp0, num_elements);
    }
    {
    unsigned offset__tmp0 = 504;
    unsigned offset_godunovMatrix = godunovMatrix_offset;
    unsigned offset_godunovState = godunovState_offset;
    
    
    device.gemm(device::ColMajor, device::NoTrans, device::NoTrans, 56, 9, 9, 1.0, _tmp0, 56, godunovMatrix, 9, 1.0, godunovState, 56, offset__tmp0, offset_godunovMatrix, offset_godunovState, num_elements);
    }
    device.api->popStackMemory();
  }
  constexpr unsigned long const kernel::nodalFlux::NonZeroFlops[];
  constexpr unsigned long const kernel::nodalFlux::HardwareFlops[];
  constexpr kernel::nodalFlux::member_function_ptr kernel::nodalFlux::ExecutePtrs[];
  void kernel::nodalFlux::execute0() {
    assert(Q != nullptr);
    assert(V3mTo2nTWDivM(0,0) != nullptr);
    assert(fluxSolver != nullptr);
    assert(godunovState != nullptr);
    assert(kernel::nodalFlux::num_elements != 0);
    device::DeviceInstance& device = device::DeviceInstance::getInstance();
    float *_tmp0;
    float *d_buffer0 = (float*)device.api->getStackMemory(504 * num_elements * sizeof(float));
    _tmp0 = d_buffer0;
    {
    unsigned offset_godunovState = godunovState_offset;
    unsigned offset_fluxSolver = fluxSolver_offset;
    unsigned offset__tmp0 = 504;
    
    
    device.gemm(device::ColMajor, device::NoTrans, device::NoTrans, 56, 9, 9, 1.0, godunovState, 56, fluxSolver, 9, 0.0, _tmp0, 56, offset_godunovState, offset_fluxSolver, offset__tmp0, num_elements);
    }
    {
    unsigned offset_V3mTo2nTWDivM = 0;
    unsigned offset__tmp0 = 504;
    unsigned offset_Q = Q_offset;
    
    
    device.gemm(device::ColMajor, device::NoTrans, device::NoTrans, 56, 9, 49, 1.0, V3mTo2nTWDivM(0,0), 56, _tmp0, 56, 1.0, Q, 56, offset_V3mTo2nTWDivM, offset__tmp0, offset_Q, num_elements);
    }
    device.api->popStackMemory();
  }
  void kernel::nodalFlux::execute1() {
    assert(Q != nullptr);
    assert(V3mTo2nTWDivM(1,0) != nullptr);
    assert(fluxSolver != nullptr);
    assert(godunovState != nullptr);
    assert(kernel::nodalFlux::num_elements != 0);
    device::DeviceInstance& device = device::DeviceInstance::getInstance();
    float *_tmp0;
    float *d_buffer0 = (float*)device.api->getStackMemory(504 * num_elements * sizeof(float));
    _tmp0 = d_buffer0;
    {
    unsigned offset_godunovState = godunovState_offset;
    unsigned offset_fluxSolver = fluxSolver_offset;
    unsigned offset__tmp0 = 504;
    
    
    device.gemm(device::ColMajor, device::NoTrans, device::NoTrans, 56, 9, 9, 1.0, godunovState, 56, fluxSolver, 9, 0.0, _tmp0, 56, offset_godunovState, offset_fluxSolver, offset__tmp0, num_elements);
    }
    {
    unsigned offset_V3mTo2nTWDivM = 0;
    unsigned offset__tmp0 = 504;
    unsigned offset_Q = Q_offset;
    
    
    device.gemm(device::ColMajor, device::NoTrans, device::NoTrans, 56, 9, 49, 1.0, V3mTo2nTWDivM(1,0), 56, _tmp0, 56, 1.0, Q, 56, offset_V3mTo2nTWDivM, offset__tmp0, offset_Q, num_elements);
    }
    device.api->popStackMemory();
  }
  void kernel::nodalFlux::execute2() {
    assert(Q != nullptr);
    assert(V3mTo2nTWDivM(2,0) != nullptr);
    assert(fluxSolver != nullptr);
    assert(godunovState != nullptr);
    assert(kernel::nodalFlux::num_elements != 0);
    device::DeviceInstance& device = device::DeviceInstance::getInstance();
    float *_tmp0;
    float *d_buffer0 = (float*)device.api->getStackMemory(504 * num_elements * sizeof(float));
    _tmp0 = d_buffer0;
    {
    unsigned offset_godunovState = godunovState_offset;
    unsigned offset_fluxSolver = fluxSolver_offset;
    unsigned offset__tmp0 = 504;
    
    
    device.gemm(device::ColMajor, device::NoTrans, device::NoTrans, 56, 9, 9, 1.0, godunovState, 56, fluxSolver, 9, 0.0, _tmp0, 56, offset_godunovState, offset_fluxSolver, offset__tmp0, num_elements);
    }
    {
    unsigned offset_V3mTo2nTWDivM = 0;
    unsigned offset__tmp0 = 504;
    unsigned offset_Q = Q_offset;
    
    
    device.gemm(device::ColMajor, device::NoTrans, device::NoTrans, 56, 9, 49, 1.0, V3mTo2nTWDivM(2,0), 56, _tmp0, 56, 1.0, Q, 56, offset_V3mTo2nTWDivM, offset__tmp0, offset_Q, num_elements);
    }
    device.api->popStackMemory();
  }
  void kernel::nodalFlux::execute3() {
    assert(Q != nullptr);
    assert(V3mTo2nTWDivM(3,0) != nullptr);
    assert(fluxSolver != nullptr);
    assert(godunovState != nullptr);
    assert(kernel::nodalFlux::num_elements != 0);
    device::DeviceInstance& device = device::DeviceInstance::getInstance();
    float *_tmp0;
    float *d_buffer0 = (float*)device.api->getStackMemory(504 * num_elements * sizeof(float));
    _tmp0 = d_buffer0;
    {
    unsigned offset_godunovState = godunovState_offset;
    unsigned offset_fluxSolver = fluxSolver_offset;
    unsigned offset__tmp0 = 504;
    
    
    device.gemm(device::ColMajor, device::NoTrans, device::NoTrans, 56, 9, 9, 1.0, godunovState, 56, fluxSolver, 9, 0.0, _tmp0, 56, offset_godunovState, offset_fluxSolver, offset__tmp0, num_elements);
    }
    {
    unsigned offset_V3mTo2nTWDivM = 0;
    unsigned offset__tmp0 = 504;
    unsigned offset_Q = Q_offset;
    
    
    device.gemm(device::ColMajor, device::NoTrans, device::NoTrans, 56, 9, 49, 1.0, V3mTo2nTWDivM(3,0), 56, _tmp0, 56, 1.0, Q, 56, offset_V3mTo2nTWDivM, offset__tmp0, offset_Q, num_elements);
    }
    device.api->popStackMemory();
  }
  void kernel::nodalFlux::execute4() {
    assert(Q != nullptr);
    assert(V3mTo2nTWDivM(0,1) != nullptr);
    assert(fluxSolver != nullptr);
    assert(godunovState != nullptr);
    assert(kernel::nodalFlux::num_elements != 0);
    device::DeviceInstance& device = device::DeviceInstance::getInstance();
    float *_tmp0;
    float *d_buffer0 = (float*)device.api->getStackMemory(504 * num_elements * sizeof(float));
    _tmp0 = d_buffer0;
    {
    unsigned offset_godunovState = godunovState_offset;
    unsigned offset_fluxSolver = fluxSolver_offset;
    unsigned offset__tmp0 = 504;
    
    
    device.gemm(device::ColMajor, device::NoTrans, device::NoTrans, 56, 9, 9, 1.0, godunovState, 56, fluxSolver, 9, 0.0, _tmp0, 56, offset_godunovState, offset_fluxSolver, offset__tmp0, num_elements);
    }
    {
    unsigned offset_V3mTo2nTWDivM = 0;
    unsigned offset__tmp0 = 504;
    unsigned offset_Q = Q_offset;
    
    
    device.gemm(device::ColMajor, device::NoTrans, device::NoTrans, 56, 9, 49, 1.0, V3mTo2nTWDivM(0,1), 56, _tmp0, 56, 1.0, Q, 56, offset_V3mTo2nTWDivM, offset__tmp0, offset_Q, num_elements);
    }
    device.api->popStackMemory();
  }
  void kernel::nodalFlux::execute5() {
    assert(Q != nullptr);
    assert(V3mTo2nTWDivM(1,1) != nullptr);
    assert(fluxSolver != nullptr);
    assert(godunovState != nullptr);
    assert(kernel::nodalFlux::num_elements != 0);
    device::DeviceInstance& device = device::DeviceInstance::getInstance();
    float *_tmp0;
    float *d_buffer0 = (float*)device.api->getStackMemory(504 * num_elements * sizeof(float));
    _tmp0 = d_buffer0;
    {
    unsigned offset_godunovState = godunovState_offset;
    unsigned offset_fluxSolver = fluxSolver_offset;
    unsigned offset__tmp0 = 504;
    
    
    device.gemm(device::ColMajor, device::NoTrans, device::NoTrans, 56, 9, 9, 1.0, godunovState, 56, fluxSolver, 9, 0.0, _tmp0, 56, offset_godunovState, offset_fluxSolver, offset__tmp0, num_elements);
    }
    {
    unsigned offset_V3mTo2nTWDivM = 0;
    unsigned offset__tmp0 = 504;
    unsigned offset_Q = Q_offset;
    
    
    device.gemm(device::ColMajor, device::NoTrans, device::NoTrans, 56, 9, 49, 1.0, V3mTo2nTWDivM(1,1), 56, _tmp0, 56, 1.0, Q, 56, offset_V3mTo2nTWDivM, offset__tmp0, offset_Q, num_elements);
    }
    device.api->popStackMemory();
  }
  void kernel::nodalFlux::execute6() {
    assert(Q != nullptr);
    assert(V3mTo2nTWDivM(2,1) != nullptr);
    assert(fluxSolver != nullptr);
    assert(godunovState != nullptr);
    assert(kernel::nodalFlux::num_elements != 0);
    device::DeviceInstance& device = device::DeviceInstance::getInstance();
    float *_tmp0;
    float *d_buffer0 = (float*)device.api->getStackMemory(504 * num_elements * sizeof(float));
    _tmp0 = d_buffer0;
    {
    unsigned offset_godunovState = godunovState_offset;
    unsigned offset_fluxSolver = fluxSolver_offset;
    unsigned offset__tmp0 = 504;
    
    
    device.gemm(device::ColMajor, device::NoTrans, device::NoTrans, 56, 9, 9, 1.0, godunovState, 56, fluxSolver, 9, 0.0, _tmp0, 56, offset_godunovState, offset_fluxSolver, offset__tmp0, num_elements);
    }
    {
    unsigned offset_V3mTo2nTWDivM = 0;
    unsigned offset__tmp0 = 504;
    unsigned offset_Q = Q_offset;
    
    
    device.gemm(device::ColMajor, device::NoTrans, device::NoTrans, 56, 9, 49, 1.0, V3mTo2nTWDivM(2,1), 56, _tmp0, 56, 1.0, Q, 56, offset_V3mTo2nTWDivM, offset__tmp0, offset_Q, num_elements);
    }
    device.api->popStackMemory();
  }
  void kernel::nodalFlux::execute7() {
    assert(Q != nullptr);
    assert(V3mTo2nTWDivM(3,1) != nullptr);
    assert(fluxSolver != nullptr);
    assert(godunovState != nullptr);
    assert(kernel::nodalFlux::num_elements != 0);
    device::DeviceInstance& device = device::DeviceInstance::getInstance();
    float *_tmp0;
    float *d_buffer0 = (float*)device.api->getStackMemory(504 * num_elements * sizeof(float));
    _tmp0 = d_buffer0;
    {
    unsigned offset_godunovState = godunovState_offset;
    unsigned offset_fluxSolver = fluxSolver_offset;
    unsigned offset__tmp0 = 504;
    
    
    device.gemm(device::ColMajor, device::NoTrans, device::NoTrans, 56, 9, 9, 1.0, godunovState, 56, fluxSolver, 9, 0.0, _tmp0, 56, offset_godunovState, offset_fluxSolver, offset__tmp0, num_elements);
    }
    {
    unsigned offset_V3mTo2nTWDivM = 0;
    unsigned offset__tmp0 = 504;
    unsigned offset_Q = Q_offset;
    
    
    device.gemm(device::ColMajor, device::NoTrans, device::NoTrans, 56, 9, 49, 1.0, V3mTo2nTWDivM(3,1), 56, _tmp0, 56, 1.0, Q, 56, offset_V3mTo2nTWDivM, offset__tmp0, offset_Q, num_elements);
    }
    device.api->popStackMemory();
  }
  void kernel::nodalFlux::execute8() {
    assert(Q != nullptr);
    assert(V3mTo2nTWDivM(0,2) != nullptr);
    assert(fluxSolver != nullptr);
    assert(godunovState != nullptr);
    assert(kernel::nodalFlux::num_elements != 0);
    device::DeviceInstance& device = device::DeviceInstance::getInstance();
    float *_tmp0;
    float *d_buffer0 = (float*)device.api->getStackMemory(504 * num_elements * sizeof(float));
    _tmp0 = d_buffer0;
    {
    unsigned offset_godunovState = godunovState_offset;
    unsigned offset_fluxSolver = fluxSolver_offset;
    unsigned offset__tmp0 = 504;
    
    
    device.gemm(device::ColMajor, device::NoTrans, device::NoTrans, 56, 9, 9, 1.0, godunovState, 56, fluxSolver, 9, 0.0, _tmp0, 56, offset_godunovState, offset_fluxSolver, offset__tmp0, num_elements);
    }
    {
    unsigned offset_V3mTo2nTWDivM = 0;
    unsigned offset__tmp0 = 504;
    unsigned offset_Q = Q_offset;
    
    
    device.gemm(device::ColMajor, device::NoTrans, device::NoTrans, 56, 9, 49, 1.0, V3mTo2nTWDivM(0,2), 56, _tmp0, 56, 1.0, Q, 56, offset_V3mTo2nTWDivM, offset__tmp0, offset_Q, num_elements);
    }
    device.api->popStackMemory();
  }
  void kernel::nodalFlux::execute9() {
    assert(Q != nullptr);
    assert(V3mTo2nTWDivM(1,2) != nullptr);
    assert(fluxSolver != nullptr);
    assert(godunovState != nullptr);
    assert(kernel::nodalFlux::num_elements != 0);
    device::DeviceInstance& device = device::DeviceInstance::getInstance();
    float *_tmp0;
    float *d_buffer0 = (float*)device.api->getStackMemory(504 * num_elements * sizeof(float));
    _tmp0 = d_buffer0;
    {
    unsigned offset_godunovState = godunovState_offset;
    unsigned offset_fluxSolver = fluxSolver_offset;
    unsigned offset__tmp0 = 504;
    
    
    device.gemm(device::ColMajor, device::NoTrans, device::NoTrans, 56, 9, 9, 1.0, godunovState, 56, fluxSolver, 9, 0.0, _tmp0, 56, offset_godunovState, offset_fluxSolver, offset__tmp0, num_elements);
    }
    {
    unsigned offset_V3mTo2nTWDivM = 0;
    unsigned offset__tmp0 = 504;
    unsigned offset_Q = Q_offset;
    
    
    device.gemm(device::ColMajor, device::NoTrans, device::NoTrans, 56, 9, 49, 1.0, V3mTo2nTWDivM(1,2), 56, _tmp0, 56, 1.0, Q, 56, offset_V3mTo2nTWDivM, offset__tmp0, offset_Q, num_elements);
    }
    device.api->popStackMemory();
  }
  void kernel::nodalFlux::execute10() {
    assert(Q != nullptr);
    assert(V3mTo2nTWDivM(2,2) != nullptr);
    assert(fluxSolver != nullptr);
    assert(godunovState != nullptr);
    assert(kernel::nodalFlux::num_elements != 0);
    device::DeviceInstance& device = device::DeviceInstance::getInstance();
    float *_tmp0;
    float *d_buffer0 = (float*)device.api->getStackMemory(504 * num_elements * sizeof(float));
    _tmp0 = d_buffer0;
    {
    unsigned offset_godunovState = godunovState_offset;
    unsigned offset_fluxSolver = fluxSolver_offset;
    unsigned offset__tmp0 = 504;
    
    
    device.gemm(device::ColMajor, device::NoTrans, device::NoTrans, 56, 9, 9, 1.0, godunovState, 56, fluxSolver, 9, 0.0, _tmp0, 56, offset_godunovState, offset_fluxSolver, offset__tmp0, num_elements);
    }
    {
    unsigned offset_V3mTo2nTWDivM = 0;
    unsigned offset__tmp0 = 504;
    unsigned offset_Q = Q_offset;
    
    
    device.gemm(device::ColMajor, device::NoTrans, device::NoTrans, 56, 9, 49, 1.0, V3mTo2nTWDivM(2,2), 56, _tmp0, 56, 1.0, Q, 56, offset_V3mTo2nTWDivM, offset__tmp0, offset_Q, num_elements);
    }
    device.api->popStackMemory();
  }
  void kernel::nodalFlux::execute11() {
    assert(Q != nullptr);
    assert(V3mTo2nTWDivM(3,2) != nullptr);
    assert(fluxSolver != nullptr);
    assert(godunovState != nullptr);
    assert(kernel::nodalFlux::num_elements != 0);
    device::DeviceInstance& device = device::DeviceInstance::getInstance();
    float *_tmp0;
    float *d_buffer0 = (float*)device.api->getStackMemory(504 * num_elements * sizeof(float));
    _tmp0 = d_buffer0;
    {
    unsigned offset_godunovState = godunovState_offset;
    unsigned offset_fluxSolver = fluxSolver_offset;
    unsigned offset__tmp0 = 504;
    
    
    device.gemm(device::ColMajor, device::NoTrans, device::NoTrans, 56, 9, 9, 1.0, godunovState, 56, fluxSolver, 9, 0.0, _tmp0, 56, offset_godunovState, offset_fluxSolver, offset__tmp0, num_elements);
    }
    {
    unsigned offset_V3mTo2nTWDivM = 0;
    unsigned offset__tmp0 = 504;
    unsigned offset_Q = Q_offset;
    
    
    device.gemm(device::ColMajor, device::NoTrans, device::NoTrans, 56, 9, 49, 1.0, V3mTo2nTWDivM(3,2), 56, _tmp0, 56, 1.0, Q, 56, offset_V3mTo2nTWDivM, offset__tmp0, offset_Q, num_elements);
    }
    device.api->popStackMemory();
  }
  void kernel::nodalFlux::execute12() {
    assert(Q != nullptr);
    assert(V3mTo2nTWDivM(0,3) != nullptr);
    assert(fluxSolver != nullptr);
    assert(godunovState != nullptr);
    assert(kernel::nodalFlux::num_elements != 0);
    device::DeviceInstance& device = device::DeviceInstance::getInstance();
    float *_tmp0;
    float *d_buffer0 = (float*)device.api->getStackMemory(504 * num_elements * sizeof(float));
    _tmp0 = d_buffer0;
    {
    unsigned offset_godunovState = godunovState_offset;
    unsigned offset_fluxSolver = fluxSolver_offset;
    unsigned offset__tmp0 = 504;
    
    
    device.gemm(device::ColMajor, device::NoTrans, device::NoTrans, 56, 9, 9, 1.0, godunovState, 56, fluxSolver, 9, 0.0, _tmp0, 56, offset_godunovState, offset_fluxSolver, offset__tmp0, num_elements);
    }
    {
    unsigned offset_V3mTo2nTWDivM = 0;
    unsigned offset__tmp0 = 504;
    unsigned offset_Q = Q_offset;
    
    
    device.gemm(device::ColMajor, device::NoTrans, device::NoTrans, 56, 9, 49, 1.0, V3mTo2nTWDivM(0,3), 56, _tmp0, 56, 1.0, Q, 56, offset_V3mTo2nTWDivM, offset__tmp0, offset_Q, num_elements);
    }
    device.api->popStackMemory();
  }
  void kernel::nodalFlux::execute13() {
    assert(Q != nullptr);
    assert(V3mTo2nTWDivM(1,3) != nullptr);
    assert(fluxSolver != nullptr);
    assert(godunovState != nullptr);
    assert(kernel::nodalFlux::num_elements != 0);
    device::DeviceInstance& device = device::DeviceInstance::getInstance();
    float *_tmp0;
    float *d_buffer0 = (float*)device.api->getStackMemory(504 * num_elements * sizeof(float));
    _tmp0 = d_buffer0;
    {
    unsigned offset_godunovState = godunovState_offset;
    unsigned offset_fluxSolver = fluxSolver_offset;
    unsigned offset__tmp0 = 504;
    
    
    device.gemm(device::ColMajor, device::NoTrans, device::NoTrans, 56, 9, 9, 1.0, godunovState, 56, fluxSolver, 9, 0.0, _tmp0, 56, offset_godunovState, offset_fluxSolver, offset__tmp0, num_elements);
    }
    {
    unsigned offset_V3mTo2nTWDivM = 0;
    unsigned offset__tmp0 = 504;
    unsigned offset_Q = Q_offset;
    
    
    device.gemm(device::ColMajor, device::NoTrans, device::NoTrans, 56, 9, 49, 1.0, V3mTo2nTWDivM(1,3), 56, _tmp0, 56, 1.0, Q, 56, offset_V3mTo2nTWDivM, offset__tmp0, offset_Q, num_elements);
    }
    device.api->popStackMemory();
  }
  void kernel::nodalFlux::execute14() {
    assert(Q != nullptr);
    assert(V3mTo2nTWDivM(2,3) != nullptr);
    assert(fluxSolver != nullptr);
    assert(godunovState != nullptr);
    assert(kernel::nodalFlux::num_elements != 0);
    device::DeviceInstance& device = device::DeviceInstance::getInstance();
    float *_tmp0;
    float *d_buffer0 = (float*)device.api->getStackMemory(504 * num_elements * sizeof(float));
    _tmp0 = d_buffer0;
    {
    unsigned offset_godunovState = godunovState_offset;
    unsigned offset_fluxSolver = fluxSolver_offset;
    unsigned offset__tmp0 = 504;
    
    
    device.gemm(device::ColMajor, device::NoTrans, device::NoTrans, 56, 9, 9, 1.0, godunovState, 56, fluxSolver, 9, 0.0, _tmp0, 56, offset_godunovState, offset_fluxSolver, offset__tmp0, num_elements);
    }
    {
    unsigned offset_V3mTo2nTWDivM = 0;
    unsigned offset__tmp0 = 504;
    unsigned offset_Q = Q_offset;
    
    
    device.gemm(device::ColMajor, device::NoTrans, device::NoTrans, 56, 9, 49, 1.0, V3mTo2nTWDivM(2,3), 56, _tmp0, 56, 1.0, Q, 56, offset_V3mTo2nTWDivM, offset__tmp0, offset_Q, num_elements);
    }
    device.api->popStackMemory();
  }
  void kernel::nodalFlux::execute15() {
    assert(Q != nullptr);
    assert(V3mTo2nTWDivM(3,3) != nullptr);
    assert(fluxSolver != nullptr);
    assert(godunovState != nullptr);
    assert(kernel::nodalFlux::num_elements != 0);
    device::DeviceInstance& device = device::DeviceInstance::getInstance();
    float *_tmp0;
    float *d_buffer0 = (float*)device.api->getStackMemory(504 * num_elements * sizeof(float));
    _tmp0 = d_buffer0;
    {
    unsigned offset_godunovState = godunovState_offset;
    unsigned offset_fluxSolver = fluxSolver_offset;
    unsigned offset__tmp0 = 504;
    
    
    device.gemm(device::ColMajor, device::NoTrans, device::NoTrans, 56, 9, 9, 1.0, godunovState, 56, fluxSolver, 9, 0.0, _tmp0, 56, offset_godunovState, offset_fluxSolver, offset__tmp0, num_elements);
    }
    {
    unsigned offset_V3mTo2nTWDivM = 0;
    unsigned offset__tmp0 = 504;
    unsigned offset_Q = Q_offset;
    
    
    device.gemm(device::ColMajor, device::NoTrans, device::NoTrans, 56, 9, 49, 1.0, V3mTo2nTWDivM(3,3), 56, _tmp0, 56, 1.0, Q, 56, offset_V3mTo2nTWDivM, offset__tmp0, offset_Q, num_elements);
    }
    device.api->popStackMemory();
  }
}
