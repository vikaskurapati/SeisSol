#ifdef ACL_DEVICE
#include "device_utils.h"
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
    // get device temporary memory menager
    DeviceTemporaryMemoryMenager &tmp_manager = DeviceTemporaryMemoryMenager::get_instance();
    double *_tmp0, *_tmp2, *_tmp4;
    double *d_buffer0 = (double*)tmp_manager.get_mem(324 * num_elements * sizeof(double));
    _tmp0 = d_buffer0;
    {
    unsigned offsets_I = I_offset;
    unsigned offsets_star = star_offset(0);
    unsigned offsets__tmp0 = 324;
    
    
    device_gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 36, 9, 9, 1.0, I, 56, star(0), 9, 0.0, _tmp0, 36, offsets_I, offsets_star, offsets__tmp0, num_elements);
    }
    {
    unsigned offsets_kDivM = 0;
    unsigned offsets__tmp0 = 324;
    unsigned offsets_Q = Q_offset;
    
    
    device_gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 56, 9, 35, 1.0, kDivM(0), 56, _tmp0, 36, 1.0, Q, 56, offsets_kDivM, offsets__tmp0, offsets_Q, num_elements);
    }
    _tmp2 = d_buffer0;
    {
    unsigned offsets_I = I_offset;
    unsigned offsets_star = star_offset(1);
    unsigned offsets__tmp2 = 324;
    
    
    device_gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 36, 9, 9, 1.0, I, 56, star(1), 9, 0.0, _tmp2, 36, offsets_I, offsets_star, offsets__tmp2, num_elements);
    }
    {
    unsigned offsets_kDivM = 0;
    unsigned offsets__tmp2 = 324;
    unsigned offsets_Q = Q_offset;
    
    
    device_gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 56, 9, 35, 1.0, kDivM(1), 56, _tmp2, 36, 1.0, Q, 56, offsets_kDivM, offsets__tmp2, offsets_Q, num_elements);
    }
    _tmp4 = d_buffer0;
    {
    unsigned offsets_I = I_offset;
    unsigned offsets_star = star_offset(2);
    unsigned offsets__tmp4 = 324;
    
    
    device_gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 36, 9, 9, 1.0, I, 56, star(2), 9, 0.0, _tmp4, 36, offsets_I, offsets_star, offsets__tmp4, num_elements);
    }
    {
    unsigned offsets_kDivM = 0;
    unsigned offsets__tmp4 = 324;
    unsigned offsets_Q = Q_offset;
    
    
    device_gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 56, 9, 35, 1.0, kDivM(2), 56, _tmp4, 36, 1.0, Q, 56, offsets_kDivM, offsets__tmp4, offsets_Q, num_elements);
    }
    tmp_manager.free();
  }
  constexpr unsigned long const kernel::rotateGodunovStateLocal::NonZeroFlops;
  constexpr unsigned long const kernel::rotateGodunovStateLocal::HardwareFlops;
  void kernel::rotateGodunovStateLocal::execute() {
    assert(QgodLocal != nullptr);
    assert(Tinv != nullptr);
    assert(godunovMatrix != nullptr);
    assert(kernel::rotateGodunovStateLocal::num_elements != 0);
    {
      unsigned offsets_Tinv = Tinv_offset;
      unsigned offsets_QgodLocal = QgodLocal_offset;
      unsigned offsets_godunovMatrix = godunovMatrix_offset;


      device_gemm(CblasColMajor, CblasTrans, CblasNoTrans, 9, 9, 9, 1.0, Tinv, 9, QgodLocal, 9, 0.0, godunovMatrix, 9, offsets_Tinv, offsets_QgodLocal, offsets_godunovMatrix, num_elements);
    }
  }
  constexpr unsigned long const kernel::rotateGodunovStateNeighbor::NonZeroFlops;
  constexpr unsigned long const kernel::rotateGodunovStateNeighbor::HardwareFlops;
  void kernel::rotateGodunovStateNeighbor::execute() {
    assert(QgodNeighbor != nullptr);
    assert(Tinv != nullptr);
    assert(godunovMatrix != nullptr);
    assert(kernel::rotateGodunovStateNeighbor::num_elements != 0);
    {
      unsigned offsets_Tinv = Tinv_offset;
      unsigned offsets_QgodNeighbor = QgodNeighbor_offset;
      unsigned offsets_godunovMatrix = godunovMatrix_offset;


      device_gemm(CblasColMajor, CblasTrans, CblasNoTrans, 9, 9, 9, 1.0, Tinv, 9, QgodNeighbor, 9, 0.0, godunovMatrix, 9, offsets_Tinv, offsets_QgodNeighbor, offsets_godunovMatrix, num_elements);
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
    {
      unsigned offsets_star = star_offset(0);
      unsigned offsets_T = T_offset;
      unsigned offsets_fluxSolver = fluxSolver_offset;


      device_gemm(CblasColMajor, CblasNoTrans, CblasTrans, 9, 9, 9, fluxScale, star(0), 9, T, 9, 0.0, fluxSolver, 9, offsets_star, offsets_T, offsets_fluxSolver, num_elements);
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
    // get device temporary memory menager
    DeviceTemporaryMemoryMenager &tmp_manager = DeviceTemporaryMemoryMenager::get_instance();
    double *_tmp0, *_tmp1;
    double *d_buffer0 = (double*)tmp_manager.get_mem(216 * num_elements * sizeof(double));
    double *d_buffer1 = (double*)tmp_manager.get_mem(216 * num_elements * sizeof(double));
    _tmp0 = d_buffer0;
    {
    unsigned offsets_fMrT = 0;
    unsigned offsets_I = I_offset;
    unsigned offsets__tmp0 = 216;
    
    
    device_gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 24, 9, 56, 1.0, fMrT(0), 24, I, 56, 0.0, _tmp0, 24, offsets_fMrT, offsets_I, offsets__tmp0, num_elements);
    }
    _tmp1 = d_buffer1;
    {
    unsigned offsets__tmp0 = 216;
    unsigned offsets_AplusT = AplusT_offset;
    unsigned offsets__tmp1 = 216;
    
    
    device_gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 24, 9, 9, 1.0, _tmp0, 24, AplusT, 9, 0.0, _tmp1, 24, offsets__tmp0, offsets_AplusT, offsets__tmp1, num_elements);
    }
    {
    unsigned offsets_rDivM = 0;
    unsigned offsets__tmp1 = 216;
    unsigned offsets_Q = Q_offset;
    
    
    device_gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 56, 9, 21, 1.0, rDivM(0), 56, _tmp1, 24, 1.0, Q, 56, offsets_rDivM, offsets__tmp1, offsets_Q, num_elements);
    }
    tmp_manager.free();
    tmp_manager.free();
  }
  void kernel::localFlux::execute1() {
    assert(AplusT != nullptr);
    assert(I != nullptr);
    assert(Q != nullptr);
    assert(fMrT(1) != nullptr);
    assert(rDivM(1) != nullptr);
    assert(kernel::localFlux::num_elements != 0);
    // get device temporary memory menager
    DeviceTemporaryMemoryMenager &tmp_manager = DeviceTemporaryMemoryMenager::get_instance();
    double *_tmp0, *_tmp1;
    double *d_buffer0 = (double*)tmp_manager.get_mem(216 * num_elements * sizeof(double));
    double *d_buffer1 = (double*)tmp_manager.get_mem(216 * num_elements * sizeof(double));
    _tmp0 = d_buffer0;
    {
    unsigned offsets_fMrT = 0;
    unsigned offsets_I = I_offset;
    unsigned offsets__tmp0 = 216;
    
    
    device_gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 24, 9, 56, 1.0, fMrT(1), 24, I, 56, 0.0, _tmp0, 24, offsets_fMrT, offsets_I, offsets__tmp0, num_elements);
    }
    _tmp1 = d_buffer1;
    {
    unsigned offsets__tmp0 = 216;
    unsigned offsets_AplusT = AplusT_offset;
    unsigned offsets__tmp1 = 216;
    
    
    device_gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 24, 9, 9, 1.0, _tmp0, 24, AplusT, 9, 0.0, _tmp1, 24, offsets__tmp0, offsets_AplusT, offsets__tmp1, num_elements);
    }
    {
    unsigned offsets_rDivM = 0;
    unsigned offsets__tmp1 = 216;
    unsigned offsets_Q = Q_offset;
    
    
    device_gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 56, 9, 21, 1.0, rDivM(1), 56, _tmp1, 24, 1.0, Q, 56, offsets_rDivM, offsets__tmp1, offsets_Q, num_elements);
    }
    tmp_manager.free();
    tmp_manager.free();
  }
  void kernel::localFlux::execute2() {
    assert(AplusT != nullptr);
    assert(I != nullptr);
    assert(Q != nullptr);
    assert(fMrT(2) != nullptr);
    assert(rDivM(2) != nullptr);
    assert(kernel::localFlux::num_elements != 0);
    // get device temporary memory menager
    DeviceTemporaryMemoryMenager &tmp_manager = DeviceTemporaryMemoryMenager::get_instance();
    double *_tmp0, *_tmp1;
    double *d_buffer0 = (double*)tmp_manager.get_mem(216 * num_elements * sizeof(double));
    double *d_buffer1 = (double*)tmp_manager.get_mem(216 * num_elements * sizeof(double));
    _tmp0 = d_buffer0;
    {
    unsigned offsets_fMrT = 0;
    unsigned offsets_I = I_offset;
    unsigned offsets__tmp0 = 216;
    
    
    device_gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 24, 9, 56, 1.0, fMrT(2), 24, I, 56, 0.0, _tmp0, 24, offsets_fMrT, offsets_I, offsets__tmp0, num_elements);
    }
    _tmp1 = d_buffer1;
    {
    unsigned offsets__tmp0 = 216;
    unsigned offsets_AplusT = AplusT_offset;
    unsigned offsets__tmp1 = 216;
    
    
    device_gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 24, 9, 9, 1.0, _tmp0, 24, AplusT, 9, 0.0, _tmp1, 24, offsets__tmp0, offsets_AplusT, offsets__tmp1, num_elements);
    }
    {
    unsigned offsets_rDivM = 0;
    unsigned offsets__tmp1 = 216;
    unsigned offsets_Q = Q_offset;
    
    
    device_gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 56, 9, 21, 1.0, rDivM(2), 56, _tmp1, 24, 1.0, Q, 56, offsets_rDivM, offsets__tmp1, offsets_Q, num_elements);
    }
    tmp_manager.free();
    tmp_manager.free();
  }
  void kernel::localFlux::execute3() {
    assert(AplusT != nullptr);
    assert(I != nullptr);
    assert(Q != nullptr);
    assert(fMrT(3) != nullptr);
    assert(rDivM(3) != nullptr);
    assert(kernel::localFlux::num_elements != 0);
    // get device temporary memory menager
    DeviceTemporaryMemoryMenager &tmp_manager = DeviceTemporaryMemoryMenager::get_instance();
    double *_tmp0, *_tmp1;
    double *d_buffer0 = (double*)tmp_manager.get_mem(216 * num_elements * sizeof(double));
    double *d_buffer1 = (double*)tmp_manager.get_mem(216 * num_elements * sizeof(double));
    _tmp0 = d_buffer0;
    {
    unsigned offsets_fMrT = 0;
    unsigned offsets_I = I_offset;
    unsigned offsets__tmp0 = 216;
    
    
    device_gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 24, 9, 56, 1.0, fMrT(3), 24, I, 56, 0.0, _tmp0, 24, offsets_fMrT, offsets_I, offsets__tmp0, num_elements);
    }
    _tmp1 = d_buffer1;
    {
    unsigned offsets__tmp0 = 216;
    unsigned offsets_AplusT = AplusT_offset;
    unsigned offsets__tmp1 = 216;
    
    
    device_gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 24, 9, 9, 1.0, _tmp0, 24, AplusT, 9, 0.0, _tmp1, 24, offsets__tmp0, offsets_AplusT, offsets__tmp1, num_elements);
    }
    {
    unsigned offsets_rDivM = 0;
    unsigned offsets__tmp1 = 216;
    unsigned offsets_Q = Q_offset;
    
    
    device_gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 56, 9, 21, 1.0, rDivM(3), 56, _tmp1, 24, 1.0, Q, 56, offsets_rDivM, offsets__tmp1, offsets_Q, num_elements);
    }
    tmp_manager.free();
    tmp_manager.free();
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
    // get device temporary memory menager
    DeviceTemporaryMemoryMenager &tmp_manager = DeviceTemporaryMemoryMenager::get_instance();
    double *_tmp0, *_tmp1, *_tmp2;
    double *d_buffer0 = (double*)tmp_manager.get_mem(216 * num_elements * sizeof(double));
    double *d_buffer1 = (double*)tmp_manager.get_mem(216 * num_elements * sizeof(double));
    _tmp0 = d_buffer0;
    {
      unsigned offsets_rT = 0;
      unsigned offsets_I = I_offset;
      unsigned offsets__tmp0 = 216;


      device_gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 24, 9, 56, 1.0, rT(0), 24, I, 56, 0.0, _tmp0, 24, offsets_rT, offsets_I, offsets__tmp0, num_elements);
    }
    _tmp1 = d_buffer1;
    {
      unsigned offsets_fP = 0;
      unsigned offsets__tmp0 = 216;
      unsigned offsets__tmp1 = 216;


      device_gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 24, 9, 21, 1.0, fP(0), 24, _tmp0, 24, 0.0, _tmp1, 24, offsets_fP, offsets__tmp0, offsets__tmp1, num_elements);
    }
    _tmp2 = d_buffer0;
    {
      unsigned offsets__tmp1 = 216;
      unsigned offsets_AminusT = AminusT_offset;
      unsigned offsets__tmp2 = 216;


      device_gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 24, 9, 9, 1.0, _tmp1, 24, AminusT, 9, 0.0, _tmp2, 24, offsets__tmp1, offsets_AminusT, offsets__tmp2, num_elements);
    }
    {
      unsigned offsets_rDivM = 0;
      unsigned offsets__tmp2 = 216;
      unsigned offsets_Q = Q_offset;


      device_gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 56, 9, 21, 1.0, rDivM(0), 56, _tmp2, 24, 1.0, Q, 56, offsets_rDivM, offsets__tmp2, offsets_Q, num_elements);
    }
    tmp_manager.free();
    tmp_manager.free();
  }
  void kernel::neighboringFlux::execute1() {
    assert(AminusT != nullptr);
    assert(I != nullptr);
    assert(Q != nullptr);
    assert(fP(1) != nullptr);
    assert(rDivM(0) != nullptr);
    assert(rT(0) != nullptr);
    assert(kernel::neighboringFlux::num_elements != 0);
    // get device temporary memory menager
    DeviceTemporaryMemoryMenager &tmp_manager = DeviceTemporaryMemoryMenager::get_instance();
    double *_tmp0, *_tmp1, *_tmp2;
    double *d_buffer0 = (double*)tmp_manager.get_mem(216 * num_elements * sizeof(double));
    double *d_buffer1 = (double*)tmp_manager.get_mem(216 * num_elements * sizeof(double));
    _tmp0 = d_buffer0;
    {
      unsigned offsets_rT = 0;
      unsigned offsets_I = I_offset;
      unsigned offsets__tmp0 = 216;


      device_gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 24, 9, 56, 1.0, rT(0), 24, I, 56, 0.0, _tmp0, 24, offsets_rT, offsets_I, offsets__tmp0, num_elements);
    }
    _tmp1 = d_buffer1;
    {
      unsigned offsets_fP = 0;
      unsigned offsets__tmp0 = 216;
      unsigned offsets__tmp1 = 216;


      device_gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 24, 9, 21, 1.0, fP(1), 24, _tmp0, 24, 0.0, _tmp1, 24, offsets_fP, offsets__tmp0, offsets__tmp1, num_elements);
    }
    _tmp2 = d_buffer0;
    {
      unsigned offsets__tmp1 = 216;
      unsigned offsets_AminusT = AminusT_offset;
      unsigned offsets__tmp2 = 216;


      device_gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 24, 9, 9, 1.0, _tmp1, 24, AminusT, 9, 0.0, _tmp2, 24, offsets__tmp1, offsets_AminusT, offsets__tmp2, num_elements);
    }
    {
      unsigned offsets_rDivM = 0;
      unsigned offsets__tmp2 = 216;
      unsigned offsets_Q = Q_offset;


      device_gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 56, 9, 21, 1.0, rDivM(0), 56, _tmp2, 24, 1.0, Q, 56, offsets_rDivM, offsets__tmp2, offsets_Q, num_elements);
    }
    tmp_manager.free();
    tmp_manager.free();
  }
  void kernel::neighboringFlux::execute2() {
    assert(AminusT != nullptr);
    assert(I != nullptr);
    assert(Q != nullptr);
    assert(fP(2) != nullptr);
    assert(rDivM(0) != nullptr);
    assert(rT(0) != nullptr);
    assert(kernel::neighboringFlux::num_elements != 0);
    // get device temporary memory menager
    DeviceTemporaryMemoryMenager &tmp_manager = DeviceTemporaryMemoryMenager::get_instance();
    double *_tmp0, *_tmp1, *_tmp2;
    double *d_buffer0 = (double*)tmp_manager.get_mem(216 * num_elements * sizeof(double));
    double *d_buffer1 = (double*)tmp_manager.get_mem(216 * num_elements * sizeof(double));
    _tmp0 = d_buffer0;
    {
      unsigned offsets_rT = 0;
      unsigned offsets_I = I_offset;
      unsigned offsets__tmp0 = 216;


      device_gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 24, 9, 56, 1.0, rT(0), 24, I, 56, 0.0, _tmp0, 24, offsets_rT, offsets_I, offsets__tmp0, num_elements);
    }
    _tmp1 = d_buffer1;
    {
      unsigned offsets_fP = 0;
      unsigned offsets__tmp0 = 216;
      unsigned offsets__tmp1 = 216;


      device_gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 24, 9, 21, 1.0, fP(2), 24, _tmp0, 24, 0.0, _tmp1, 24, offsets_fP, offsets__tmp0, offsets__tmp1, num_elements);
    }
    _tmp2 = d_buffer0;
    {
      unsigned offsets__tmp1 = 216;
      unsigned offsets_AminusT = AminusT_offset;
      unsigned offsets__tmp2 = 216;


      device_gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 24, 9, 9, 1.0, _tmp1, 24, AminusT, 9, 0.0, _tmp2, 24, offsets__tmp1, offsets_AminusT, offsets__tmp2, num_elements);
    }
    {
      unsigned offsets_rDivM = 0;
      unsigned offsets__tmp2 = 216;
      unsigned offsets_Q = Q_offset;


      device_gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 56, 9, 21, 1.0, rDivM(0), 56, _tmp2, 24, 1.0, Q, 56, offsets_rDivM, offsets__tmp2, offsets_Q, num_elements);
    }
    tmp_manager.free();
    tmp_manager.free();
  }
  void kernel::neighboringFlux::execute3() {
    assert(AminusT != nullptr);
    assert(I != nullptr);
    assert(Q != nullptr);
    assert(fP(0) != nullptr);
    assert(rDivM(0) != nullptr);
    assert(rT(1) != nullptr);
    assert(kernel::neighboringFlux::num_elements != 0);
    // get device temporary memory menager
    DeviceTemporaryMemoryMenager &tmp_manager = DeviceTemporaryMemoryMenager::get_instance();
    double *_tmp0, *_tmp1, *_tmp2;
    double *d_buffer0 = (double*)tmp_manager.get_mem(216 * num_elements * sizeof(double));
    double *d_buffer1 = (double*)tmp_manager.get_mem(216 * num_elements * sizeof(double));
    _tmp0 = d_buffer0;
    {
      unsigned offsets_rT = 0;
      unsigned offsets_I = I_offset;
      unsigned offsets__tmp0 = 216;


      device_gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 24, 9, 56, 1.0, rT(1), 24, I, 56, 0.0, _tmp0, 24, offsets_rT, offsets_I, offsets__tmp0, num_elements);
    }
    _tmp1 = d_buffer1;
    {
      unsigned offsets_fP = 0;
      unsigned offsets__tmp0 = 216;
      unsigned offsets__tmp1 = 216;


      device_gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 24, 9, 21, 1.0, fP(0), 24, _tmp0, 24, 0.0, _tmp1, 24, offsets_fP, offsets__tmp0, offsets__tmp1, num_elements);
    }
    _tmp2 = d_buffer0;
    {
      unsigned offsets__tmp1 = 216;
      unsigned offsets_AminusT = AminusT_offset;
      unsigned offsets__tmp2 = 216;


      device_gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 24, 9, 9, 1.0, _tmp1, 24, AminusT, 9, 0.0, _tmp2, 24, offsets__tmp1, offsets_AminusT, offsets__tmp2, num_elements);
    }
    {
      unsigned offsets_rDivM = 0;
      unsigned offsets__tmp2 = 216;
      unsigned offsets_Q = Q_offset;


      device_gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 56, 9, 21, 1.0, rDivM(0), 56, _tmp2, 24, 1.0, Q, 56, offsets_rDivM, offsets__tmp2, offsets_Q, num_elements);
    }
    tmp_manager.free();
    tmp_manager.free();
  }
  void kernel::neighboringFlux::execute4() {
    assert(AminusT != nullptr);
    assert(I != nullptr);
    assert(Q != nullptr);
    assert(fP(1) != nullptr);
    assert(rDivM(0) != nullptr);
    assert(rT(1) != nullptr);
    assert(kernel::neighboringFlux::num_elements != 0);
    // get device temporary memory menager
    DeviceTemporaryMemoryMenager &tmp_manager = DeviceTemporaryMemoryMenager::get_instance();
    double *_tmp0, *_tmp1, *_tmp2;
    double *d_buffer0 = (double*)tmp_manager.get_mem(216 * num_elements * sizeof(double));
    double *d_buffer1 = (double*)tmp_manager.get_mem(216 * num_elements * sizeof(double));
    _tmp0 = d_buffer0;
    {
      unsigned offsets_rT = 0;
      unsigned offsets_I = I_offset;
      unsigned offsets__tmp0 = 216;


      device_gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 24, 9, 56, 1.0, rT(1), 24, I, 56, 0.0, _tmp0, 24, offsets_rT, offsets_I, offsets__tmp0, num_elements);
    }
    _tmp1 = d_buffer1;
    {
      unsigned offsets_fP = 0;
      unsigned offsets__tmp0 = 216;
      unsigned offsets__tmp1 = 216;


      device_gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 24, 9, 21, 1.0, fP(1), 24, _tmp0, 24, 0.0, _tmp1, 24, offsets_fP, offsets__tmp0, offsets__tmp1, num_elements);
    }
    _tmp2 = d_buffer0;
    {
      unsigned offsets__tmp1 = 216;
      unsigned offsets_AminusT = AminusT_offset;
      unsigned offsets__tmp2 = 216;


      device_gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 24, 9, 9, 1.0, _tmp1, 24, AminusT, 9, 0.0, _tmp2, 24, offsets__tmp1, offsets_AminusT, offsets__tmp2, num_elements);
    }
    {
      unsigned offsets_rDivM = 0;
      unsigned offsets__tmp2 = 216;
      unsigned offsets_Q = Q_offset;


      device_gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 56, 9, 21, 1.0, rDivM(0), 56, _tmp2, 24, 1.0, Q, 56, offsets_rDivM, offsets__tmp2, offsets_Q, num_elements);
    }
    tmp_manager.free();
    tmp_manager.free();
  }
  void kernel::neighboringFlux::execute5() {
    assert(AminusT != nullptr);
    assert(I != nullptr);
    assert(Q != nullptr);
    assert(fP(2) != nullptr);
    assert(rDivM(0) != nullptr);
    assert(rT(1) != nullptr);
    assert(kernel::neighboringFlux::num_elements != 0);
    // get device temporary memory menager
    DeviceTemporaryMemoryMenager &tmp_manager = DeviceTemporaryMemoryMenager::get_instance();
    double *_tmp0, *_tmp1, *_tmp2;
    double *d_buffer0 = (double*)tmp_manager.get_mem(216 * num_elements * sizeof(double));
    double *d_buffer1 = (double*)tmp_manager.get_mem(216 * num_elements * sizeof(double));
    _tmp0 = d_buffer0;
    {
      unsigned offsets_rT = 0;
      unsigned offsets_I = I_offset;
      unsigned offsets__tmp0 = 216;


      device_gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 24, 9, 56, 1.0, rT(1), 24, I, 56, 0.0, _tmp0, 24, offsets_rT, offsets_I, offsets__tmp0, num_elements);
    }
    _tmp1 = d_buffer1;
    {
      unsigned offsets_fP = 0;
      unsigned offsets__tmp0 = 216;
      unsigned offsets__tmp1 = 216;


      device_gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 24, 9, 21, 1.0, fP(2), 24, _tmp0, 24, 0.0, _tmp1, 24, offsets_fP, offsets__tmp0, offsets__tmp1, num_elements);
    }
    _tmp2 = d_buffer0;
    {
      unsigned offsets__tmp1 = 216;
      unsigned offsets_AminusT = AminusT_offset;
      unsigned offsets__tmp2 = 216;


      device_gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 24, 9, 9, 1.0, _tmp1, 24, AminusT, 9, 0.0, _tmp2, 24, offsets__tmp1, offsets_AminusT, offsets__tmp2, num_elements);
    }
    {
      unsigned offsets_rDivM = 0;
      unsigned offsets__tmp2 = 216;
      unsigned offsets_Q = Q_offset;


      device_gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 56, 9, 21, 1.0, rDivM(0), 56, _tmp2, 24, 1.0, Q, 56, offsets_rDivM, offsets__tmp2, offsets_Q, num_elements);
    }
    tmp_manager.free();
    tmp_manager.free();
  }
  void kernel::neighboringFlux::execute6() {
    assert(AminusT != nullptr);
    assert(I != nullptr);
    assert(Q != nullptr);
    assert(fP(0) != nullptr);
    assert(rDivM(0) != nullptr);
    assert(rT(2) != nullptr);
    assert(kernel::neighboringFlux::num_elements != 0);
    // get device temporary memory menager
    DeviceTemporaryMemoryMenager &tmp_manager = DeviceTemporaryMemoryMenager::get_instance();
    double *_tmp0, *_tmp1, *_tmp2;
    double *d_buffer0 = (double*)tmp_manager.get_mem(216 * num_elements * sizeof(double));
    double *d_buffer1 = (double*)tmp_manager.get_mem(216 * num_elements * sizeof(double));
    _tmp0 = d_buffer0;
    {
      unsigned offsets_rT = 0;
      unsigned offsets_I = I_offset;
      unsigned offsets__tmp0 = 216;


      device_gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 24, 9, 56, 1.0, rT(2), 24, I, 56, 0.0, _tmp0, 24, offsets_rT, offsets_I, offsets__tmp0, num_elements);
    }
    _tmp1 = d_buffer1;
    {
      unsigned offsets_fP = 0;
      unsigned offsets__tmp0 = 216;
      unsigned offsets__tmp1 = 216;


      device_gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 24, 9, 21, 1.0, fP(0), 24, _tmp0, 24, 0.0, _tmp1, 24, offsets_fP, offsets__tmp0, offsets__tmp1, num_elements);
    }
    _tmp2 = d_buffer0;
    {
      unsigned offsets__tmp1 = 216;
      unsigned offsets_AminusT = AminusT_offset;
      unsigned offsets__tmp2 = 216;


      device_gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 24, 9, 9, 1.0, _tmp1, 24, AminusT, 9, 0.0, _tmp2, 24, offsets__tmp1, offsets_AminusT, offsets__tmp2, num_elements);
    }
    {
      unsigned offsets_rDivM = 0;
      unsigned offsets__tmp2 = 216;
      unsigned offsets_Q = Q_offset;


      device_gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 56, 9, 21, 1.0, rDivM(0), 56, _tmp2, 24, 1.0, Q, 56, offsets_rDivM, offsets__tmp2, offsets_Q, num_elements);
    }
    tmp_manager.free();
    tmp_manager.free();
  }
  void kernel::neighboringFlux::execute7() {
    assert(AminusT != nullptr);
    assert(I != nullptr);
    assert(Q != nullptr);
    assert(fP(1) != nullptr);
    assert(rDivM(0) != nullptr);
    assert(rT(2) != nullptr);
    assert(kernel::neighboringFlux::num_elements != 0);
    // get device temporary memory menager
    DeviceTemporaryMemoryMenager &tmp_manager = DeviceTemporaryMemoryMenager::get_instance();
    double *_tmp0, *_tmp1, *_tmp2;
    double *d_buffer0 = (double*)tmp_manager.get_mem(216 * num_elements * sizeof(double));
    double *d_buffer1 = (double*)tmp_manager.get_mem(216 * num_elements * sizeof(double));
    _tmp0 = d_buffer0;
    {
      unsigned offsets_rT = 0;
      unsigned offsets_I = I_offset;
      unsigned offsets__tmp0 = 216;


      device_gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 24, 9, 56, 1.0, rT(2), 24, I, 56, 0.0, _tmp0, 24, offsets_rT, offsets_I, offsets__tmp0, num_elements);
    }
    _tmp1 = d_buffer1;
    {
      unsigned offsets_fP = 0;
      unsigned offsets__tmp0 = 216;
      unsigned offsets__tmp1 = 216;


      device_gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 24, 9, 21, 1.0, fP(1), 24, _tmp0, 24, 0.0, _tmp1, 24, offsets_fP, offsets__tmp0, offsets__tmp1, num_elements);
    }
    _tmp2 = d_buffer0;
    {
      unsigned offsets__tmp1 = 216;
      unsigned offsets_AminusT = AminusT_offset;
      unsigned offsets__tmp2 = 216;


      device_gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 24, 9, 9, 1.0, _tmp1, 24, AminusT, 9, 0.0, _tmp2, 24, offsets__tmp1, offsets_AminusT, offsets__tmp2, num_elements);
    }
    {
      unsigned offsets_rDivM = 0;
      unsigned offsets__tmp2 = 216;
      unsigned offsets_Q = Q_offset;


      device_gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 56, 9, 21, 1.0, rDivM(0), 56, _tmp2, 24, 1.0, Q, 56, offsets_rDivM, offsets__tmp2, offsets_Q, num_elements);
    }
    tmp_manager.free();
    tmp_manager.free();
  }
  void kernel::neighboringFlux::execute8() {
    assert(AminusT != nullptr);
    assert(I != nullptr);
    assert(Q != nullptr);
    assert(fP(2) != nullptr);
    assert(rDivM(0) != nullptr);
    assert(rT(2) != nullptr);
    assert(kernel::neighboringFlux::num_elements != 0);
    // get device temporary memory menager
    DeviceTemporaryMemoryMenager &tmp_manager = DeviceTemporaryMemoryMenager::get_instance();
    double *_tmp0, *_tmp1, *_tmp2;
    double *d_buffer0 = (double*)tmp_manager.get_mem(216 * num_elements * sizeof(double));
    double *d_buffer1 = (double*)tmp_manager.get_mem(216 * num_elements * sizeof(double));
    _tmp0 = d_buffer0;
    {
      unsigned offsets_rT = 0;
      unsigned offsets_I = I_offset;
      unsigned offsets__tmp0 = 216;


      device_gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 24, 9, 56, 1.0, rT(2), 24, I, 56, 0.0, _tmp0, 24, offsets_rT, offsets_I, offsets__tmp0, num_elements);
    }
    _tmp1 = d_buffer1;
    {
      unsigned offsets_fP = 0;
      unsigned offsets__tmp0 = 216;
      unsigned offsets__tmp1 = 216;


      device_gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 24, 9, 21, 1.0, fP(2), 24, _tmp0, 24, 0.0, _tmp1, 24, offsets_fP, offsets__tmp0, offsets__tmp1, num_elements);
    }
    _tmp2 = d_buffer0;
    {
      unsigned offsets__tmp1 = 216;
      unsigned offsets_AminusT = AminusT_offset;
      unsigned offsets__tmp2 = 216;


      device_gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 24, 9, 9, 1.0, _tmp1, 24, AminusT, 9, 0.0, _tmp2, 24, offsets__tmp1, offsets_AminusT, offsets__tmp2, num_elements);
    }
    {
      unsigned offsets_rDivM = 0;
      unsigned offsets__tmp2 = 216;
      unsigned offsets_Q = Q_offset;


      device_gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 56, 9, 21, 1.0, rDivM(0), 56, _tmp2, 24, 1.0, Q, 56, offsets_rDivM, offsets__tmp2, offsets_Q, num_elements);
    }
    tmp_manager.free();
    tmp_manager.free();
  }
  void kernel::neighboringFlux::execute9() {
    assert(AminusT != nullptr);
    assert(I != nullptr);
    assert(Q != nullptr);
    assert(fP(0) != nullptr);
    assert(rDivM(0) != nullptr);
    assert(rT(3) != nullptr);
    assert(kernel::neighboringFlux::num_elements != 0);
    // get device temporary memory menager
    DeviceTemporaryMemoryMenager &tmp_manager = DeviceTemporaryMemoryMenager::get_instance();
    double *_tmp0, *_tmp1, *_tmp2;
    double *d_buffer0 = (double*)tmp_manager.get_mem(216 * num_elements * sizeof(double));
    double *d_buffer1 = (double*)tmp_manager.get_mem(216 * num_elements * sizeof(double));
    _tmp0 = d_buffer0;
    {
      unsigned offsets_rT = 0;
      unsigned offsets_I = I_offset;
      unsigned offsets__tmp0 = 216;


      device_gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 24, 9, 56, 1.0, rT(3), 24, I, 56, 0.0, _tmp0, 24, offsets_rT, offsets_I, offsets__tmp0, num_elements);
    }
    _tmp1 = d_buffer1;
    {
      unsigned offsets_fP = 0;
      unsigned offsets__tmp0 = 216;
      unsigned offsets__tmp1 = 216;


      device_gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 24, 9, 21, 1.0, fP(0), 24, _tmp0, 24, 0.0, _tmp1, 24, offsets_fP, offsets__tmp0, offsets__tmp1, num_elements);
    }
    _tmp2 = d_buffer0;
    {
      unsigned offsets__tmp1 = 216;
      unsigned offsets_AminusT = AminusT_offset;
      unsigned offsets__tmp2 = 216;


      device_gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 24, 9, 9, 1.0, _tmp1, 24, AminusT, 9, 0.0, _tmp2, 24, offsets__tmp1, offsets_AminusT, offsets__tmp2, num_elements);
    }
    {
      unsigned offsets_rDivM = 0;
      unsigned offsets__tmp2 = 216;
      unsigned offsets_Q = Q_offset;


      device_gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 56, 9, 21, 1.0, rDivM(0), 56, _tmp2, 24, 1.0, Q, 56, offsets_rDivM, offsets__tmp2, offsets_Q, num_elements);
    }
    tmp_manager.free();
    tmp_manager.free();
  }
  void kernel::neighboringFlux::execute10() {
    assert(AminusT != nullptr);
    assert(I != nullptr);
    assert(Q != nullptr);
    assert(fP(1) != nullptr);
    assert(rDivM(0) != nullptr);
    assert(rT(3) != nullptr);
    assert(kernel::neighboringFlux::num_elements != 0);
    // get device temporary memory menager
    DeviceTemporaryMemoryMenager &tmp_manager = DeviceTemporaryMemoryMenager::get_instance();
    double *_tmp0, *_tmp1, *_tmp2;
    double *d_buffer0 = (double*)tmp_manager.get_mem(216 * num_elements * sizeof(double));
    double *d_buffer1 = (double*)tmp_manager.get_mem(216 * num_elements * sizeof(double));
    _tmp0 = d_buffer0;
    {
      unsigned offsets_rT = 0;
      unsigned offsets_I = I_offset;
      unsigned offsets__tmp0 = 216;


      device_gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 24, 9, 56, 1.0, rT(3), 24, I, 56, 0.0, _tmp0, 24, offsets_rT, offsets_I, offsets__tmp0, num_elements);
    }
    _tmp1 = d_buffer1;
    {
      unsigned offsets_fP = 0;
      unsigned offsets__tmp0 = 216;
      unsigned offsets__tmp1 = 216;


      device_gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 24, 9, 21, 1.0, fP(1), 24, _tmp0, 24, 0.0, _tmp1, 24, offsets_fP, offsets__tmp0, offsets__tmp1, num_elements);
    }
    _tmp2 = d_buffer0;
    {
      unsigned offsets__tmp1 = 216;
      unsigned offsets_AminusT = AminusT_offset;
      unsigned offsets__tmp2 = 216;


      device_gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 24, 9, 9, 1.0, _tmp1, 24, AminusT, 9, 0.0, _tmp2, 24, offsets__tmp1, offsets_AminusT, offsets__tmp2, num_elements);
    }
    {
      unsigned offsets_rDivM = 0;
      unsigned offsets__tmp2 = 216;
      unsigned offsets_Q = Q_offset;


      device_gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 56, 9, 21, 1.0, rDivM(0), 56, _tmp2, 24, 1.0, Q, 56, offsets_rDivM, offsets__tmp2, offsets_Q, num_elements);
    }
    tmp_manager.free();
    tmp_manager.free();
  }
  void kernel::neighboringFlux::execute11() {
    assert(AminusT != nullptr);
    assert(I != nullptr);
    assert(Q != nullptr);
    assert(fP(2) != nullptr);
    assert(rDivM(0) != nullptr);
    assert(rT(3) != nullptr);
    assert(kernel::neighboringFlux::num_elements != 0);
    // get device temporary memory menager
    DeviceTemporaryMemoryMenager &tmp_manager = DeviceTemporaryMemoryMenager::get_instance();
    double *_tmp0, *_tmp1, *_tmp2;
    double *d_buffer0 = (double*)tmp_manager.get_mem(216 * num_elements * sizeof(double));
    double *d_buffer1 = (double*)tmp_manager.get_mem(216 * num_elements * sizeof(double));
    _tmp0 = d_buffer0;
    {
      unsigned offsets_rT = 0;
      unsigned offsets_I = I_offset;
      unsigned offsets__tmp0 = 216;


      device_gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 24, 9, 56, 1.0, rT(3), 24, I, 56, 0.0, _tmp0, 24, offsets_rT, offsets_I, offsets__tmp0, num_elements);
    }
    _tmp1 = d_buffer1;
    {
      unsigned offsets_fP = 0;
      unsigned offsets__tmp0 = 216;
      unsigned offsets__tmp1 = 216;


      device_gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 24, 9, 21, 1.0, fP(2), 24, _tmp0, 24, 0.0, _tmp1, 24, offsets_fP, offsets__tmp0, offsets__tmp1, num_elements);
    }
    _tmp2 = d_buffer0;
    {
      unsigned offsets__tmp1 = 216;
      unsigned offsets_AminusT = AminusT_offset;
      unsigned offsets__tmp2 = 216;


      device_gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 24, 9, 9, 1.0, _tmp1, 24, AminusT, 9, 0.0, _tmp2, 24, offsets__tmp1, offsets_AminusT, offsets__tmp2, num_elements);
    }
    {
      unsigned offsets_rDivM = 0;
      unsigned offsets__tmp2 = 216;
      unsigned offsets_Q = Q_offset;


      device_gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 56, 9, 21, 1.0, rDivM(0), 56, _tmp2, 24, 1.0, Q, 56, offsets_rDivM, offsets__tmp2, offsets_Q, num_elements);
    }
    tmp_manager.free();
    tmp_manager.free();
  }
  void kernel::neighboringFlux::execute12() {
    assert(AminusT != nullptr);
    assert(I != nullptr);
    assert(Q != nullptr);
    assert(fP(0) != nullptr);
    assert(rDivM(1) != nullptr);
    assert(rT(0) != nullptr);
    assert(kernel::neighboringFlux::num_elements != 0);
    // get device temporary memory menager
    DeviceTemporaryMemoryMenager &tmp_manager = DeviceTemporaryMemoryMenager::get_instance();
    double *_tmp0, *_tmp1, *_tmp2;
    double *d_buffer0 = (double*)tmp_manager.get_mem(216 * num_elements * sizeof(double));
    double *d_buffer1 = (double*)tmp_manager.get_mem(216 * num_elements * sizeof(double));
    _tmp0 = d_buffer0;
    {
      unsigned offsets_rT = 0;
      unsigned offsets_I = I_offset;
      unsigned offsets__tmp0 = 216;


      device_gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 24, 9, 56, 1.0, rT(0), 24, I, 56, 0.0, _tmp0, 24, offsets_rT, offsets_I, offsets__tmp0, num_elements);
    }
    _tmp1 = d_buffer1;
    {
      unsigned offsets_fP = 0;
      unsigned offsets__tmp0 = 216;
      unsigned offsets__tmp1 = 216;


      device_gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 24, 9, 21, 1.0, fP(0), 24, _tmp0, 24, 0.0, _tmp1, 24, offsets_fP, offsets__tmp0, offsets__tmp1, num_elements);
    }
    _tmp2 = d_buffer0;
    {
      unsigned offsets__tmp1 = 216;
      unsigned offsets_AminusT = AminusT_offset;
      unsigned offsets__tmp2 = 216;


      device_gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 24, 9, 9, 1.0, _tmp1, 24, AminusT, 9, 0.0, _tmp2, 24, offsets__tmp1, offsets_AminusT, offsets__tmp2, num_elements);
    }
    {
      unsigned offsets_rDivM = 0;
      unsigned offsets__tmp2 = 216;
      unsigned offsets_Q = Q_offset;


      device_gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 56, 9, 21, 1.0, rDivM(1), 56, _tmp2, 24, 1.0, Q, 56, offsets_rDivM, offsets__tmp2, offsets_Q, num_elements);
    }
    tmp_manager.free();
    tmp_manager.free();
  }
  void kernel::neighboringFlux::execute13() {
    assert(AminusT != nullptr);
    assert(I != nullptr);
    assert(Q != nullptr);
    assert(fP(1) != nullptr);
    assert(rDivM(1) != nullptr);
    assert(rT(0) != nullptr);
    assert(kernel::neighboringFlux::num_elements != 0);
    // get device temporary memory menager
    DeviceTemporaryMemoryMenager &tmp_manager = DeviceTemporaryMemoryMenager::get_instance();
    double *_tmp0, *_tmp1, *_tmp2;
    double *d_buffer0 = (double*)tmp_manager.get_mem(216 * num_elements * sizeof(double));
    double *d_buffer1 = (double*)tmp_manager.get_mem(216 * num_elements * sizeof(double));
    _tmp0 = d_buffer0;
    {
      unsigned offsets_rT = 0;
      unsigned offsets_I = I_offset;
      unsigned offsets__tmp0 = 216;


      device_gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 24, 9, 56, 1.0, rT(0), 24, I, 56, 0.0, _tmp0, 24, offsets_rT, offsets_I, offsets__tmp0, num_elements);
    }
    _tmp1 = d_buffer1;
    {
      unsigned offsets_fP = 0;
      unsigned offsets__tmp0 = 216;
      unsigned offsets__tmp1 = 216;


      device_gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 24, 9, 21, 1.0, fP(1), 24, _tmp0, 24, 0.0, _tmp1, 24, offsets_fP, offsets__tmp0, offsets__tmp1, num_elements);
    }
    _tmp2 = d_buffer0;
    {
      unsigned offsets__tmp1 = 216;
      unsigned offsets_AminusT = AminusT_offset;
      unsigned offsets__tmp2 = 216;


      device_gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 24, 9, 9, 1.0, _tmp1, 24, AminusT, 9, 0.0, _tmp2, 24, offsets__tmp1, offsets_AminusT, offsets__tmp2, num_elements);
    }
    {
      unsigned offsets_rDivM = 0;
      unsigned offsets__tmp2 = 216;
      unsigned offsets_Q = Q_offset;


      device_gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 56, 9, 21, 1.0, rDivM(1), 56, _tmp2, 24, 1.0, Q, 56, offsets_rDivM, offsets__tmp2, offsets_Q, num_elements);
    }
    tmp_manager.free();
    tmp_manager.free();
  }
  void kernel::neighboringFlux::execute14() {
    assert(AminusT != nullptr);
    assert(I != nullptr);
    assert(Q != nullptr);
    assert(fP(2) != nullptr);
    assert(rDivM(1) != nullptr);
    assert(rT(0) != nullptr);
    assert(kernel::neighboringFlux::num_elements != 0);
    // get device temporary memory menager
    DeviceTemporaryMemoryMenager &tmp_manager = DeviceTemporaryMemoryMenager::get_instance();
    double *_tmp0, *_tmp1, *_tmp2;
    double *d_buffer0 = (double*)tmp_manager.get_mem(216 * num_elements * sizeof(double));
    double *d_buffer1 = (double*)tmp_manager.get_mem(216 * num_elements * sizeof(double));
    _tmp0 = d_buffer0;
    {
      unsigned offsets_rT = 0;
      unsigned offsets_I = I_offset;
      unsigned offsets__tmp0 = 216;


      device_gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 24, 9, 56, 1.0, rT(0), 24, I, 56, 0.0, _tmp0, 24, offsets_rT, offsets_I, offsets__tmp0, num_elements);
    }
    _tmp1 = d_buffer1;
    {
      unsigned offsets_fP = 0;
      unsigned offsets__tmp0 = 216;
      unsigned offsets__tmp1 = 216;


      device_gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 24, 9, 21, 1.0, fP(2), 24, _tmp0, 24, 0.0, _tmp1, 24, offsets_fP, offsets__tmp0, offsets__tmp1, num_elements);
    }
    _tmp2 = d_buffer0;
    {
      unsigned offsets__tmp1 = 216;
      unsigned offsets_AminusT = AminusT_offset;
      unsigned offsets__tmp2 = 216;


      device_gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 24, 9, 9, 1.0, _tmp1, 24, AminusT, 9, 0.0, _tmp2, 24, offsets__tmp1, offsets_AminusT, offsets__tmp2, num_elements);
    }
    {
      unsigned offsets_rDivM = 0;
      unsigned offsets__tmp2 = 216;
      unsigned offsets_Q = Q_offset;


      device_gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 56, 9, 21, 1.0, rDivM(1), 56, _tmp2, 24, 1.0, Q, 56, offsets_rDivM, offsets__tmp2, offsets_Q, num_elements);
    }
    tmp_manager.free();
    tmp_manager.free();
  }
  void kernel::neighboringFlux::execute15() {
    assert(AminusT != nullptr);
    assert(I != nullptr);
    assert(Q != nullptr);
    assert(fP(0) != nullptr);
    assert(rDivM(1) != nullptr);
    assert(rT(1) != nullptr);
    assert(kernel::neighboringFlux::num_elements != 0);
    // get device temporary memory menager
    DeviceTemporaryMemoryMenager &tmp_manager = DeviceTemporaryMemoryMenager::get_instance();
    double *_tmp0, *_tmp1, *_tmp2;
    double *d_buffer0 = (double*)tmp_manager.get_mem(216 * num_elements * sizeof(double));
    double *d_buffer1 = (double*)tmp_manager.get_mem(216 * num_elements * sizeof(double));
    _tmp0 = d_buffer0;
    {
      unsigned offsets_rT = 0;
      unsigned offsets_I = I_offset;
      unsigned offsets__tmp0 = 216;


      device_gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 24, 9, 56, 1.0, rT(1), 24, I, 56, 0.0, _tmp0, 24, offsets_rT, offsets_I, offsets__tmp0, num_elements);
    }
    _tmp1 = d_buffer1;
    {
      unsigned offsets_fP = 0;
      unsigned offsets__tmp0 = 216;
      unsigned offsets__tmp1 = 216;


      device_gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 24, 9, 21, 1.0, fP(0), 24, _tmp0, 24, 0.0, _tmp1, 24, offsets_fP, offsets__tmp0, offsets__tmp1, num_elements);
    }
    _tmp2 = d_buffer0;
    {
      unsigned offsets__tmp1 = 216;
      unsigned offsets_AminusT = AminusT_offset;
      unsigned offsets__tmp2 = 216;


      device_gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 24, 9, 9, 1.0, _tmp1, 24, AminusT, 9, 0.0, _tmp2, 24, offsets__tmp1, offsets_AminusT, offsets__tmp2, num_elements);
    }
    {
      unsigned offsets_rDivM = 0;
      unsigned offsets__tmp2 = 216;
      unsigned offsets_Q = Q_offset;


      device_gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 56, 9, 21, 1.0, rDivM(1), 56, _tmp2, 24, 1.0, Q, 56, offsets_rDivM, offsets__tmp2, offsets_Q, num_elements);
    }
    tmp_manager.free();
    tmp_manager.free();
  }
  void kernel::neighboringFlux::execute16() {
    assert(AminusT != nullptr);
    assert(I != nullptr);
    assert(Q != nullptr);
    assert(fP(1) != nullptr);
    assert(rDivM(1) != nullptr);
    assert(rT(1) != nullptr);
    assert(kernel::neighboringFlux::num_elements != 0);
    // get device temporary memory menager
    DeviceTemporaryMemoryMenager &tmp_manager = DeviceTemporaryMemoryMenager::get_instance();
    double *_tmp0, *_tmp1, *_tmp2;
    double *d_buffer0 = (double*)tmp_manager.get_mem(216 * num_elements * sizeof(double));
    double *d_buffer1 = (double*)tmp_manager.get_mem(216 * num_elements * sizeof(double));
    _tmp0 = d_buffer0;
    {
      unsigned offsets_rT = 0;
      unsigned offsets_I = I_offset;
      unsigned offsets__tmp0 = 216;


      device_gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 24, 9, 56, 1.0, rT(1), 24, I, 56, 0.0, _tmp0, 24, offsets_rT, offsets_I, offsets__tmp0, num_elements);
    }
    _tmp1 = d_buffer1;
    {
      unsigned offsets_fP = 0;
      unsigned offsets__tmp0 = 216;
      unsigned offsets__tmp1 = 216;


      device_gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 24, 9, 21, 1.0, fP(1), 24, _tmp0, 24, 0.0, _tmp1, 24, offsets_fP, offsets__tmp0, offsets__tmp1, num_elements);
    }
    _tmp2 = d_buffer0;
    {
      unsigned offsets__tmp1 = 216;
      unsigned offsets_AminusT = AminusT_offset;
      unsigned offsets__tmp2 = 216;


      device_gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 24, 9, 9, 1.0, _tmp1, 24, AminusT, 9, 0.0, _tmp2, 24, offsets__tmp1, offsets_AminusT, offsets__tmp2, num_elements);
    }
    {
      unsigned offsets_rDivM = 0;
      unsigned offsets__tmp2 = 216;
      unsigned offsets_Q = Q_offset;


      device_gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 56, 9, 21, 1.0, rDivM(1), 56, _tmp2, 24, 1.0, Q, 56, offsets_rDivM, offsets__tmp2, offsets_Q, num_elements);
    }
    tmp_manager.free();
    tmp_manager.free();
  }
  void kernel::neighboringFlux::execute17() {
    assert(AminusT != nullptr);
    assert(I != nullptr);
    assert(Q != nullptr);
    assert(fP(2) != nullptr);
    assert(rDivM(1) != nullptr);
    assert(rT(1) != nullptr);
    assert(kernel::neighboringFlux::num_elements != 0);
    // get device temporary memory menager
    DeviceTemporaryMemoryMenager &tmp_manager = DeviceTemporaryMemoryMenager::get_instance();
    double *_tmp0, *_tmp1, *_tmp2;
    double *d_buffer0 = (double*)tmp_manager.get_mem(216 * num_elements * sizeof(double));
    double *d_buffer1 = (double*)tmp_manager.get_mem(216 * num_elements * sizeof(double));
    _tmp0 = d_buffer0;
    {
      unsigned offsets_rT = 0;
      unsigned offsets_I = I_offset;
      unsigned offsets__tmp0 = 216;


      device_gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 24, 9, 56, 1.0, rT(1), 24, I, 56, 0.0, _tmp0, 24, offsets_rT, offsets_I, offsets__tmp0, num_elements);
    }
    _tmp1 = d_buffer1;
    {
      unsigned offsets_fP = 0;
      unsigned offsets__tmp0 = 216;
      unsigned offsets__tmp1 = 216;


      device_gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 24, 9, 21, 1.0, fP(2), 24, _tmp0, 24, 0.0, _tmp1, 24, offsets_fP, offsets__tmp0, offsets__tmp1, num_elements);
    }
    _tmp2 = d_buffer0;
    {
      unsigned offsets__tmp1 = 216;
      unsigned offsets_AminusT = AminusT_offset;
      unsigned offsets__tmp2 = 216;


      device_gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 24, 9, 9, 1.0, _tmp1, 24, AminusT, 9, 0.0, _tmp2, 24, offsets__tmp1, offsets_AminusT, offsets__tmp2, num_elements);
    }
    {
      unsigned offsets_rDivM = 0;
      unsigned offsets__tmp2 = 216;
      unsigned offsets_Q = Q_offset;


      device_gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 56, 9, 21, 1.0, rDivM(1), 56, _tmp2, 24, 1.0, Q, 56, offsets_rDivM, offsets__tmp2, offsets_Q, num_elements);
    }
    tmp_manager.free();
    tmp_manager.free();
  }
  void kernel::neighboringFlux::execute18() {
    assert(AminusT != nullptr);
    assert(I != nullptr);
    assert(Q != nullptr);
    assert(fP(0) != nullptr);
    assert(rDivM(1) != nullptr);
    assert(rT(2) != nullptr);
    assert(kernel::neighboringFlux::num_elements != 0);
    // get device temporary memory menager
    DeviceTemporaryMemoryMenager &tmp_manager = DeviceTemporaryMemoryMenager::get_instance();
    double *_tmp0, *_tmp1, *_tmp2;
    double *d_buffer0 = (double*)tmp_manager.get_mem(216 * num_elements * sizeof(double));
    double *d_buffer1 = (double*)tmp_manager.get_mem(216 * num_elements * sizeof(double));
    _tmp0 = d_buffer0;
    {
      unsigned offsets_rT = 0;
      unsigned offsets_I = I_offset;
      unsigned offsets__tmp0 = 216;


      device_gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 24, 9, 56, 1.0, rT(2), 24, I, 56, 0.0, _tmp0, 24, offsets_rT, offsets_I, offsets__tmp0, num_elements);
    }
    _tmp1 = d_buffer1;
    {
      unsigned offsets_fP = 0;
      unsigned offsets__tmp0 = 216;
      unsigned offsets__tmp1 = 216;


      device_gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 24, 9, 21, 1.0, fP(0), 24, _tmp0, 24, 0.0, _tmp1, 24, offsets_fP, offsets__tmp0, offsets__tmp1, num_elements);
    }
    _tmp2 = d_buffer0;
    {
      unsigned offsets__tmp1 = 216;
      unsigned offsets_AminusT = AminusT_offset;
      unsigned offsets__tmp2 = 216;


      device_gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 24, 9, 9, 1.0, _tmp1, 24, AminusT, 9, 0.0, _tmp2, 24, offsets__tmp1, offsets_AminusT, offsets__tmp2, num_elements);
    }
    {
      unsigned offsets_rDivM = 0;
      unsigned offsets__tmp2 = 216;
      unsigned offsets_Q = Q_offset;


      device_gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 56, 9, 21, 1.0, rDivM(1), 56, _tmp2, 24, 1.0, Q, 56, offsets_rDivM, offsets__tmp2, offsets_Q, num_elements);
    }
    tmp_manager.free();
    tmp_manager.free();
  }
  void kernel::neighboringFlux::execute19() {
    assert(AminusT != nullptr);
    assert(I != nullptr);
    assert(Q != nullptr);
    assert(fP(1) != nullptr);
    assert(rDivM(1) != nullptr);
    assert(rT(2) != nullptr);
    assert(kernel::neighboringFlux::num_elements != 0);
    // get device temporary memory menager
    DeviceTemporaryMemoryMenager &tmp_manager = DeviceTemporaryMemoryMenager::get_instance();
    double *_tmp0, *_tmp1, *_tmp2;
    double *d_buffer0 = (double*)tmp_manager.get_mem(216 * num_elements * sizeof(double));
    double *d_buffer1 = (double*)tmp_manager.get_mem(216 * num_elements * sizeof(double));
    _tmp0 = d_buffer0;
    {
      unsigned offsets_rT = 0;
      unsigned offsets_I = I_offset;
      unsigned offsets__tmp0 = 216;


      device_gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 24, 9, 56, 1.0, rT(2), 24, I, 56, 0.0, _tmp0, 24, offsets_rT, offsets_I, offsets__tmp0, num_elements);
    }
    _tmp1 = d_buffer1;
    {
      unsigned offsets_fP = 0;
      unsigned offsets__tmp0 = 216;
      unsigned offsets__tmp1 = 216;


      device_gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 24, 9, 21, 1.0, fP(1), 24, _tmp0, 24, 0.0, _tmp1, 24, offsets_fP, offsets__tmp0, offsets__tmp1, num_elements);
    }
    _tmp2 = d_buffer0;
    {
      unsigned offsets__tmp1 = 216;
      unsigned offsets_AminusT = AminusT_offset;
      unsigned offsets__tmp2 = 216;


      device_gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 24, 9, 9, 1.0, _tmp1, 24, AminusT, 9, 0.0, _tmp2, 24, offsets__tmp1, offsets_AminusT, offsets__tmp2, num_elements);
    }
    {
      unsigned offsets_rDivM = 0;
      unsigned offsets__tmp2 = 216;
      unsigned offsets_Q = Q_offset;


      device_gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 56, 9, 21, 1.0, rDivM(1), 56, _tmp2, 24, 1.0, Q, 56, offsets_rDivM, offsets__tmp2, offsets_Q, num_elements);
    }
    tmp_manager.free();
    tmp_manager.free();
  }
  void kernel::neighboringFlux::execute20() {
    assert(AminusT != nullptr);
    assert(I != nullptr);
    assert(Q != nullptr);
    assert(fP(2) != nullptr);
    assert(rDivM(1) != nullptr);
    assert(rT(2) != nullptr);
    assert(kernel::neighboringFlux::num_elements != 0);
    // get device temporary memory menager
    DeviceTemporaryMemoryMenager &tmp_manager = DeviceTemporaryMemoryMenager::get_instance();
    double *_tmp0, *_tmp1, *_tmp2;
    double *d_buffer0 = (double*)tmp_manager.get_mem(216 * num_elements * sizeof(double));
    double *d_buffer1 = (double*)tmp_manager.get_mem(216 * num_elements * sizeof(double));
    _tmp0 = d_buffer0;
    {
      unsigned offsets_rT = 0;
      unsigned offsets_I = I_offset;
      unsigned offsets__tmp0 = 216;


      device_gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 24, 9, 56, 1.0, rT(2), 24, I, 56, 0.0, _tmp0, 24, offsets_rT, offsets_I, offsets__tmp0, num_elements);
    }
    _tmp1 = d_buffer1;
    {
      unsigned offsets_fP = 0;
      unsigned offsets__tmp0 = 216;
      unsigned offsets__tmp1 = 216;


      device_gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 24, 9, 21, 1.0, fP(2), 24, _tmp0, 24, 0.0, _tmp1, 24, offsets_fP, offsets__tmp0, offsets__tmp1, num_elements);
    }
    _tmp2 = d_buffer0;
    {
      unsigned offsets__tmp1 = 216;
      unsigned offsets_AminusT = AminusT_offset;
      unsigned offsets__tmp2 = 216;


      device_gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 24, 9, 9, 1.0, _tmp1, 24, AminusT, 9, 0.0, _tmp2, 24, offsets__tmp1, offsets_AminusT, offsets__tmp2, num_elements);
    }
    {
      unsigned offsets_rDivM = 0;
      unsigned offsets__tmp2 = 216;
      unsigned offsets_Q = Q_offset;


      device_gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 56, 9, 21, 1.0, rDivM(1), 56, _tmp2, 24, 1.0, Q, 56, offsets_rDivM, offsets__tmp2, offsets_Q, num_elements);
    }
    tmp_manager.free();
    tmp_manager.free();
  }
  void kernel::neighboringFlux::execute21() {
    assert(AminusT != nullptr);
    assert(I != nullptr);
    assert(Q != nullptr);
    assert(fP(0) != nullptr);
    assert(rDivM(1) != nullptr);
    assert(rT(3) != nullptr);
    assert(kernel::neighboringFlux::num_elements != 0);
    // get device temporary memory menager
    DeviceTemporaryMemoryMenager &tmp_manager = DeviceTemporaryMemoryMenager::get_instance();
    double *_tmp0, *_tmp1, *_tmp2;
    double *d_buffer0 = (double*)tmp_manager.get_mem(216 * num_elements * sizeof(double));
    double *d_buffer1 = (double*)tmp_manager.get_mem(216 * num_elements * sizeof(double));
    _tmp0 = d_buffer0;
    {
      unsigned offsets_rT = 0;
      unsigned offsets_I = I_offset;
      unsigned offsets__tmp0 = 216;


      device_gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 24, 9, 56, 1.0, rT(3), 24, I, 56, 0.0, _tmp0, 24, offsets_rT, offsets_I, offsets__tmp0, num_elements);
    }
    _tmp1 = d_buffer1;
    {
      unsigned offsets_fP = 0;
      unsigned offsets__tmp0 = 216;
      unsigned offsets__tmp1 = 216;


      device_gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 24, 9, 21, 1.0, fP(0), 24, _tmp0, 24, 0.0, _tmp1, 24, offsets_fP, offsets__tmp0, offsets__tmp1, num_elements);
    }
    _tmp2 = d_buffer0;
    {
      unsigned offsets__tmp1 = 216;
      unsigned offsets_AminusT = AminusT_offset;
      unsigned offsets__tmp2 = 216;


      device_gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 24, 9, 9, 1.0, _tmp1, 24, AminusT, 9, 0.0, _tmp2, 24, offsets__tmp1, offsets_AminusT, offsets__tmp2, num_elements);
    }
    {
      unsigned offsets_rDivM = 0;
      unsigned offsets__tmp2 = 216;
      unsigned offsets_Q = Q_offset;


      device_gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 56, 9, 21, 1.0, rDivM(1), 56, _tmp2, 24, 1.0, Q, 56, offsets_rDivM, offsets__tmp2, offsets_Q, num_elements);
    }
    tmp_manager.free();
    tmp_manager.free();
  }
  void kernel::neighboringFlux::execute22() {
    assert(AminusT != nullptr);
    assert(I != nullptr);
    assert(Q != nullptr);
    assert(fP(1) != nullptr);
    assert(rDivM(1) != nullptr);
    assert(rT(3) != nullptr);
    assert(kernel::neighboringFlux::num_elements != 0);
    // get device temporary memory menager
    DeviceTemporaryMemoryMenager &tmp_manager = DeviceTemporaryMemoryMenager::get_instance();
    double *_tmp0, *_tmp1, *_tmp2;
    double *d_buffer0 = (double*)tmp_manager.get_mem(216 * num_elements * sizeof(double));
    double *d_buffer1 = (double*)tmp_manager.get_mem(216 * num_elements * sizeof(double));
    _tmp0 = d_buffer0;
    {
      unsigned offsets_rT = 0;
      unsigned offsets_I = I_offset;
      unsigned offsets__tmp0 = 216;


      device_gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 24, 9, 56, 1.0, rT(3), 24, I, 56, 0.0, _tmp0, 24, offsets_rT, offsets_I, offsets__tmp0, num_elements);
    }
    _tmp1 = d_buffer1;
    {
      unsigned offsets_fP = 0;
      unsigned offsets__tmp0 = 216;
      unsigned offsets__tmp1 = 216;


      device_gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 24, 9, 21, 1.0, fP(1), 24, _tmp0, 24, 0.0, _tmp1, 24, offsets_fP, offsets__tmp0, offsets__tmp1, num_elements);
    }
    _tmp2 = d_buffer0;
    {
      unsigned offsets__tmp1 = 216;
      unsigned offsets_AminusT = AminusT_offset;
      unsigned offsets__tmp2 = 216;


      device_gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 24, 9, 9, 1.0, _tmp1, 24, AminusT, 9, 0.0, _tmp2, 24, offsets__tmp1, offsets_AminusT, offsets__tmp2, num_elements);
    }
    {
      unsigned offsets_rDivM = 0;
      unsigned offsets__tmp2 = 216;
      unsigned offsets_Q = Q_offset;


      device_gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 56, 9, 21, 1.0, rDivM(1), 56, _tmp2, 24, 1.0, Q, 56, offsets_rDivM, offsets__tmp2, offsets_Q, num_elements);
    }
    tmp_manager.free();
    tmp_manager.free();
  }
  void kernel::neighboringFlux::execute23() {
    assert(AminusT != nullptr);
    assert(I != nullptr);
    assert(Q != nullptr);
    assert(fP(2) != nullptr);
    assert(rDivM(1) != nullptr);
    assert(rT(3) != nullptr);
    assert(kernel::neighboringFlux::num_elements != 0);
    // get device temporary memory menager
    DeviceTemporaryMemoryMenager &tmp_manager = DeviceTemporaryMemoryMenager::get_instance();
    double *_tmp0, *_tmp1, *_tmp2;
    double *d_buffer0 = (double*)tmp_manager.get_mem(216 * num_elements * sizeof(double));
    double *d_buffer1 = (double*)tmp_manager.get_mem(216 * num_elements * sizeof(double));
    _tmp0 = d_buffer0;
    {
      unsigned offsets_rT = 0;
      unsigned offsets_I = I_offset;
      unsigned offsets__tmp0 = 216;


      device_gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 24, 9, 56, 1.0, rT(3), 24, I, 56, 0.0, _tmp0, 24, offsets_rT, offsets_I, offsets__tmp0, num_elements);
    }
    _tmp1 = d_buffer1;
    {
      unsigned offsets_fP = 0;
      unsigned offsets__tmp0 = 216;
      unsigned offsets__tmp1 = 216;


      device_gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 24, 9, 21, 1.0, fP(2), 24, _tmp0, 24, 0.0, _tmp1, 24, offsets_fP, offsets__tmp0, offsets__tmp1, num_elements);
    }
    _tmp2 = d_buffer0;
    {
      unsigned offsets__tmp1 = 216;
      unsigned offsets_AminusT = AminusT_offset;
      unsigned offsets__tmp2 = 216;


      device_gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 24, 9, 9, 1.0, _tmp1, 24, AminusT, 9, 0.0, _tmp2, 24, offsets__tmp1, offsets_AminusT, offsets__tmp2, num_elements);
    }
    {
      unsigned offsets_rDivM = 0;
      unsigned offsets__tmp2 = 216;
      unsigned offsets_Q = Q_offset;


      device_gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 56, 9, 21, 1.0, rDivM(1), 56, _tmp2, 24, 1.0, Q, 56, offsets_rDivM, offsets__tmp2, offsets_Q, num_elements);
    }
    tmp_manager.free();
    tmp_manager.free();
  }
  void kernel::neighboringFlux::execute24() {
    assert(AminusT != nullptr);
    assert(I != nullptr);
    assert(Q != nullptr);
    assert(fP(0) != nullptr);
    assert(rDivM(2) != nullptr);
    assert(rT(0) != nullptr);
    assert(kernel::neighboringFlux::num_elements != 0);
    // get device temporary memory menager
    DeviceTemporaryMemoryMenager &tmp_manager = DeviceTemporaryMemoryMenager::get_instance();
    double *_tmp0, *_tmp1, *_tmp2;
    double *d_buffer0 = (double*)tmp_manager.get_mem(216 * num_elements * sizeof(double));
    double *d_buffer1 = (double*)tmp_manager.get_mem(216 * num_elements * sizeof(double));
    _tmp0 = d_buffer0;
    {
      unsigned offsets_rT = 0;
      unsigned offsets_I = I_offset;
      unsigned offsets__tmp0 = 216;


      device_gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 24, 9, 56, 1.0, rT(0), 24, I, 56, 0.0, _tmp0, 24, offsets_rT, offsets_I, offsets__tmp0, num_elements);
    }
    _tmp1 = d_buffer1;
    {
      unsigned offsets_fP = 0;
      unsigned offsets__tmp0 = 216;
      unsigned offsets__tmp1 = 216;


      device_gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 24, 9, 21, 1.0, fP(0), 24, _tmp0, 24, 0.0, _tmp1, 24, offsets_fP, offsets__tmp0, offsets__tmp1, num_elements);
    }
    _tmp2 = d_buffer0;
    {
      unsigned offsets__tmp1 = 216;
      unsigned offsets_AminusT = AminusT_offset;
      unsigned offsets__tmp2 = 216;


      device_gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 24, 9, 9, 1.0, _tmp1, 24, AminusT, 9, 0.0, _tmp2, 24, offsets__tmp1, offsets_AminusT, offsets__tmp2, num_elements);
    }
    {
      unsigned offsets_rDivM = 0;
      unsigned offsets__tmp2 = 216;
      unsigned offsets_Q = Q_offset;


      device_gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 56, 9, 21, 1.0, rDivM(2), 56, _tmp2, 24, 1.0, Q, 56, offsets_rDivM, offsets__tmp2, offsets_Q, num_elements);
    }
    tmp_manager.free();
    tmp_manager.free();
  }
  void kernel::neighboringFlux::execute25() {
    assert(AminusT != nullptr);
    assert(I != nullptr);
    assert(Q != nullptr);
    assert(fP(1) != nullptr);
    assert(rDivM(2) != nullptr);
    assert(rT(0) != nullptr);
    assert(kernel::neighboringFlux::num_elements != 0);
    // get device temporary memory menager
    DeviceTemporaryMemoryMenager &tmp_manager = DeviceTemporaryMemoryMenager::get_instance();
    double *_tmp0, *_tmp1, *_tmp2;
    double *d_buffer0 = (double*)tmp_manager.get_mem(216 * num_elements * sizeof(double));
    double *d_buffer1 = (double*)tmp_manager.get_mem(216 * num_elements * sizeof(double));
    _tmp0 = d_buffer0;
    {
      unsigned offsets_rT = 0;
      unsigned offsets_I = I_offset;
      unsigned offsets__tmp0 = 216;


      device_gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 24, 9, 56, 1.0, rT(0), 24, I, 56, 0.0, _tmp0, 24, offsets_rT, offsets_I, offsets__tmp0, num_elements);
    }
    _tmp1 = d_buffer1;
    {
      unsigned offsets_fP = 0;
      unsigned offsets__tmp0 = 216;
      unsigned offsets__tmp1 = 216;


      device_gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 24, 9, 21, 1.0, fP(1), 24, _tmp0, 24, 0.0, _tmp1, 24, offsets_fP, offsets__tmp0, offsets__tmp1, num_elements);
    }
    _tmp2 = d_buffer0;
    {
      unsigned offsets__tmp1 = 216;
      unsigned offsets_AminusT = AminusT_offset;
      unsigned offsets__tmp2 = 216;


      device_gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 24, 9, 9, 1.0, _tmp1, 24, AminusT, 9, 0.0, _tmp2, 24, offsets__tmp1, offsets_AminusT, offsets__tmp2, num_elements);
    }
    {
      unsigned offsets_rDivM = 0;
      unsigned offsets__tmp2 = 216;
      unsigned offsets_Q = Q_offset;


      device_gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 56, 9, 21, 1.0, rDivM(2), 56, _tmp2, 24, 1.0, Q, 56, offsets_rDivM, offsets__tmp2, offsets_Q, num_elements);
    }
    tmp_manager.free();
    tmp_manager.free();
  }
  void kernel::neighboringFlux::execute26() {
    assert(AminusT != nullptr);
    assert(I != nullptr);
    assert(Q != nullptr);
    assert(fP(2) != nullptr);
    assert(rDivM(2) != nullptr);
    assert(rT(0) != nullptr);
    assert(kernel::neighboringFlux::num_elements != 0);
    // get device temporary memory menager
    DeviceTemporaryMemoryMenager &tmp_manager = DeviceTemporaryMemoryMenager::get_instance();
    double *_tmp0, *_tmp1, *_tmp2;
    double *d_buffer0 = (double*)tmp_manager.get_mem(216 * num_elements * sizeof(double));
    double *d_buffer1 = (double*)tmp_manager.get_mem(216 * num_elements * sizeof(double));
    _tmp0 = d_buffer0;
    {
      unsigned offsets_rT = 0;
      unsigned offsets_I = I_offset;
      unsigned offsets__tmp0 = 216;


      device_gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 24, 9, 56, 1.0, rT(0), 24, I, 56, 0.0, _tmp0, 24, offsets_rT, offsets_I, offsets__tmp0, num_elements);
    }
    _tmp1 = d_buffer1;
    {
      unsigned offsets_fP = 0;
      unsigned offsets__tmp0 = 216;
      unsigned offsets__tmp1 = 216;


      device_gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 24, 9, 21, 1.0, fP(2), 24, _tmp0, 24, 0.0, _tmp1, 24, offsets_fP, offsets__tmp0, offsets__tmp1, num_elements);
    }
    _tmp2 = d_buffer0;
    {
      unsigned offsets__tmp1 = 216;
      unsigned offsets_AminusT = AminusT_offset;
      unsigned offsets__tmp2 = 216;


      device_gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 24, 9, 9, 1.0, _tmp1, 24, AminusT, 9, 0.0, _tmp2, 24, offsets__tmp1, offsets_AminusT, offsets__tmp2, num_elements);
    }
    {
      unsigned offsets_rDivM = 0;
      unsigned offsets__tmp2 = 216;
      unsigned offsets_Q = Q_offset;


      device_gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 56, 9, 21, 1.0, rDivM(2), 56, _tmp2, 24, 1.0, Q, 56, offsets_rDivM, offsets__tmp2, offsets_Q, num_elements);
    }
    tmp_manager.free();
    tmp_manager.free();
  }
  void kernel::neighboringFlux::execute27() {
    assert(AminusT != nullptr);
    assert(I != nullptr);
    assert(Q != nullptr);
    assert(fP(0) != nullptr);
    assert(rDivM(2) != nullptr);
    assert(rT(1) != nullptr);
    assert(kernel::neighboringFlux::num_elements != 0);
    // get device temporary memory menager
    DeviceTemporaryMemoryMenager &tmp_manager = DeviceTemporaryMemoryMenager::get_instance();
    double *_tmp0, *_tmp1, *_tmp2;
    double *d_buffer0 = (double*)tmp_manager.get_mem(216 * num_elements * sizeof(double));
    double *d_buffer1 = (double*)tmp_manager.get_mem(216 * num_elements * sizeof(double));
    _tmp0 = d_buffer0;
    {
      unsigned offsets_rT = 0;
      unsigned offsets_I = I_offset;
      unsigned offsets__tmp0 = 216;


      device_gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 24, 9, 56, 1.0, rT(1), 24, I, 56, 0.0, _tmp0, 24, offsets_rT, offsets_I, offsets__tmp0, num_elements);
    }
    _tmp1 = d_buffer1;
    {
      unsigned offsets_fP = 0;
      unsigned offsets__tmp0 = 216;
      unsigned offsets__tmp1 = 216;


      device_gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 24, 9, 21, 1.0, fP(0), 24, _tmp0, 24, 0.0, _tmp1, 24, offsets_fP, offsets__tmp0, offsets__tmp1, num_elements);
    }
    _tmp2 = d_buffer0;
    {
      unsigned offsets__tmp1 = 216;
      unsigned offsets_AminusT = AminusT_offset;
      unsigned offsets__tmp2 = 216;


      device_gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 24, 9, 9, 1.0, _tmp1, 24, AminusT, 9, 0.0, _tmp2, 24, offsets__tmp1, offsets_AminusT, offsets__tmp2, num_elements);
    }
    {
      unsigned offsets_rDivM = 0;
      unsigned offsets__tmp2 = 216;
      unsigned offsets_Q = Q_offset;


      device_gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 56, 9, 21, 1.0, rDivM(2), 56, _tmp2, 24, 1.0, Q, 56, offsets_rDivM, offsets__tmp2, offsets_Q, num_elements);
    }
    tmp_manager.free();
    tmp_manager.free();
  }
  void kernel::neighboringFlux::execute28() {
    assert(AminusT != nullptr);
    assert(I != nullptr);
    assert(Q != nullptr);
    assert(fP(1) != nullptr);
    assert(rDivM(2) != nullptr);
    assert(rT(1) != nullptr);
    assert(kernel::neighboringFlux::num_elements != 0);
    // get device temporary memory menager
    DeviceTemporaryMemoryMenager &tmp_manager = DeviceTemporaryMemoryMenager::get_instance();
    double *_tmp0, *_tmp1, *_tmp2;
    double *d_buffer0 = (double*)tmp_manager.get_mem(216 * num_elements * sizeof(double));
    double *d_buffer1 = (double*)tmp_manager.get_mem(216 * num_elements * sizeof(double));
    _tmp0 = d_buffer0;
    {
      unsigned offsets_rT = 0;
      unsigned offsets_I = I_offset;
      unsigned offsets__tmp0 = 216;


      device_gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 24, 9, 56, 1.0, rT(1), 24, I, 56, 0.0, _tmp0, 24, offsets_rT, offsets_I, offsets__tmp0, num_elements);
    }
    _tmp1 = d_buffer1;
    {
      unsigned offsets_fP = 0;
      unsigned offsets__tmp0 = 216;
      unsigned offsets__tmp1 = 216;


      device_gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 24, 9, 21, 1.0, fP(1), 24, _tmp0, 24, 0.0, _tmp1, 24, offsets_fP, offsets__tmp0, offsets__tmp1, num_elements);
    }
    _tmp2 = d_buffer0;
    {
      unsigned offsets__tmp1 = 216;
      unsigned offsets_AminusT = AminusT_offset;
      unsigned offsets__tmp2 = 216;


      device_gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 24, 9, 9, 1.0, _tmp1, 24, AminusT, 9, 0.0, _tmp2, 24, offsets__tmp1, offsets_AminusT, offsets__tmp2, num_elements);
    }
    {
      unsigned offsets_rDivM = 0;
      unsigned offsets__tmp2 = 216;
      unsigned offsets_Q = Q_offset;


      device_gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 56, 9, 21, 1.0, rDivM(2), 56, _tmp2, 24, 1.0, Q, 56, offsets_rDivM, offsets__tmp2, offsets_Q, num_elements);
    }
    tmp_manager.free();
    tmp_manager.free();
  }
  void kernel::neighboringFlux::execute29() {
    assert(AminusT != nullptr);
    assert(I != nullptr);
    assert(Q != nullptr);
    assert(fP(2) != nullptr);
    assert(rDivM(2) != nullptr);
    assert(rT(1) != nullptr);
    assert(kernel::neighboringFlux::num_elements != 0);
    // get device temporary memory menager
    DeviceTemporaryMemoryMenager &tmp_manager = DeviceTemporaryMemoryMenager::get_instance();
    double *_tmp0, *_tmp1, *_tmp2;
    double *d_buffer0 = (double*)tmp_manager.get_mem(216 * num_elements * sizeof(double));
    double *d_buffer1 = (double*)tmp_manager.get_mem(216 * num_elements * sizeof(double));
    _tmp0 = d_buffer0;
    {
      unsigned offsets_rT = 0;
      unsigned offsets_I = I_offset;
      unsigned offsets__tmp0 = 216;


      device_gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 24, 9, 56, 1.0, rT(1), 24, I, 56, 0.0, _tmp0, 24, offsets_rT, offsets_I, offsets__tmp0, num_elements);
    }
    _tmp1 = d_buffer1;
    {
      unsigned offsets_fP = 0;
      unsigned offsets__tmp0 = 216;
      unsigned offsets__tmp1 = 216;


      device_gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 24, 9, 21, 1.0, fP(2), 24, _tmp0, 24, 0.0, _tmp1, 24, offsets_fP, offsets__tmp0, offsets__tmp1, num_elements);
    }
    _tmp2 = d_buffer0;
    {
      unsigned offsets__tmp1 = 216;
      unsigned offsets_AminusT = AminusT_offset;
      unsigned offsets__tmp2 = 216;


      device_gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 24, 9, 9, 1.0, _tmp1, 24, AminusT, 9, 0.0, _tmp2, 24, offsets__tmp1, offsets_AminusT, offsets__tmp2, num_elements);
    }
    {
      unsigned offsets_rDivM = 0;
      unsigned offsets__tmp2 = 216;
      unsigned offsets_Q = Q_offset;


      device_gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 56, 9, 21, 1.0, rDivM(2), 56, _tmp2, 24, 1.0, Q, 56, offsets_rDivM, offsets__tmp2, offsets_Q, num_elements);
    }
    tmp_manager.free();
    tmp_manager.free();
  }
  void kernel::neighboringFlux::execute30() {
    assert(AminusT != nullptr);
    assert(I != nullptr);
    assert(Q != nullptr);
    assert(fP(0) != nullptr);
    assert(rDivM(2) != nullptr);
    assert(rT(2) != nullptr);
    assert(kernel::neighboringFlux::num_elements != 0);
    // get device temporary memory menager
    DeviceTemporaryMemoryMenager &tmp_manager = DeviceTemporaryMemoryMenager::get_instance();
    double *_tmp0, *_tmp1, *_tmp2;
    double *d_buffer0 = (double*)tmp_manager.get_mem(216 * num_elements * sizeof(double));
    double *d_buffer1 = (double*)tmp_manager.get_mem(216 * num_elements * sizeof(double));
    _tmp0 = d_buffer0;
    {
      unsigned offsets_rT = 0;
      unsigned offsets_I = I_offset;
      unsigned offsets__tmp0 = 216;


      device_gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 24, 9, 56, 1.0, rT(2), 24, I, 56, 0.0, _tmp0, 24, offsets_rT, offsets_I, offsets__tmp0, num_elements);
    }
    _tmp1 = d_buffer1;
    {
      unsigned offsets_fP = 0;
      unsigned offsets__tmp0 = 216;
      unsigned offsets__tmp1 = 216;


      device_gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 24, 9, 21, 1.0, fP(0), 24, _tmp0, 24, 0.0, _tmp1, 24, offsets_fP, offsets__tmp0, offsets__tmp1, num_elements);
    }
    _tmp2 = d_buffer0;
    {
      unsigned offsets__tmp1 = 216;
      unsigned offsets_AminusT = AminusT_offset;
      unsigned offsets__tmp2 = 216;


      device_gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 24, 9, 9, 1.0, _tmp1, 24, AminusT, 9, 0.0, _tmp2, 24, offsets__tmp1, offsets_AminusT, offsets__tmp2, num_elements);
    }
    {
      unsigned offsets_rDivM = 0;
      unsigned offsets__tmp2 = 216;
      unsigned offsets_Q = Q_offset;


      device_gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 56, 9, 21, 1.0, rDivM(2), 56, _tmp2, 24, 1.0, Q, 56, offsets_rDivM, offsets__tmp2, offsets_Q, num_elements);
    }
    tmp_manager.free();
    tmp_manager.free();
  }
  void kernel::neighboringFlux::execute31() {
    assert(AminusT != nullptr);
    assert(I != nullptr);
    assert(Q != nullptr);
    assert(fP(1) != nullptr);
    assert(rDivM(2) != nullptr);
    assert(rT(2) != nullptr);
    assert(kernel::neighboringFlux::num_elements != 0);
    // get device temporary memory menager
    DeviceTemporaryMemoryMenager &tmp_manager = DeviceTemporaryMemoryMenager::get_instance();
    double *_tmp0, *_tmp1, *_tmp2;
    double *d_buffer0 = (double*)tmp_manager.get_mem(216 * num_elements * sizeof(double));
    double *d_buffer1 = (double*)tmp_manager.get_mem(216 * num_elements * sizeof(double));
    _tmp0 = d_buffer0;
    {
      unsigned offsets_rT = 0;
      unsigned offsets_I = I_offset;
      unsigned offsets__tmp0 = 216;


      device_gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 24, 9, 56, 1.0, rT(2), 24, I, 56, 0.0, _tmp0, 24, offsets_rT, offsets_I, offsets__tmp0, num_elements);
    }
    _tmp1 = d_buffer1;
    {
      unsigned offsets_fP = 0;
      unsigned offsets__tmp0 = 216;
      unsigned offsets__tmp1 = 216;


      device_gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 24, 9, 21, 1.0, fP(1), 24, _tmp0, 24, 0.0, _tmp1, 24, offsets_fP, offsets__tmp0, offsets__tmp1, num_elements);
    }
    _tmp2 = d_buffer0;
    {
      unsigned offsets__tmp1 = 216;
      unsigned offsets_AminusT = AminusT_offset;
      unsigned offsets__tmp2 = 216;


      device_gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 24, 9, 9, 1.0, _tmp1, 24, AminusT, 9, 0.0, _tmp2, 24, offsets__tmp1, offsets_AminusT, offsets__tmp2, num_elements);
    }
    {
      unsigned offsets_rDivM = 0;
      unsigned offsets__tmp2 = 216;
      unsigned offsets_Q = Q_offset;


      device_gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 56, 9, 21, 1.0, rDivM(2), 56, _tmp2, 24, 1.0, Q, 56, offsets_rDivM, offsets__tmp2, offsets_Q, num_elements);
    }
    tmp_manager.free();
    tmp_manager.free();
  }
  void kernel::neighboringFlux::execute32() {
    assert(AminusT != nullptr);
    assert(I != nullptr);
    assert(Q != nullptr);
    assert(fP(2) != nullptr);
    assert(rDivM(2) != nullptr);
    assert(rT(2) != nullptr);
    assert(kernel::neighboringFlux::num_elements != 0);
    // get device temporary memory menager
    DeviceTemporaryMemoryMenager &tmp_manager = DeviceTemporaryMemoryMenager::get_instance();
    double *_tmp0, *_tmp1, *_tmp2;
    double *d_buffer0 = (double*)tmp_manager.get_mem(216 * num_elements * sizeof(double));
    double *d_buffer1 = (double*)tmp_manager.get_mem(216 * num_elements * sizeof(double));
    _tmp0 = d_buffer0;
    {
      unsigned offsets_rT = 0;
      unsigned offsets_I = I_offset;
      unsigned offsets__tmp0 = 216;


      device_gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 24, 9, 56, 1.0, rT(2), 24, I, 56, 0.0, _tmp0, 24, offsets_rT, offsets_I, offsets__tmp0, num_elements);
    }
    _tmp1 = d_buffer1;
    {
      unsigned offsets_fP = 0;
      unsigned offsets__tmp0 = 216;
      unsigned offsets__tmp1 = 216;


      device_gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 24, 9, 21, 1.0, fP(2), 24, _tmp0, 24, 0.0, _tmp1, 24, offsets_fP, offsets__tmp0, offsets__tmp1, num_elements);
    }
    _tmp2 = d_buffer0;
    {
      unsigned offsets__tmp1 = 216;
      unsigned offsets_AminusT = AminusT_offset;
      unsigned offsets__tmp2 = 216;


      device_gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 24, 9, 9, 1.0, _tmp1, 24, AminusT, 9, 0.0, _tmp2, 24, offsets__tmp1, offsets_AminusT, offsets__tmp2, num_elements);
    }
    {
      unsigned offsets_rDivM = 0;
      unsigned offsets__tmp2 = 216;
      unsigned offsets_Q = Q_offset;


      device_gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 56, 9, 21, 1.0, rDivM(2), 56, _tmp2, 24, 1.0, Q, 56, offsets_rDivM, offsets__tmp2, offsets_Q, num_elements);
    }
    tmp_manager.free();
    tmp_manager.free();
  }
  void kernel::neighboringFlux::execute33() {
    assert(AminusT != nullptr);
    assert(I != nullptr);
    assert(Q != nullptr);
    assert(fP(0) != nullptr);
    assert(rDivM(2) != nullptr);
    assert(rT(3) != nullptr);
    assert(kernel::neighboringFlux::num_elements != 0);
    // get device temporary memory menager
    DeviceTemporaryMemoryMenager &tmp_manager = DeviceTemporaryMemoryMenager::get_instance();
    double *_tmp0, *_tmp1, *_tmp2;
    double *d_buffer0 = (double*)tmp_manager.get_mem(216 * num_elements * sizeof(double));
    double *d_buffer1 = (double*)tmp_manager.get_mem(216 * num_elements * sizeof(double));
    _tmp0 = d_buffer0;
    {
      unsigned offsets_rT = 0;
      unsigned offsets_I = I_offset;
      unsigned offsets__tmp0 = 216;


      device_gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 24, 9, 56, 1.0, rT(3), 24, I, 56, 0.0, _tmp0, 24, offsets_rT, offsets_I, offsets__tmp0, num_elements);
    }
    _tmp1 = d_buffer1;
    {
      unsigned offsets_fP = 0;
      unsigned offsets__tmp0 = 216;
      unsigned offsets__tmp1 = 216;


      device_gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 24, 9, 21, 1.0, fP(0), 24, _tmp0, 24, 0.0, _tmp1, 24, offsets_fP, offsets__tmp0, offsets__tmp1, num_elements);
    }
    _tmp2 = d_buffer0;
    {
      unsigned offsets__tmp1 = 216;
      unsigned offsets_AminusT = AminusT_offset;
      unsigned offsets__tmp2 = 216;


      device_gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 24, 9, 9, 1.0, _tmp1, 24, AminusT, 9, 0.0, _tmp2, 24, offsets__tmp1, offsets_AminusT, offsets__tmp2, num_elements);
    }
    {
      unsigned offsets_rDivM = 0;
      unsigned offsets__tmp2 = 216;
      unsigned offsets_Q = Q_offset;


      device_gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 56, 9, 21, 1.0, rDivM(2), 56, _tmp2, 24, 1.0, Q, 56, offsets_rDivM, offsets__tmp2, offsets_Q, num_elements);
    }
    tmp_manager.free();
    tmp_manager.free();
  }
  void kernel::neighboringFlux::execute34() {
    assert(AminusT != nullptr);
    assert(I != nullptr);
    assert(Q != nullptr);
    assert(fP(1) != nullptr);
    assert(rDivM(2) != nullptr);
    assert(rT(3) != nullptr);
    assert(kernel::neighboringFlux::num_elements != 0);
    // get device temporary memory menager
    DeviceTemporaryMemoryMenager &tmp_manager = DeviceTemporaryMemoryMenager::get_instance();
    double *_tmp0, *_tmp1, *_tmp2;
    double *d_buffer0 = (double*)tmp_manager.get_mem(216 * num_elements * sizeof(double));
    double *d_buffer1 = (double*)tmp_manager.get_mem(216 * num_elements * sizeof(double));
    _tmp0 = d_buffer0;
    {
      unsigned offsets_rT = 0;
      unsigned offsets_I = I_offset;
      unsigned offsets__tmp0 = 216;


      device_gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 24, 9, 56, 1.0, rT(3), 24, I, 56, 0.0, _tmp0, 24, offsets_rT, offsets_I, offsets__tmp0, num_elements);
    }
    _tmp1 = d_buffer1;
    {
      unsigned offsets_fP = 0;
      unsigned offsets__tmp0 = 216;
      unsigned offsets__tmp1 = 216;


      device_gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 24, 9, 21, 1.0, fP(1), 24, _tmp0, 24, 0.0, _tmp1, 24, offsets_fP, offsets__tmp0, offsets__tmp1, num_elements);
    }
    _tmp2 = d_buffer0;
    {
      unsigned offsets__tmp1 = 216;
      unsigned offsets_AminusT = AminusT_offset;
      unsigned offsets__tmp2 = 216;


      device_gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 24, 9, 9, 1.0, _tmp1, 24, AminusT, 9, 0.0, _tmp2, 24, offsets__tmp1, offsets_AminusT, offsets__tmp2, num_elements);
    }
    {
      unsigned offsets_rDivM = 0;
      unsigned offsets__tmp2 = 216;
      unsigned offsets_Q = Q_offset;


      device_gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 56, 9, 21, 1.0, rDivM(2), 56, _tmp2, 24, 1.0, Q, 56, offsets_rDivM, offsets__tmp2, offsets_Q, num_elements);
    }
    tmp_manager.free();
    tmp_manager.free();
  }
  void kernel::neighboringFlux::execute35() {
    assert(AminusT != nullptr);
    assert(I != nullptr);
    assert(Q != nullptr);
    assert(fP(2) != nullptr);
    assert(rDivM(2) != nullptr);
    assert(rT(3) != nullptr);
    assert(kernel::neighboringFlux::num_elements != 0);
    // get device temporary memory menager
    DeviceTemporaryMemoryMenager &tmp_manager = DeviceTemporaryMemoryMenager::get_instance();
    double *_tmp0, *_tmp1, *_tmp2;
    double *d_buffer0 = (double*)tmp_manager.get_mem(216 * num_elements * sizeof(double));
    double *d_buffer1 = (double*)tmp_manager.get_mem(216 * num_elements * sizeof(double));
    _tmp0 = d_buffer0;
    {
      unsigned offsets_rT = 0;
      unsigned offsets_I = I_offset;
      unsigned offsets__tmp0 = 216;


      device_gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 24, 9, 56, 1.0, rT(3), 24, I, 56, 0.0, _tmp0, 24, offsets_rT, offsets_I, offsets__tmp0, num_elements);
    }
    _tmp1 = d_buffer1;
    {
      unsigned offsets_fP = 0;
      unsigned offsets__tmp0 = 216;
      unsigned offsets__tmp1 = 216;


      device_gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 24, 9, 21, 1.0, fP(2), 24, _tmp0, 24, 0.0, _tmp1, 24, offsets_fP, offsets__tmp0, offsets__tmp1, num_elements);
    }
    _tmp2 = d_buffer0;
    {
      unsigned offsets__tmp1 = 216;
      unsigned offsets_AminusT = AminusT_offset;
      unsigned offsets__tmp2 = 216;


      device_gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 24, 9, 9, 1.0, _tmp1, 24, AminusT, 9, 0.0, _tmp2, 24, offsets__tmp1, offsets_AminusT, offsets__tmp2, num_elements);
    }
    {
      unsigned offsets_rDivM = 0;
      unsigned offsets__tmp2 = 216;
      unsigned offsets_Q = Q_offset;


      device_gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 56, 9, 21, 1.0, rDivM(2), 56, _tmp2, 24, 1.0, Q, 56, offsets_rDivM, offsets__tmp2, offsets_Q, num_elements);
    }
    tmp_manager.free();
    tmp_manager.free();
  }
  void kernel::neighboringFlux::execute36() {
    assert(AminusT != nullptr);
    assert(I != nullptr);
    assert(Q != nullptr);
    assert(fP(0) != nullptr);
    assert(rDivM(3) != nullptr);
    assert(rT(0) != nullptr);
    assert(kernel::neighboringFlux::num_elements != 0);
    // get device temporary memory menager
    DeviceTemporaryMemoryMenager &tmp_manager = DeviceTemporaryMemoryMenager::get_instance();
    double *_tmp0, *_tmp1, *_tmp2;
    double *d_buffer0 = (double*)tmp_manager.get_mem(216 * num_elements * sizeof(double));
    double *d_buffer1 = (double*)tmp_manager.get_mem(216 * num_elements * sizeof(double));
    _tmp0 = d_buffer0;
    {
      unsigned offsets_rT = 0;
      unsigned offsets_I = I_offset;
      unsigned offsets__tmp0 = 216;


      device_gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 24, 9, 56, 1.0, rT(0), 24, I, 56, 0.0, _tmp0, 24, offsets_rT, offsets_I, offsets__tmp0, num_elements);
    }
    _tmp1 = d_buffer1;
    {
      unsigned offsets_fP = 0;
      unsigned offsets__tmp0 = 216;
      unsigned offsets__tmp1 = 216;


      device_gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 24, 9, 21, 1.0, fP(0), 24, _tmp0, 24, 0.0, _tmp1, 24, offsets_fP, offsets__tmp0, offsets__tmp1, num_elements);
    }
    _tmp2 = d_buffer0;
    {
      unsigned offsets__tmp1 = 216;
      unsigned offsets_AminusT = AminusT_offset;
      unsigned offsets__tmp2 = 216;


      device_gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 24, 9, 9, 1.0, _tmp1, 24, AminusT, 9, 0.0, _tmp2, 24, offsets__tmp1, offsets_AminusT, offsets__tmp2, num_elements);
    }
    {
      unsigned offsets_rDivM = 0;
      unsigned offsets__tmp2 = 216;
      unsigned offsets_Q = Q_offset;


      device_gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 56, 9, 21, 1.0, rDivM(3), 56, _tmp2, 24, 1.0, Q, 56, offsets_rDivM, offsets__tmp2, offsets_Q, num_elements);
    }
    tmp_manager.free();
    tmp_manager.free();
  }
  void kernel::neighboringFlux::execute37() {
    assert(AminusT != nullptr);
    assert(I != nullptr);
    assert(Q != nullptr);
    assert(fP(1) != nullptr);
    assert(rDivM(3) != nullptr);
    assert(rT(0) != nullptr);
    assert(kernel::neighboringFlux::num_elements != 0);
    // get device temporary memory menager
    DeviceTemporaryMemoryMenager &tmp_manager = DeviceTemporaryMemoryMenager::get_instance();
    double *_tmp0, *_tmp1, *_tmp2;
    double *d_buffer0 = (double*)tmp_manager.get_mem(216 * num_elements * sizeof(double));
    double *d_buffer1 = (double*)tmp_manager.get_mem(216 * num_elements * sizeof(double));
    _tmp0 = d_buffer0;
    {
      unsigned offsets_rT = 0;
      unsigned offsets_I = I_offset;
      unsigned offsets__tmp0 = 216;


      device_gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 24, 9, 56, 1.0, rT(0), 24, I, 56, 0.0, _tmp0, 24, offsets_rT, offsets_I, offsets__tmp0, num_elements);
    }
    _tmp1 = d_buffer1;
    {
      unsigned offsets_fP = 0;
      unsigned offsets__tmp0 = 216;
      unsigned offsets__tmp1 = 216;


      device_gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 24, 9, 21, 1.0, fP(1), 24, _tmp0, 24, 0.0, _tmp1, 24, offsets_fP, offsets__tmp0, offsets__tmp1, num_elements);
    }
    _tmp2 = d_buffer0;
    {
      unsigned offsets__tmp1 = 216;
      unsigned offsets_AminusT = AminusT_offset;
      unsigned offsets__tmp2 = 216;


      device_gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 24, 9, 9, 1.0, _tmp1, 24, AminusT, 9, 0.0, _tmp2, 24, offsets__tmp1, offsets_AminusT, offsets__tmp2, num_elements);
    }
    {
      unsigned offsets_rDivM = 0;
      unsigned offsets__tmp2 = 216;
      unsigned offsets_Q = Q_offset;


      device_gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 56, 9, 21, 1.0, rDivM(3), 56, _tmp2, 24, 1.0, Q, 56, offsets_rDivM, offsets__tmp2, offsets_Q, num_elements);
    }
    tmp_manager.free();
    tmp_manager.free();
  }
  void kernel::neighboringFlux::execute38() {
    assert(AminusT != nullptr);
    assert(I != nullptr);
    assert(Q != nullptr);
    assert(fP(2) != nullptr);
    assert(rDivM(3) != nullptr);
    assert(rT(0) != nullptr);
    assert(kernel::neighboringFlux::num_elements != 0);
    // get device temporary memory menager
    DeviceTemporaryMemoryMenager &tmp_manager = DeviceTemporaryMemoryMenager::get_instance();
    double *_tmp0, *_tmp1, *_tmp2;
    double *d_buffer0 = (double*)tmp_manager.get_mem(216 * num_elements * sizeof(double));
    double *d_buffer1 = (double*)tmp_manager.get_mem(216 * num_elements * sizeof(double));
    _tmp0 = d_buffer0;
    {
      unsigned offsets_rT = 0;
      unsigned offsets_I = I_offset;
      unsigned offsets__tmp0 = 216;


      device_gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 24, 9, 56, 1.0, rT(0), 24, I, 56, 0.0, _tmp0, 24, offsets_rT, offsets_I, offsets__tmp0, num_elements);
    }
    _tmp1 = d_buffer1;
    {
      unsigned offsets_fP = 0;
      unsigned offsets__tmp0 = 216;
      unsigned offsets__tmp1 = 216;


      device_gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 24, 9, 21, 1.0, fP(2), 24, _tmp0, 24, 0.0, _tmp1, 24, offsets_fP, offsets__tmp0, offsets__tmp1, num_elements);
    }
    _tmp2 = d_buffer0;
    {
      unsigned offsets__tmp1 = 216;
      unsigned offsets_AminusT = AminusT_offset;
      unsigned offsets__tmp2 = 216;


      device_gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 24, 9, 9, 1.0, _tmp1, 24, AminusT, 9, 0.0, _tmp2, 24, offsets__tmp1, offsets_AminusT, offsets__tmp2, num_elements);
    }
    {
      unsigned offsets_rDivM = 0;
      unsigned offsets__tmp2 = 216;
      unsigned offsets_Q = Q_offset;


      device_gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 56, 9, 21, 1.0, rDivM(3), 56, _tmp2, 24, 1.0, Q, 56, offsets_rDivM, offsets__tmp2, offsets_Q, num_elements);
    }
    tmp_manager.free();
    tmp_manager.free();
  }
  void kernel::neighboringFlux::execute39() {
    assert(AminusT != nullptr);
    assert(I != nullptr);
    assert(Q != nullptr);
    assert(fP(0) != nullptr);
    assert(rDivM(3) != nullptr);
    assert(rT(1) != nullptr);
    assert(kernel::neighboringFlux::num_elements != 0);
    // get device temporary memory menager
    DeviceTemporaryMemoryMenager &tmp_manager = DeviceTemporaryMemoryMenager::get_instance();
    double *_tmp0, *_tmp1, *_tmp2;
    double *d_buffer0 = (double*)tmp_manager.get_mem(216 * num_elements * sizeof(double));
    double *d_buffer1 = (double*)tmp_manager.get_mem(216 * num_elements * sizeof(double));
    _tmp0 = d_buffer0;
    {
      unsigned offsets_rT = 0;
      unsigned offsets_I = I_offset;
      unsigned offsets__tmp0 = 216;


      device_gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 24, 9, 56, 1.0, rT(1), 24, I, 56, 0.0, _tmp0, 24, offsets_rT, offsets_I, offsets__tmp0, num_elements);
    }
    _tmp1 = d_buffer1;
    {
      unsigned offsets_fP = 0;
      unsigned offsets__tmp0 = 216;
      unsigned offsets__tmp1 = 216;


      device_gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 24, 9, 21, 1.0, fP(0), 24, _tmp0, 24, 0.0, _tmp1, 24, offsets_fP, offsets__tmp0, offsets__tmp1, num_elements);
    }
    _tmp2 = d_buffer0;
    {
      unsigned offsets__tmp1 = 216;
      unsigned offsets_AminusT = AminusT_offset;
      unsigned offsets__tmp2 = 216;


      device_gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 24, 9, 9, 1.0, _tmp1, 24, AminusT, 9, 0.0, _tmp2, 24, offsets__tmp1, offsets_AminusT, offsets__tmp2, num_elements);
    }
    {
      unsigned offsets_rDivM = 0;
      unsigned offsets__tmp2 = 216;
      unsigned offsets_Q = Q_offset;


      device_gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 56, 9, 21, 1.0, rDivM(3), 56, _tmp2, 24, 1.0, Q, 56, offsets_rDivM, offsets__tmp2, offsets_Q, num_elements);
    }
    tmp_manager.free();
    tmp_manager.free();
  }
  void kernel::neighboringFlux::execute40() {
    assert(AminusT != nullptr);
    assert(I != nullptr);
    assert(Q != nullptr);
    assert(fP(1) != nullptr);
    assert(rDivM(3) != nullptr);
    assert(rT(1) != nullptr);
    assert(kernel::neighboringFlux::num_elements != 0);
    // get device temporary memory menager
    DeviceTemporaryMemoryMenager &tmp_manager = DeviceTemporaryMemoryMenager::get_instance();
    double *_tmp0, *_tmp1, *_tmp2;
    double *d_buffer0 = (double*)tmp_manager.get_mem(216 * num_elements * sizeof(double));
    double *d_buffer1 = (double*)tmp_manager.get_mem(216 * num_elements * sizeof(double));
    _tmp0 = d_buffer0;
    {
      unsigned offsets_rT = 0;
      unsigned offsets_I = I_offset;
      unsigned offsets__tmp0 = 216;


      device_gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 24, 9, 56, 1.0, rT(1), 24, I, 56, 0.0, _tmp0, 24, offsets_rT, offsets_I, offsets__tmp0, num_elements);
    }
    _tmp1 = d_buffer1;
    {
      unsigned offsets_fP = 0;
      unsigned offsets__tmp0 = 216;
      unsigned offsets__tmp1 = 216;


      device_gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 24, 9, 21, 1.0, fP(1), 24, _tmp0, 24, 0.0, _tmp1, 24, offsets_fP, offsets__tmp0, offsets__tmp1, num_elements);
    }
    _tmp2 = d_buffer0;
    {
      unsigned offsets__tmp1 = 216;
      unsigned offsets_AminusT = AminusT_offset;
      unsigned offsets__tmp2 = 216;


      device_gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 24, 9, 9, 1.0, _tmp1, 24, AminusT, 9, 0.0, _tmp2, 24, offsets__tmp1, offsets_AminusT, offsets__tmp2, num_elements);
    }
    {
      unsigned offsets_rDivM = 0;
      unsigned offsets__tmp2 = 216;
      unsigned offsets_Q = Q_offset;


      device_gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 56, 9, 21, 1.0, rDivM(3), 56, _tmp2, 24, 1.0, Q, 56, offsets_rDivM, offsets__tmp2, offsets_Q, num_elements);
    }
    tmp_manager.free();
    tmp_manager.free();
  }
  void kernel::neighboringFlux::execute41() {
    assert(AminusT != nullptr);
    assert(I != nullptr);
    assert(Q != nullptr);
    assert(fP(2) != nullptr);
    assert(rDivM(3) != nullptr);
    assert(rT(1) != nullptr);
    assert(kernel::neighboringFlux::num_elements != 0);
    // get device temporary memory menager
    DeviceTemporaryMemoryMenager &tmp_manager = DeviceTemporaryMemoryMenager::get_instance();
    double *_tmp0, *_tmp1, *_tmp2;
    double *d_buffer0 = (double*)tmp_manager.get_mem(216 * num_elements * sizeof(double));
    double *d_buffer1 = (double*)tmp_manager.get_mem(216 * num_elements * sizeof(double));
    _tmp0 = d_buffer0;
    {
      unsigned offsets_rT = 0;
      unsigned offsets_I = I_offset;
      unsigned offsets__tmp0 = 216;


      device_gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 24, 9, 56, 1.0, rT(1), 24, I, 56, 0.0, _tmp0, 24, offsets_rT, offsets_I, offsets__tmp0, num_elements);
    }
    _tmp1 = d_buffer1;
    {
      unsigned offsets_fP = 0;
      unsigned offsets__tmp0 = 216;
      unsigned offsets__tmp1 = 216;


      device_gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 24, 9, 21, 1.0, fP(2), 24, _tmp0, 24, 0.0, _tmp1, 24, offsets_fP, offsets__tmp0, offsets__tmp1, num_elements);
    }
    _tmp2 = d_buffer0;
    {
      unsigned offsets__tmp1 = 216;
      unsigned offsets_AminusT = AminusT_offset;
      unsigned offsets__tmp2 = 216;


      device_gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 24, 9, 9, 1.0, _tmp1, 24, AminusT, 9, 0.0, _tmp2, 24, offsets__tmp1, offsets_AminusT, offsets__tmp2, num_elements);
    }
    {
      unsigned offsets_rDivM = 0;
      unsigned offsets__tmp2 = 216;
      unsigned offsets_Q = Q_offset;


      device_gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 56, 9, 21, 1.0, rDivM(3), 56, _tmp2, 24, 1.0, Q, 56, offsets_rDivM, offsets__tmp2, offsets_Q, num_elements);
    }
    tmp_manager.free();
    tmp_manager.free();
  }
  void kernel::neighboringFlux::execute42() {
    assert(AminusT != nullptr);
    assert(I != nullptr);
    assert(Q != nullptr);
    assert(fP(0) != nullptr);
    assert(rDivM(3) != nullptr);
    assert(rT(2) != nullptr);
    assert(kernel::neighboringFlux::num_elements != 0);
    // get device temporary memory menager
    DeviceTemporaryMemoryMenager &tmp_manager = DeviceTemporaryMemoryMenager::get_instance();
    double *_tmp0, *_tmp1, *_tmp2;
    double *d_buffer0 = (double*)tmp_manager.get_mem(216 * num_elements * sizeof(double));
    double *d_buffer1 = (double*)tmp_manager.get_mem(216 * num_elements * sizeof(double));
    _tmp0 = d_buffer0;
    {
      unsigned offsets_rT = 0;
      unsigned offsets_I = I_offset;
      unsigned offsets__tmp0 = 216;


      device_gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 24, 9, 56, 1.0, rT(2), 24, I, 56, 0.0, _tmp0, 24, offsets_rT, offsets_I, offsets__tmp0, num_elements);
    }
    _tmp1 = d_buffer1;
    {
      unsigned offsets_fP = 0;
      unsigned offsets__tmp0 = 216;
      unsigned offsets__tmp1 = 216;


      device_gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 24, 9, 21, 1.0, fP(0), 24, _tmp0, 24, 0.0, _tmp1, 24, offsets_fP, offsets__tmp0, offsets__tmp1, num_elements);
    }
    _tmp2 = d_buffer0;
    {
      unsigned offsets__tmp1 = 216;
      unsigned offsets_AminusT = AminusT_offset;
      unsigned offsets__tmp2 = 216;


      device_gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 24, 9, 9, 1.0, _tmp1, 24, AminusT, 9, 0.0, _tmp2, 24, offsets__tmp1, offsets_AminusT, offsets__tmp2, num_elements);
    }
    {
      unsigned offsets_rDivM = 0;
      unsigned offsets__tmp2 = 216;
      unsigned offsets_Q = Q_offset;


      device_gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 56, 9, 21, 1.0, rDivM(3), 56, _tmp2, 24, 1.0, Q, 56, offsets_rDivM, offsets__tmp2, offsets_Q, num_elements);
    }
    tmp_manager.free();
    tmp_manager.free();
  }
  void kernel::neighboringFlux::execute43() {
    assert(AminusT != nullptr);
    assert(I != nullptr);
    assert(Q != nullptr);
    assert(fP(1) != nullptr);
    assert(rDivM(3) != nullptr);
    assert(rT(2) != nullptr);
    assert(kernel::neighboringFlux::num_elements != 0);
    // get device temporary memory menager
    DeviceTemporaryMemoryMenager &tmp_manager = DeviceTemporaryMemoryMenager::get_instance();
    double *_tmp0, *_tmp1, *_tmp2;
    double *d_buffer0 = (double*)tmp_manager.get_mem(216 * num_elements * sizeof(double));
    double *d_buffer1 = (double*)tmp_manager.get_mem(216 * num_elements * sizeof(double));
    _tmp0 = d_buffer0;
    {
      unsigned offsets_rT = 0;
      unsigned offsets_I = I_offset;
      unsigned offsets__tmp0 = 216;


      device_gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 24, 9, 56, 1.0, rT(2), 24, I, 56, 0.0, _tmp0, 24, offsets_rT, offsets_I, offsets__tmp0, num_elements);
    }
    _tmp1 = d_buffer1;
    {
      unsigned offsets_fP = 0;
      unsigned offsets__tmp0 = 216;
      unsigned offsets__tmp1 = 216;


      device_gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 24, 9, 21, 1.0, fP(1), 24, _tmp0, 24, 0.0, _tmp1, 24, offsets_fP, offsets__tmp0, offsets__tmp1, num_elements);
    }
    _tmp2 = d_buffer0;
    {
      unsigned offsets__tmp1 = 216;
      unsigned offsets_AminusT = AminusT_offset;
      unsigned offsets__tmp2 = 216;


      device_gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 24, 9, 9, 1.0, _tmp1, 24, AminusT, 9, 0.0, _tmp2, 24, offsets__tmp1, offsets_AminusT, offsets__tmp2, num_elements);
    }
    {
      unsigned offsets_rDivM = 0;
      unsigned offsets__tmp2 = 216;
      unsigned offsets_Q = Q_offset;


      device_gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 56, 9, 21, 1.0, rDivM(3), 56, _tmp2, 24, 1.0, Q, 56, offsets_rDivM, offsets__tmp2, offsets_Q, num_elements);
    }
    tmp_manager.free();
    tmp_manager.free();
  }
  void kernel::neighboringFlux::execute44() {
    assert(AminusT != nullptr);
    assert(I != nullptr);
    assert(Q != nullptr);
    assert(fP(2) != nullptr);
    assert(rDivM(3) != nullptr);
    assert(rT(2) != nullptr);
    assert(kernel::neighboringFlux::num_elements != 0);
    // get device temporary memory menager
    DeviceTemporaryMemoryMenager &tmp_manager = DeviceTemporaryMemoryMenager::get_instance();
    double *_tmp0, *_tmp1, *_tmp2;
    double *d_buffer0 = (double*)tmp_manager.get_mem(216 * num_elements * sizeof(double));
    double *d_buffer1 = (double*)tmp_manager.get_mem(216 * num_elements * sizeof(double));
    _tmp0 = d_buffer0;
    {
      unsigned offsets_rT = 0;
      unsigned offsets_I = I_offset;
      unsigned offsets__tmp0 = 216;


      device_gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 24, 9, 56, 1.0, rT(2), 24, I, 56, 0.0, _tmp0, 24, offsets_rT, offsets_I, offsets__tmp0, num_elements);
    }
    _tmp1 = d_buffer1;
    {
      unsigned offsets_fP = 0;
      unsigned offsets__tmp0 = 216;
      unsigned offsets__tmp1 = 216;


      device_gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 24, 9, 21, 1.0, fP(2), 24, _tmp0, 24, 0.0, _tmp1, 24, offsets_fP, offsets__tmp0, offsets__tmp1, num_elements);
    }
    _tmp2 = d_buffer0;
    {
      unsigned offsets__tmp1 = 216;
      unsigned offsets_AminusT = AminusT_offset;
      unsigned offsets__tmp2 = 216;


      device_gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 24, 9, 9, 1.0, _tmp1, 24, AminusT, 9, 0.0, _tmp2, 24, offsets__tmp1, offsets_AminusT, offsets__tmp2, num_elements);
    }
    {
      unsigned offsets_rDivM = 0;
      unsigned offsets__tmp2 = 216;
      unsigned offsets_Q = Q_offset;


      device_gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 56, 9, 21, 1.0, rDivM(3), 56, _tmp2, 24, 1.0, Q, 56, offsets_rDivM, offsets__tmp2, offsets_Q, num_elements);
    }
    tmp_manager.free();
    tmp_manager.free();
  }
  void kernel::neighboringFlux::execute45() {
    assert(AminusT != nullptr);
    assert(I != nullptr);
    assert(Q != nullptr);
    assert(fP(0) != nullptr);
    assert(rDivM(3) != nullptr);
    assert(rT(3) != nullptr);
    assert(kernel::neighboringFlux::num_elements != 0);
    // get device temporary memory menager
    DeviceTemporaryMemoryMenager &tmp_manager = DeviceTemporaryMemoryMenager::get_instance();
    double *_tmp0, *_tmp1, *_tmp2;
    double *d_buffer0 = (double*)tmp_manager.get_mem(216 * num_elements * sizeof(double));
    double *d_buffer1 = (double*)tmp_manager.get_mem(216 * num_elements * sizeof(double));
    _tmp0 = d_buffer0;
    {
      unsigned offsets_rT = 0;
      unsigned offsets_I = I_offset;
      unsigned offsets__tmp0 = 216;


      device_gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 24, 9, 56, 1.0, rT(3), 24, I, 56, 0.0, _tmp0, 24, offsets_rT, offsets_I, offsets__tmp0, num_elements);
    }
    _tmp1 = d_buffer1;
    {
      unsigned offsets_fP = 0;
      unsigned offsets__tmp0 = 216;
      unsigned offsets__tmp1 = 216;


      device_gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 24, 9, 21, 1.0, fP(0), 24, _tmp0, 24, 0.0, _tmp1, 24, offsets_fP, offsets__tmp0, offsets__tmp1, num_elements);
    }
    _tmp2 = d_buffer0;
    {
      unsigned offsets__tmp1 = 216;
      unsigned offsets_AminusT = AminusT_offset;
      unsigned offsets__tmp2 = 216;


      device_gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 24, 9, 9, 1.0, _tmp1, 24, AminusT, 9, 0.0, _tmp2, 24, offsets__tmp1, offsets_AminusT, offsets__tmp2, num_elements);
    }
    {
      unsigned offsets_rDivM = 0;
      unsigned offsets__tmp2 = 216;
      unsigned offsets_Q = Q_offset;


      device_gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 56, 9, 21, 1.0, rDivM(3), 56, _tmp2, 24, 1.0, Q, 56, offsets_rDivM, offsets__tmp2, offsets_Q, num_elements);
    }
    tmp_manager.free();
    tmp_manager.free();
  }
  void kernel::neighboringFlux::execute46() {
    assert(AminusT != nullptr);
    assert(I != nullptr);
    assert(Q != nullptr);
    assert(fP(1) != nullptr);
    assert(rDivM(3) != nullptr);
    assert(rT(3) != nullptr);
    assert(kernel::neighboringFlux::num_elements != 0);
    // get device temporary memory menager
    DeviceTemporaryMemoryMenager &tmp_manager = DeviceTemporaryMemoryMenager::get_instance();
    double *_tmp0, *_tmp1, *_tmp2;
    double *d_buffer0 = (double*)tmp_manager.get_mem(216 * num_elements * sizeof(double));
    double *d_buffer1 = (double*)tmp_manager.get_mem(216 * num_elements * sizeof(double));
    _tmp0 = d_buffer0;
    {
      unsigned offsets_rT = 0;
      unsigned offsets_I = I_offset;
      unsigned offsets__tmp0 = 216;


      device_gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 24, 9, 56, 1.0, rT(3), 24, I, 56, 0.0, _tmp0, 24, offsets_rT, offsets_I, offsets__tmp0, num_elements);
    }
    _tmp1 = d_buffer1;
    {
      unsigned offsets_fP = 0;
      unsigned offsets__tmp0 = 216;
      unsigned offsets__tmp1 = 216;


      device_gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 24, 9, 21, 1.0, fP(1), 24, _tmp0, 24, 0.0, _tmp1, 24, offsets_fP, offsets__tmp0, offsets__tmp1, num_elements);
    }
    _tmp2 = d_buffer0;
    {
      unsigned offsets__tmp1 = 216;
      unsigned offsets_AminusT = AminusT_offset;
      unsigned offsets__tmp2 = 216;


      device_gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 24, 9, 9, 1.0, _tmp1, 24, AminusT, 9, 0.0, _tmp2, 24, offsets__tmp1, offsets_AminusT, offsets__tmp2, num_elements);
    }
    {
      unsigned offsets_rDivM = 0;
      unsigned offsets__tmp2 = 216;
      unsigned offsets_Q = Q_offset;


      device_gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 56, 9, 21, 1.0, rDivM(3), 56, _tmp2, 24, 1.0, Q, 56, offsets_rDivM, offsets__tmp2, offsets_Q, num_elements);
    }
    tmp_manager.free();
    tmp_manager.free();
  }
  void kernel::neighboringFlux::execute47() {
    assert(AminusT != nullptr);
    assert(I != nullptr);
    assert(Q != nullptr);
    assert(fP(2) != nullptr);
    assert(rDivM(3) != nullptr);
    assert(rT(3) != nullptr);
    assert(kernel::neighboringFlux::num_elements != 0);
    // get device temporary memory menager
    DeviceTemporaryMemoryMenager &tmp_manager = DeviceTemporaryMemoryMenager::get_instance();
    double *_tmp0, *_tmp1, *_tmp2;
    double *d_buffer0 = (double*)tmp_manager.get_mem(216 * num_elements * sizeof(double));
    double *d_buffer1 = (double*)tmp_manager.get_mem(216 * num_elements * sizeof(double));
    _tmp0 = d_buffer0;
    {
      unsigned offsets_rT = 0;
      unsigned offsets_I = I_offset;
      unsigned offsets__tmp0 = 216;


      device_gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 24, 9, 56, 1.0, rT(3), 24, I, 56, 0.0, _tmp0, 24, offsets_rT, offsets_I, offsets__tmp0, num_elements);
    }
    _tmp1 = d_buffer1;
    {
      unsigned offsets_fP = 0;
      unsigned offsets__tmp0 = 216;
      unsigned offsets__tmp1 = 216;


      device_gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 24, 9, 21, 1.0, fP(2), 24, _tmp0, 24, 0.0, _tmp1, 24, offsets_fP, offsets__tmp0, offsets__tmp1, num_elements);
    }
    _tmp2 = d_buffer0;
    {
      unsigned offsets__tmp1 = 216;
      unsigned offsets_AminusT = AminusT_offset;
      unsigned offsets__tmp2 = 216;


      device_gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 24, 9, 9, 1.0, _tmp1, 24, AminusT, 9, 0.0, _tmp2, 24, offsets__tmp1, offsets_AminusT, offsets__tmp2, num_elements);
    }
    {
      unsigned offsets_rDivM = 0;
      unsigned offsets__tmp2 = 216;
      unsigned offsets_Q = Q_offset;


      device_gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 56, 9, 21, 1.0, rDivM(3), 56, _tmp2, 24, 1.0, Q, 56, offsets_rDivM, offsets__tmp2, offsets_Q, num_elements);
    }
    tmp_manager.free();
    tmp_manager.free();
  }
  constexpr unsigned long const kernel::derivativeTaylorExpansion::NonZeroFlops[];
  constexpr unsigned long const kernel::derivativeTaylorExpansion::HardwareFlops[];
  constexpr kernel::derivativeTaylorExpansion::member_function_ptr kernel::derivativeTaylorExpansion::ExecutePtrs[];
  void kernel::derivativeTaylorExpansion::execute0() {
    assert(!std::isnan(power));
    assert(I != nullptr);
    assert(dQ(0) != nullptr);
    assert(kernel::derivativeTaylorExpansion::num_elements != 0);
    device_copy_add_scale(56, 9, power, dQ(0), 56, 0.0, I, 56, dQ_offset(0) + 0, I_offset + 0, num_elements);
  }
  void kernel::derivativeTaylorExpansion::execute1() {
    assert(!std::isnan(power));
    assert(I != nullptr);
    assert(dQ(1) != nullptr);
    assert(kernel::derivativeTaylorExpansion::num_elements != 0);
    device_copy_add_scale(35, 9, power, dQ(1), 36, 1.0, I, 56, dQ_offset(1) + 0, I_offset + 0, num_elements);
  }
  void kernel::derivativeTaylorExpansion::execute2() {
    assert(!std::isnan(power));
    assert(I != nullptr);
    assert(dQ(2) != nullptr);
    assert(kernel::derivativeTaylorExpansion::num_elements != 0);
    device_copy_add_scale(20, 9, power, dQ(2), 20, 1.0, I, 56, dQ_offset(2) + 0, I_offset + 0, num_elements);
  }
  void kernel::derivativeTaylorExpansion::execute3() {
    assert(!std::isnan(power));
    assert(I != nullptr);
    assert(dQ(3) != nullptr);
    assert(kernel::derivativeTaylorExpansion::num_elements != 0);
    device_copy_add_scale(10, 9, power, dQ(3), 12, 1.0, I, 56, dQ_offset(3) + 0, I_offset + 0, num_elements);
  }
  void kernel::derivativeTaylorExpansion::execute4() {
    assert(!std::isnan(power));
    assert(I != nullptr);
    assert(dQ(4) != nullptr);
    assert(kernel::derivativeTaylorExpansion::num_elements != 0);
    device_copy_add_scale(4, 9, power, dQ(4), 4, 1.0, I, 56, dQ_offset(4) + 0, I_offset + 0, num_elements);
  }
  void kernel::derivativeTaylorExpansion::execute5() {
    assert(!std::isnan(power));
    assert(I != nullptr);
    assert(dQ(5) != nullptr);
    assert(kernel::derivativeTaylorExpansion::num_elements != 0);
    device_copy_add_scale(1, 9, power, dQ(5), 4, 1.0, I, 56, dQ_offset(5) + 0, I_offset + 0, num_elements);
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
    // get device temporary memory menager
    DeviceTemporaryMemoryMenager &tmp_manager = DeviceTemporaryMemoryMenager::get_instance();
    double *_tmp0, *_tmp2, *_tmp4;
    double *d_buffer0 = (double*)tmp_manager.get_mem(324 * num_elements * sizeof(double));
    _tmp0 = d_buffer0;
    {
    unsigned offsets_kDivMT = 0;
    unsigned offsets_dQ = dQ_offset(0);
    unsigned offsets__tmp0 = 324;
    
    
    device_gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 36, 9, 53, 1.0, kDivMT(0), 36, dQ(0), 56, 0.0, _tmp0, 36, offsets_kDivMT, offsets_dQ + 1, offsets__tmp0, num_elements);
    }
    {
    unsigned offsets__tmp0 = 324;
    unsigned offsets_star = star_offset(0);
    unsigned offsets_dQ = dQ_offset(1);
    
    
    device_gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 36, 9, 9, 1.0, _tmp0, 36, star(0), 9, 0.0, dQ(1), 36, offsets__tmp0, offsets_star, offsets_dQ, num_elements);
    }
    _tmp2 = d_buffer0;
    {
    unsigned offsets_kDivMT = 0;
    unsigned offsets_dQ = dQ_offset(0);
    unsigned offsets__tmp2 = 324;
    
    
    device_gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 36, 9, 54, 1.0, kDivMT(1), 36, dQ(0), 56, 0.0, _tmp2, 36, offsets_kDivMT, offsets_dQ + 1, offsets__tmp2, num_elements);
    }
    {
    unsigned offsets__tmp2 = 324;
    unsigned offsets_star = star_offset(1);
    unsigned offsets_dQ = dQ_offset(1);
    
    
    device_gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 36, 9, 9, 1.0, _tmp2, 36, star(1), 9, 1.0, dQ(1), 36, offsets__tmp2, offsets_star, offsets_dQ, num_elements);
    }
    _tmp4 = d_buffer0;
    {
    unsigned offsets_kDivMT = 0;
    unsigned offsets_dQ = dQ_offset(0);
    unsigned offsets__tmp4 = 324;
    
    
    device_gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 36, 9, 55, 1.0, kDivMT(2), 36, dQ(0), 56, 0.0, _tmp4, 36, offsets_kDivMT, offsets_dQ + 1, offsets__tmp4, num_elements);
    }
    {
    unsigned offsets__tmp4 = 324;
    unsigned offsets_star = star_offset(2);
    unsigned offsets_dQ = dQ_offset(1);
    
    
    device_gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 36, 9, 9, 1.0, _tmp4, 36, star(2), 9, 1.0, dQ(1), 36, offsets__tmp4, offsets_star, offsets_dQ, num_elements);
    }
    tmp_manager.free();
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
    // get device temporary memory menager
    DeviceTemporaryMemoryMenager &tmp_manager = DeviceTemporaryMemoryMenager::get_instance();
    double *_tmp0, *_tmp2, *_tmp4;
    double *d_buffer0 = (double*)tmp_manager.get_mem(180 * num_elements * sizeof(double));
    _tmp0 = d_buffer0;
    {
    unsigned offsets_kDivMT = 0;
    unsigned offsets_dQ = dQ_offset(1);
    unsigned offsets__tmp0 = 180;
    
    
    device_gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 20, 9, 32, 1.0, kDivMT(0), 36, dQ(1), 36, 0.0, _tmp0, 20, offsets_kDivMT, offsets_dQ + 1, offsets__tmp0, num_elements);
    }
    {
    unsigned offsets__tmp0 = 180;
    unsigned offsets_star = star_offset(0);
    unsigned offsets_dQ = dQ_offset(2);
    
    
    device_gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 20, 9, 9, 1.0, _tmp0, 20, star(0), 9, 0.0, dQ(2), 20, offsets__tmp0, offsets_star, offsets_dQ, num_elements);
    }
    _tmp2 = d_buffer0;
    {
    unsigned offsets_kDivMT = 0;
    unsigned offsets_dQ = dQ_offset(1);
    unsigned offsets__tmp2 = 180;
    
    
    device_gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 20, 9, 33, 1.0, kDivMT(1), 36, dQ(1), 36, 0.0, _tmp2, 20, offsets_kDivMT, offsets_dQ + 1, offsets__tmp2, num_elements);
    }
    {
    unsigned offsets__tmp2 = 180;
    unsigned offsets_star = star_offset(1);
    unsigned offsets_dQ = dQ_offset(2);
    
    
    device_gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 20, 9, 9, 1.0, _tmp2, 20, star(1), 9, 1.0, dQ(2), 20, offsets__tmp2, offsets_star, offsets_dQ, num_elements);
    }
    _tmp4 = d_buffer0;
    {
    unsigned offsets_kDivMT = 0;
    unsigned offsets_dQ = dQ_offset(1);
    unsigned offsets__tmp4 = 180;
    
    
    device_gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 20, 9, 34, 1.0, kDivMT(2), 36, dQ(1), 36, 0.0, _tmp4, 20, offsets_kDivMT, offsets_dQ + 1, offsets__tmp4, num_elements);
    }
    {
    unsigned offsets__tmp4 = 180;
    unsigned offsets_star = star_offset(2);
    unsigned offsets_dQ = dQ_offset(2);
    
    
    device_gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 20, 9, 9, 1.0, _tmp4, 20, star(2), 9, 1.0, dQ(2), 20, offsets__tmp4, offsets_star, offsets_dQ, num_elements);
    }
    tmp_manager.free();
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
    // get device temporary memory menager
    DeviceTemporaryMemoryMenager &tmp_manager = DeviceTemporaryMemoryMenager::get_instance();
    double *_tmp0, *_tmp2, *_tmp4;
    double *d_buffer0 = (double*)tmp_manager.get_mem(108 * num_elements * sizeof(double));
    _tmp0 = d_buffer0;
    {
    unsigned offsets_kDivMT = 0;
    unsigned offsets_dQ = dQ_offset(2);
    unsigned offsets__tmp0 = 108;
    
    
    device_gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 12, 9, 17, 1.0, kDivMT(0), 36, dQ(2), 20, 0.0, _tmp0, 12, offsets_kDivMT, offsets_dQ + 1, offsets__tmp0, num_elements);
    }
    {
    unsigned offsets__tmp0 = 108;
    unsigned offsets_star = star_offset(0);
    unsigned offsets_dQ = dQ_offset(3);
    
    
    device_gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 12, 9, 9, 1.0, _tmp0, 12, star(0), 9, 0.0, dQ(3), 12, offsets__tmp0, offsets_star, offsets_dQ, num_elements);
    }
    _tmp2 = d_buffer0;
    {
    unsigned offsets_kDivMT = 0;
    unsigned offsets_dQ = dQ_offset(2);
    unsigned offsets__tmp2 = 108;
    
    
    device_gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 12, 9, 18, 1.0, kDivMT(1), 36, dQ(2), 20, 0.0, _tmp2, 12, offsets_kDivMT, offsets_dQ + 1, offsets__tmp2, num_elements);
    }
    {
    unsigned offsets__tmp2 = 108;
    unsigned offsets_star = star_offset(1);
    unsigned offsets_dQ = dQ_offset(3);
    
    
    device_gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 12, 9, 9, 1.0, _tmp2, 12, star(1), 9, 1.0, dQ(3), 12, offsets__tmp2, offsets_star, offsets_dQ, num_elements);
    }
    _tmp4 = d_buffer0;
    {
    unsigned offsets_kDivMT = 0;
    unsigned offsets_dQ = dQ_offset(2);
    unsigned offsets__tmp4 = 108;
    
    
    device_gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 12, 9, 19, 1.0, kDivMT(2), 36, dQ(2), 20, 0.0, _tmp4, 12, offsets_kDivMT, offsets_dQ + 1, offsets__tmp4, num_elements);
    }
    {
    unsigned offsets__tmp4 = 108;
    unsigned offsets_star = star_offset(2);
    unsigned offsets_dQ = dQ_offset(3);
    
    
    device_gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 12, 9, 9, 1.0, _tmp4, 12, star(2), 9, 1.0, dQ(3), 12, offsets__tmp4, offsets_star, offsets_dQ, num_elements);
    }
    tmp_manager.free();
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
    // get device temporary memory menager
    DeviceTemporaryMemoryMenager &tmp_manager = DeviceTemporaryMemoryMenager::get_instance();
    double *_tmp0, *_tmp2, *_tmp4;
    double *d_buffer0 = (double*)tmp_manager.get_mem(36 * num_elements * sizeof(double));
    _tmp0 = d_buffer0;
    {
    unsigned offsets_kDivMT = 0;
    unsigned offsets_dQ = dQ_offset(3);
    unsigned offsets__tmp0 = 36;
    
    
    device_gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 4, 9, 7, 1.0, kDivMT(0), 36, dQ(3), 12, 0.0, _tmp0, 4, offsets_kDivMT, offsets_dQ + 1, offsets__tmp0, num_elements);
    }
    {
    unsigned offsets__tmp0 = 36;
    unsigned offsets_star = star_offset(0);
    unsigned offsets_dQ = dQ_offset(4);
    
    
    device_gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 4, 9, 9, 1.0, _tmp0, 4, star(0), 9, 0.0, dQ(4), 4, offsets__tmp0, offsets_star, offsets_dQ, num_elements);
    }
    _tmp2 = d_buffer0;
    {
    unsigned offsets_kDivMT = 0;
    unsigned offsets_dQ = dQ_offset(3);
    unsigned offsets__tmp2 = 36;
    
    
    device_gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 4, 9, 8, 1.0, kDivMT(1), 36, dQ(3), 12, 0.0, _tmp2, 4, offsets_kDivMT, offsets_dQ + 1, offsets__tmp2, num_elements);
    }
    {
    unsigned offsets__tmp2 = 36;
    unsigned offsets_star = star_offset(1);
    unsigned offsets_dQ = dQ_offset(4);
    
    
    device_gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 4, 9, 9, 1.0, _tmp2, 4, star(1), 9, 1.0, dQ(4), 4, offsets__tmp2, offsets_star, offsets_dQ, num_elements);
    }
    _tmp4 = d_buffer0;
    {
    unsigned offsets_kDivMT = 0;
    unsigned offsets_dQ = dQ_offset(3);
    unsigned offsets__tmp4 = 36;
    
    
    device_gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 4, 9, 9, 1.0, kDivMT(2), 36, dQ(3), 12, 0.0, _tmp4, 4, offsets_kDivMT, offsets_dQ + 1, offsets__tmp4, num_elements);
    }
    {
    unsigned offsets__tmp4 = 36;
    unsigned offsets_star = star_offset(2);
    unsigned offsets_dQ = dQ_offset(4);
    
    
    device_gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 4, 9, 9, 1.0, _tmp4, 4, star(2), 9, 1.0, dQ(4), 4, offsets__tmp4, offsets_star, offsets_dQ, num_elements);
    }
    tmp_manager.free();
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
    // get device temporary memory menager
    DeviceTemporaryMemoryMenager &tmp_manager = DeviceTemporaryMemoryMenager::get_instance();
    double *_tmp0, *_tmp2, *_tmp4;
    double *d_buffer0 = (double*)tmp_manager.get_mem(36 * num_elements * sizeof(double));
    _tmp0 = d_buffer0;
    {
    unsigned offsets_kDivMT = 0;
    unsigned offsets_dQ = dQ_offset(4);
    unsigned offsets__tmp0 = 36;
    
    
    device_gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 4, 9, 1, 1.0, kDivMT(0), 36, dQ(4), 4, 0.0, _tmp0, 4, offsets_kDivMT, offsets_dQ + 1, offsets__tmp0, num_elements);
    }
    {
    unsigned offsets__tmp0 = 36;
    unsigned offsets_star = star_offset(0);
    unsigned offsets_dQ = dQ_offset(5);
    
    
    device_gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 4, 9, 9, 1.0, _tmp0, 4, star(0), 9, 0.0, dQ(5), 4, offsets__tmp0, offsets_star, offsets_dQ, num_elements);
    }
    _tmp2 = d_buffer0;
    {
    unsigned offsets_kDivMT = 0;
    unsigned offsets_dQ = dQ_offset(4);
    unsigned offsets__tmp2 = 36;
    
    
    device_gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 4, 9, 2, 1.0, kDivMT(1), 36, dQ(4), 4, 0.0, _tmp2, 4, offsets_kDivMT, offsets_dQ + 1, offsets__tmp2, num_elements);
    }
    {
    unsigned offsets__tmp2 = 36;
    unsigned offsets_star = star_offset(1);
    unsigned offsets_dQ = dQ_offset(5);
    
    
    device_gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 4, 9, 9, 1.0, _tmp2, 4, star(1), 9, 1.0, dQ(5), 4, offsets__tmp2, offsets_star, offsets_dQ, num_elements);
    }
    _tmp4 = d_buffer0;
    {
    unsigned offsets_kDivMT = 0;
    unsigned offsets_dQ = dQ_offset(4);
    unsigned offsets__tmp4 = 36;
    
    
    device_gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 4, 9, 3, 1.0, kDivMT(2), 36, dQ(4), 4, 0.0, _tmp4, 4, offsets_kDivMT, offsets_dQ + 1, offsets__tmp4, num_elements);
    }
    {
    unsigned offsets__tmp4 = 36;
    unsigned offsets_star = star_offset(2);
    unsigned offsets_dQ = dQ_offset(5);
    
    
    device_gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 4, 9, 9, 1.0, _tmp4, 4, star(2), 9, 1.0, dQ(5), 4, offsets__tmp4, offsets_star, offsets_dQ, num_elements);
    }
    tmp_manager.free();
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
    // get device temporary memory menager
    DeviceTemporaryMemoryMenager &tmp_manager = DeviceTemporaryMemoryMenager::get_instance();
    double *_tmp0;
    double *d_buffer0 = (double*)tmp_manager.get_mem(468 * num_elements * sizeof(double));
    _tmp0 = d_buffer0;
    {
      unsigned offsets_V3mTo2n = 0;
      unsigned offsets_Q = Q_offset;
      unsigned offsets__tmp0 = 468;


      device_gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 52, 9, 56, 1.0, V3mTo2n(0,0), 52, Q, 56, 0.0, _tmp0, 52, offsets_V3mTo2n, offsets_Q, offsets__tmp0, num_elements);
    }
    {
      unsigned offsets__tmp0 = 468;
      unsigned offsets_godunovMatrix = godunovMatrix_offset;
      unsigned offsets_godunovState = godunovState_offset;


      device_gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 52, 9, 9, 1.0, _tmp0, 52, godunovMatrix, 9, 0.0, godunovState, 52, offsets__tmp0, offsets_godunovMatrix, offsets_godunovState, num_elements);
    }
    tmp_manager.free();
  }
  void kernel::godunovState::execute1() {
    assert(Q != nullptr);
    assert(V3mTo2n(1,0) != nullptr);
    assert(godunovMatrix != nullptr);
    assert(godunovState != nullptr);
    assert(kernel::godunovState::num_elements != 0);
    // get device temporary memory menager
    DeviceTemporaryMemoryMenager &tmp_manager = DeviceTemporaryMemoryMenager::get_instance();
    double *_tmp0;
    double *d_buffer0 = (double*)tmp_manager.get_mem(468 * num_elements * sizeof(double));
    _tmp0 = d_buffer0;
    {
      unsigned offsets_V3mTo2n = 0;
      unsigned offsets_Q = Q_offset;
      unsigned offsets__tmp0 = 468;


      device_gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 52, 9, 56, 1.0, V3mTo2n(1,0), 52, Q, 56, 0.0, _tmp0, 52, offsets_V3mTo2n, offsets_Q, offsets__tmp0, num_elements);
    }
    {
      unsigned offsets__tmp0 = 468;
      unsigned offsets_godunovMatrix = godunovMatrix_offset;
      unsigned offsets_godunovState = godunovState_offset;


      device_gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 52, 9, 9, 1.0, _tmp0, 52, godunovMatrix, 9, 0.0, godunovState, 52, offsets__tmp0, offsets_godunovMatrix, offsets_godunovState, num_elements);
    }
    tmp_manager.free();
  }
  void kernel::godunovState::execute2() {
    assert(Q != nullptr);
    assert(V3mTo2n(2,0) != nullptr);
    assert(godunovMatrix != nullptr);
    assert(godunovState != nullptr);
    assert(kernel::godunovState::num_elements != 0);
    // get device temporary memory menager
    DeviceTemporaryMemoryMenager &tmp_manager = DeviceTemporaryMemoryMenager::get_instance();
    double *_tmp0;
    double *d_buffer0 = (double*)tmp_manager.get_mem(468 * num_elements * sizeof(double));
    _tmp0 = d_buffer0;
    {
      unsigned offsets_V3mTo2n = 0;
      unsigned offsets_Q = Q_offset;
      unsigned offsets__tmp0 = 468;


      device_gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 52, 9, 56, 1.0, V3mTo2n(2,0), 52, Q, 56, 0.0, _tmp0, 52, offsets_V3mTo2n, offsets_Q, offsets__tmp0, num_elements);
    }
    {
      unsigned offsets__tmp0 = 468;
      unsigned offsets_godunovMatrix = godunovMatrix_offset;
      unsigned offsets_godunovState = godunovState_offset;


      device_gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 52, 9, 9, 1.0, _tmp0, 52, godunovMatrix, 9, 0.0, godunovState, 52, offsets__tmp0, offsets_godunovMatrix, offsets_godunovState, num_elements);
    }
    tmp_manager.free();
  }
  void kernel::godunovState::execute3() {
    assert(Q != nullptr);
    assert(V3mTo2n(3,0) != nullptr);
    assert(godunovMatrix != nullptr);
    assert(godunovState != nullptr);
    assert(kernel::godunovState::num_elements != 0);
    // get device temporary memory menager
    DeviceTemporaryMemoryMenager &tmp_manager = DeviceTemporaryMemoryMenager::get_instance();
    double *_tmp0;
    double *d_buffer0 = (double*)tmp_manager.get_mem(468 * num_elements * sizeof(double));
    _tmp0 = d_buffer0;
    {
      unsigned offsets_V3mTo2n = 0;
      unsigned offsets_Q = Q_offset;
      unsigned offsets__tmp0 = 468;


      device_gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 52, 9, 56, 1.0, V3mTo2n(3,0), 52, Q, 56, 0.0, _tmp0, 52, offsets_V3mTo2n, offsets_Q, offsets__tmp0, num_elements);
    }
    {
      unsigned offsets__tmp0 = 468;
      unsigned offsets_godunovMatrix = godunovMatrix_offset;
      unsigned offsets_godunovState = godunovState_offset;


      device_gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 52, 9, 9, 1.0, _tmp0, 52, godunovMatrix, 9, 0.0, godunovState, 52, offsets__tmp0, offsets_godunovMatrix, offsets_godunovState, num_elements);
    }
    tmp_manager.free();
  }
  void kernel::godunovState::execute4() {
    assert(Q != nullptr);
    assert(V3mTo2n(0,1) != nullptr);
    assert(godunovMatrix != nullptr);
    assert(godunovState != nullptr);
    assert(kernel::godunovState::num_elements != 0);
    // get device temporary memory menager
    DeviceTemporaryMemoryMenager &tmp_manager = DeviceTemporaryMemoryMenager::get_instance();
    double *_tmp0;
    double *d_buffer0 = (double*)tmp_manager.get_mem(468 * num_elements * sizeof(double));
    _tmp0 = d_buffer0;
    {
      unsigned offsets_V3mTo2n = 0;
      unsigned offsets_Q = Q_offset;
      unsigned offsets__tmp0 = 468;


      device_gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 52, 9, 56, 1.0, V3mTo2n(0,1), 52, Q, 56, 0.0, _tmp0, 52, offsets_V3mTo2n, offsets_Q, offsets__tmp0, num_elements);
    }
    {
      unsigned offsets__tmp0 = 468;
      unsigned offsets_godunovMatrix = godunovMatrix_offset;
      unsigned offsets_godunovState = godunovState_offset;


      device_gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 52, 9, 9, 1.0, _tmp0, 52, godunovMatrix, 9, 1.0, godunovState, 52, offsets__tmp0, offsets_godunovMatrix, offsets_godunovState, num_elements);
    }
    tmp_manager.free();
  }
  void kernel::godunovState::execute5() {
    assert(Q != nullptr);
    assert(V3mTo2n(1,1) != nullptr);
    assert(godunovMatrix != nullptr);
    assert(godunovState != nullptr);
    assert(kernel::godunovState::num_elements != 0);
    // get device temporary memory menager
    DeviceTemporaryMemoryMenager &tmp_manager = DeviceTemporaryMemoryMenager::get_instance();
    double *_tmp0;
    double *d_buffer0 = (double*)tmp_manager.get_mem(468 * num_elements * sizeof(double));
    _tmp0 = d_buffer0;
    {
      unsigned offsets_V3mTo2n = 0;
      unsigned offsets_Q = Q_offset;
      unsigned offsets__tmp0 = 468;


      device_gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 52, 9, 56, 1.0, V3mTo2n(1,1), 52, Q, 56, 0.0, _tmp0, 52, offsets_V3mTo2n, offsets_Q, offsets__tmp0, num_elements);
    }
    {
      unsigned offsets__tmp0 = 468;
      unsigned offsets_godunovMatrix = godunovMatrix_offset;
      unsigned offsets_godunovState = godunovState_offset;


      device_gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 52, 9, 9, 1.0, _tmp0, 52, godunovMatrix, 9, 1.0, godunovState, 52, offsets__tmp0, offsets_godunovMatrix, offsets_godunovState, num_elements);
    }
    tmp_manager.free();
  }
  void kernel::godunovState::execute6() {
    assert(Q != nullptr);
    assert(V3mTo2n(2,1) != nullptr);
    assert(godunovMatrix != nullptr);
    assert(godunovState != nullptr);
    assert(kernel::godunovState::num_elements != 0);
    // get device temporary memory menager
    DeviceTemporaryMemoryMenager &tmp_manager = DeviceTemporaryMemoryMenager::get_instance();
    double *_tmp0;
    double *d_buffer0 = (double*)tmp_manager.get_mem(468 * num_elements * sizeof(double));
    _tmp0 = d_buffer0;
    {
      unsigned offsets_V3mTo2n = 0;
      unsigned offsets_Q = Q_offset;
      unsigned offsets__tmp0 = 468;


      device_gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 52, 9, 56, 1.0, V3mTo2n(2,1), 52, Q, 56, 0.0, _tmp0, 52, offsets_V3mTo2n, offsets_Q, offsets__tmp0, num_elements);
    }
    {
      unsigned offsets__tmp0 = 468;
      unsigned offsets_godunovMatrix = godunovMatrix_offset;
      unsigned offsets_godunovState = godunovState_offset;


      device_gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 52, 9, 9, 1.0, _tmp0, 52, godunovMatrix, 9, 1.0, godunovState, 52, offsets__tmp0, offsets_godunovMatrix, offsets_godunovState, num_elements);
    }
    tmp_manager.free();
  }
  void kernel::godunovState::execute7() {
    assert(Q != nullptr);
    assert(V3mTo2n(3,1) != nullptr);
    assert(godunovMatrix != nullptr);
    assert(godunovState != nullptr);
    assert(kernel::godunovState::num_elements != 0);
    // get device temporary memory menager
    DeviceTemporaryMemoryMenager &tmp_manager = DeviceTemporaryMemoryMenager::get_instance();
    double *_tmp0;
    double *d_buffer0 = (double*)tmp_manager.get_mem(468 * num_elements * sizeof(double));
    _tmp0 = d_buffer0;
    {
      unsigned offsets_V3mTo2n = 0;
      unsigned offsets_Q = Q_offset;
      unsigned offsets__tmp0 = 468;


      device_gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 52, 9, 56, 1.0, V3mTo2n(3,1), 52, Q, 56, 0.0, _tmp0, 52, offsets_V3mTo2n, offsets_Q, offsets__tmp0, num_elements);
    }
    {
      unsigned offsets__tmp0 = 468;
      unsigned offsets_godunovMatrix = godunovMatrix_offset;
      unsigned offsets_godunovState = godunovState_offset;


      device_gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 52, 9, 9, 1.0, _tmp0, 52, godunovMatrix, 9, 1.0, godunovState, 52, offsets__tmp0, offsets_godunovMatrix, offsets_godunovState, num_elements);
    }
    tmp_manager.free();
  }
  void kernel::godunovState::execute8() {
    assert(Q != nullptr);
    assert(V3mTo2n(0,2) != nullptr);
    assert(godunovMatrix != nullptr);
    assert(godunovState != nullptr);
    assert(kernel::godunovState::num_elements != 0);
    // get device temporary memory menager
    DeviceTemporaryMemoryMenager &tmp_manager = DeviceTemporaryMemoryMenager::get_instance();
    double *_tmp0;
    double *d_buffer0 = (double*)tmp_manager.get_mem(468 * num_elements * sizeof(double));
    _tmp0 = d_buffer0;
    {
      unsigned offsets_V3mTo2n = 0;
      unsigned offsets_Q = Q_offset;
      unsigned offsets__tmp0 = 468;


      device_gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 52, 9, 56, 1.0, V3mTo2n(0,2), 52, Q, 56, 0.0, _tmp0, 52, offsets_V3mTo2n, offsets_Q, offsets__tmp0, num_elements);
    }
    {
      unsigned offsets__tmp0 = 468;
      unsigned offsets_godunovMatrix = godunovMatrix_offset;
      unsigned offsets_godunovState = godunovState_offset;


      device_gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 52, 9, 9, 1.0, _tmp0, 52, godunovMatrix, 9, 1.0, godunovState, 52, offsets__tmp0, offsets_godunovMatrix, offsets_godunovState, num_elements);
    }
    tmp_manager.free();
  }
  void kernel::godunovState::execute9() {
    assert(Q != nullptr);
    assert(V3mTo2n(1,2) != nullptr);
    assert(godunovMatrix != nullptr);
    assert(godunovState != nullptr);
    assert(kernel::godunovState::num_elements != 0);
    // get device temporary memory menager
    DeviceTemporaryMemoryMenager &tmp_manager = DeviceTemporaryMemoryMenager::get_instance();
    double *_tmp0;
    double *d_buffer0 = (double*)tmp_manager.get_mem(468 * num_elements * sizeof(double));
    _tmp0 = d_buffer0;
    {
      unsigned offsets_V3mTo2n = 0;
      unsigned offsets_Q = Q_offset;
      unsigned offsets__tmp0 = 468;


      device_gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 52, 9, 56, 1.0, V3mTo2n(1,2), 52, Q, 56, 0.0, _tmp0, 52, offsets_V3mTo2n, offsets_Q, offsets__tmp0, num_elements);
    }
    {
      unsigned offsets__tmp0 = 468;
      unsigned offsets_godunovMatrix = godunovMatrix_offset;
      unsigned offsets_godunovState = godunovState_offset;


      device_gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 52, 9, 9, 1.0, _tmp0, 52, godunovMatrix, 9, 1.0, godunovState, 52, offsets__tmp0, offsets_godunovMatrix, offsets_godunovState, num_elements);
    }
    tmp_manager.free();
  }
  void kernel::godunovState::execute10() {
    assert(Q != nullptr);
    assert(V3mTo2n(2,2) != nullptr);
    assert(godunovMatrix != nullptr);
    assert(godunovState != nullptr);
    assert(kernel::godunovState::num_elements != 0);
    // get device temporary memory menager
    DeviceTemporaryMemoryMenager &tmp_manager = DeviceTemporaryMemoryMenager::get_instance();
    double *_tmp0;
    double *d_buffer0 = (double*)tmp_manager.get_mem(468 * num_elements * sizeof(double));
    _tmp0 = d_buffer0;
    {
      unsigned offsets_V3mTo2n = 0;
      unsigned offsets_Q = Q_offset;
      unsigned offsets__tmp0 = 468;


      device_gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 52, 9, 56, 1.0, V3mTo2n(2,2), 52, Q, 56, 0.0, _tmp0, 52, offsets_V3mTo2n, offsets_Q, offsets__tmp0, num_elements);
    }
    {
      unsigned offsets__tmp0 = 468;
      unsigned offsets_godunovMatrix = godunovMatrix_offset;
      unsigned offsets_godunovState = godunovState_offset;


      device_gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 52, 9, 9, 1.0, _tmp0, 52, godunovMatrix, 9, 1.0, godunovState, 52, offsets__tmp0, offsets_godunovMatrix, offsets_godunovState, num_elements);
    }
    tmp_manager.free();
  }
  void kernel::godunovState::execute11() {
    assert(Q != nullptr);
    assert(V3mTo2n(3,2) != nullptr);
    assert(godunovMatrix != nullptr);
    assert(godunovState != nullptr);
    assert(kernel::godunovState::num_elements != 0);
    // get device temporary memory menager
    DeviceTemporaryMemoryMenager &tmp_manager = DeviceTemporaryMemoryMenager::get_instance();
    double *_tmp0;
    double *d_buffer0 = (double*)tmp_manager.get_mem(468 * num_elements * sizeof(double));
    _tmp0 = d_buffer0;
    {
      unsigned offsets_V3mTo2n = 0;
      unsigned offsets_Q = Q_offset;
      unsigned offsets__tmp0 = 468;


      device_gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 52, 9, 56, 1.0, V3mTo2n(3,2), 52, Q, 56, 0.0, _tmp0, 52, offsets_V3mTo2n, offsets_Q, offsets__tmp0, num_elements);
    }
    {
      unsigned offsets__tmp0 = 468;
      unsigned offsets_godunovMatrix = godunovMatrix_offset;
      unsigned offsets_godunovState = godunovState_offset;


      device_gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 52, 9, 9, 1.0, _tmp0, 52, godunovMatrix, 9, 1.0, godunovState, 52, offsets__tmp0, offsets_godunovMatrix, offsets_godunovState, num_elements);
    }
    tmp_manager.free();
  }
  void kernel::godunovState::execute12() {
    assert(Q != nullptr);
    assert(V3mTo2n(0,3) != nullptr);
    assert(godunovMatrix != nullptr);
    assert(godunovState != nullptr);
    assert(kernel::godunovState::num_elements != 0);
    // get device temporary memory menager
    DeviceTemporaryMemoryMenager &tmp_manager = DeviceTemporaryMemoryMenager::get_instance();
    double *_tmp0;
    double *d_buffer0 = (double*)tmp_manager.get_mem(468 * num_elements * sizeof(double));
    _tmp0 = d_buffer0;
    {
      unsigned offsets_V3mTo2n = 0;
      unsigned offsets_Q = Q_offset;
      unsigned offsets__tmp0 = 468;


      device_gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 52, 9, 56, 1.0, V3mTo2n(0,3), 52, Q, 56, 0.0, _tmp0, 52, offsets_V3mTo2n, offsets_Q, offsets__tmp0, num_elements);
    }
    {
      unsigned offsets__tmp0 = 468;
      unsigned offsets_godunovMatrix = godunovMatrix_offset;
      unsigned offsets_godunovState = godunovState_offset;


      device_gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 52, 9, 9, 1.0, _tmp0, 52, godunovMatrix, 9, 1.0, godunovState, 52, offsets__tmp0, offsets_godunovMatrix, offsets_godunovState, num_elements);
    }
    tmp_manager.free();
  }
  void kernel::godunovState::execute13() {
    assert(Q != nullptr);
    assert(V3mTo2n(1,3) != nullptr);
    assert(godunovMatrix != nullptr);
    assert(godunovState != nullptr);
    assert(kernel::godunovState::num_elements != 0);
    // get device temporary memory menager
    DeviceTemporaryMemoryMenager &tmp_manager = DeviceTemporaryMemoryMenager::get_instance();
    double *_tmp0;
    double *d_buffer0 = (double*)tmp_manager.get_mem(468 * num_elements * sizeof(double));
    _tmp0 = d_buffer0;
    {
      unsigned offsets_V3mTo2n = 0;
      unsigned offsets_Q = Q_offset;
      unsigned offsets__tmp0 = 468;


      device_gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 52, 9, 56, 1.0, V3mTo2n(1,3), 52, Q, 56, 0.0, _tmp0, 52, offsets_V3mTo2n, offsets_Q, offsets__tmp0, num_elements);
    }
    {
      unsigned offsets__tmp0 = 468;
      unsigned offsets_godunovMatrix = godunovMatrix_offset;
      unsigned offsets_godunovState = godunovState_offset;


      device_gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 52, 9, 9, 1.0, _tmp0, 52, godunovMatrix, 9, 1.0, godunovState, 52, offsets__tmp0, offsets_godunovMatrix, offsets_godunovState, num_elements);
    }
    tmp_manager.free();
  }
  void kernel::godunovState::execute14() {
    assert(Q != nullptr);
    assert(V3mTo2n(2,3) != nullptr);
    assert(godunovMatrix != nullptr);
    assert(godunovState != nullptr);
    assert(kernel::godunovState::num_elements != 0);
    // get device temporary memory menager
    DeviceTemporaryMemoryMenager &tmp_manager = DeviceTemporaryMemoryMenager::get_instance();
    double *_tmp0;
    double *d_buffer0 = (double*)tmp_manager.get_mem(468 * num_elements * sizeof(double));
    _tmp0 = d_buffer0;
    {
      unsigned offsets_V3mTo2n = 0;
      unsigned offsets_Q = Q_offset;
      unsigned offsets__tmp0 = 468;


      device_gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 52, 9, 56, 1.0, V3mTo2n(2,3), 52, Q, 56, 0.0, _tmp0, 52, offsets_V3mTo2n, offsets_Q, offsets__tmp0, num_elements);
    }
    {
      unsigned offsets__tmp0 = 468;
      unsigned offsets_godunovMatrix = godunovMatrix_offset;
      unsigned offsets_godunovState = godunovState_offset;


      device_gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 52, 9, 9, 1.0, _tmp0, 52, godunovMatrix, 9, 1.0, godunovState, 52, offsets__tmp0, offsets_godunovMatrix, offsets_godunovState, num_elements);
    }
    tmp_manager.free();
  }
  void kernel::godunovState::execute15() {
    assert(Q != nullptr);
    assert(V3mTo2n(3,3) != nullptr);
    assert(godunovMatrix != nullptr);
    assert(godunovState != nullptr);
    assert(kernel::godunovState::num_elements != 0);
    // get device temporary memory menager
    DeviceTemporaryMemoryMenager &tmp_manager = DeviceTemporaryMemoryMenager::get_instance();
    double *_tmp0;
    double *d_buffer0 = (double*)tmp_manager.get_mem(468 * num_elements * sizeof(double));
    _tmp0 = d_buffer0;
    {
      unsigned offsets_V3mTo2n = 0;
      unsigned offsets_Q = Q_offset;
      unsigned offsets__tmp0 = 468;


      device_gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 52, 9, 56, 1.0, V3mTo2n(3,3), 52, Q, 56, 0.0, _tmp0, 52, offsets_V3mTo2n, offsets_Q, offsets__tmp0, num_elements);
    }
    {
      unsigned offsets__tmp0 = 468;
      unsigned offsets_godunovMatrix = godunovMatrix_offset;
      unsigned offsets_godunovState = godunovState_offset;


      device_gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 52, 9, 9, 1.0, _tmp0, 52, godunovMatrix, 9, 1.0, godunovState, 52, offsets__tmp0, offsets_godunovMatrix, offsets_godunovState, num_elements);
    }
    tmp_manager.free();
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
    // get device temporary memory menager
    DeviceTemporaryMemoryMenager &tmp_manager = DeviceTemporaryMemoryMenager::get_instance();
    double *_tmp0;
    double *d_buffer0 = (double*)tmp_manager.get_mem(468 * num_elements * sizeof(double));
    _tmp0 = d_buffer0;
    {
      unsigned offsets_godunovState = godunovState_offset;
      unsigned offsets_fluxSolver = fluxSolver_offset;
      unsigned offsets__tmp0 = 468;


      device_gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 52, 9, 9, 1.0, godunovState, 52, fluxSolver, 9, 0.0, _tmp0, 52, offsets_godunovState, offsets_fluxSolver, offsets__tmp0, num_elements);
    }
    {
      unsigned offsets_V3mTo2nTWDivM = 0;
      unsigned offsets__tmp0 = 468;
      unsigned offsets_Q = Q_offset;


      device_gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 56, 9, 49, 1.0, V3mTo2nTWDivM(0,0), 56, _tmp0, 52, 1.0, Q, 56, offsets_V3mTo2nTWDivM, offsets__tmp0, offsets_Q, num_elements);
    }
    tmp_manager.free();
  }
  void kernel::nodalFlux::execute1() {
    assert(Q != nullptr);
    assert(V3mTo2nTWDivM(1,0) != nullptr);
    assert(fluxSolver != nullptr);
    assert(godunovState != nullptr);
    assert(kernel::nodalFlux::num_elements != 0);
    // get device temporary memory menager
    DeviceTemporaryMemoryMenager &tmp_manager = DeviceTemporaryMemoryMenager::get_instance();
    double *_tmp0;
    double *d_buffer0 = (double*)tmp_manager.get_mem(468 * num_elements * sizeof(double));
    _tmp0 = d_buffer0;
    {
      unsigned offsets_godunovState = godunovState_offset;
      unsigned offsets_fluxSolver = fluxSolver_offset;
      unsigned offsets__tmp0 = 468;


      device_gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 52, 9, 9, 1.0, godunovState, 52, fluxSolver, 9, 0.0, _tmp0, 52, offsets_godunovState, offsets_fluxSolver, offsets__tmp0, num_elements);
    }
    {
      unsigned offsets_V3mTo2nTWDivM = 0;
      unsigned offsets__tmp0 = 468;
      unsigned offsets_Q = Q_offset;


      device_gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 56, 9, 49, 1.0, V3mTo2nTWDivM(1,0), 56, _tmp0, 52, 1.0, Q, 56, offsets_V3mTo2nTWDivM, offsets__tmp0, offsets_Q, num_elements);
    }
    tmp_manager.free();
  }
  void kernel::nodalFlux::execute2() {
    assert(Q != nullptr);
    assert(V3mTo2nTWDivM(2,0) != nullptr);
    assert(fluxSolver != nullptr);
    assert(godunovState != nullptr);
    assert(kernel::nodalFlux::num_elements != 0);
    // get device temporary memory menager
    DeviceTemporaryMemoryMenager &tmp_manager = DeviceTemporaryMemoryMenager::get_instance();
    double *_tmp0;
    double *d_buffer0 = (double*)tmp_manager.get_mem(468 * num_elements * sizeof(double));
    _tmp0 = d_buffer0;
    {
      unsigned offsets_godunovState = godunovState_offset;
      unsigned offsets_fluxSolver = fluxSolver_offset;
      unsigned offsets__tmp0 = 468;


      device_gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 52, 9, 9, 1.0, godunovState, 52, fluxSolver, 9, 0.0, _tmp0, 52, offsets_godunovState, offsets_fluxSolver, offsets__tmp0, num_elements);
    }
    {
      unsigned offsets_V3mTo2nTWDivM = 0;
      unsigned offsets__tmp0 = 468;
      unsigned offsets_Q = Q_offset;


      device_gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 56, 9, 49, 1.0, V3mTo2nTWDivM(2,0), 56, _tmp0, 52, 1.0, Q, 56, offsets_V3mTo2nTWDivM, offsets__tmp0, offsets_Q, num_elements);
    }
    tmp_manager.free();
  }
  void kernel::nodalFlux::execute3() {
    assert(Q != nullptr);
    assert(V3mTo2nTWDivM(3,0) != nullptr);
    assert(fluxSolver != nullptr);
    assert(godunovState != nullptr);
    assert(kernel::nodalFlux::num_elements != 0);
    // get device temporary memory menager
    DeviceTemporaryMemoryMenager &tmp_manager = DeviceTemporaryMemoryMenager::get_instance();
    double *_tmp0;
    double *d_buffer0 = (double*)tmp_manager.get_mem(468 * num_elements * sizeof(double));
    _tmp0 = d_buffer0;
    {
      unsigned offsets_godunovState = godunovState_offset;
      unsigned offsets_fluxSolver = fluxSolver_offset;
      unsigned offsets__tmp0 = 468;


      device_gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 52, 9, 9, 1.0, godunovState, 52, fluxSolver, 9, 0.0, _tmp0, 52, offsets_godunovState, offsets_fluxSolver, offsets__tmp0, num_elements);
    }
    {
      unsigned offsets_V3mTo2nTWDivM = 0;
      unsigned offsets__tmp0 = 468;
      unsigned offsets_Q = Q_offset;


      device_gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 56, 9, 49, 1.0, V3mTo2nTWDivM(3,0), 56, _tmp0, 52, 1.0, Q, 56, offsets_V3mTo2nTWDivM, offsets__tmp0, offsets_Q, num_elements);
    }
    tmp_manager.free();
  }
  void kernel::nodalFlux::execute4() {
    assert(Q != nullptr);
    assert(V3mTo2nTWDivM(0,1) != nullptr);
    assert(fluxSolver != nullptr);
    assert(godunovState != nullptr);
    assert(kernel::nodalFlux::num_elements != 0);
    // get device temporary memory menager
    DeviceTemporaryMemoryMenager &tmp_manager = DeviceTemporaryMemoryMenager::get_instance();
    double *_tmp0;
    double *d_buffer0 = (double*)tmp_manager.get_mem(468 * num_elements * sizeof(double));
    _tmp0 = d_buffer0;
    {
      unsigned offsets_godunovState = godunovState_offset;
      unsigned offsets_fluxSolver = fluxSolver_offset;
      unsigned offsets__tmp0 = 468;


      device_gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 52, 9, 9, 1.0, godunovState, 52, fluxSolver, 9, 0.0, _tmp0, 52, offsets_godunovState, offsets_fluxSolver, offsets__tmp0, num_elements);
    }
    {
      unsigned offsets_V3mTo2nTWDivM = 0;
      unsigned offsets__tmp0 = 468;
      unsigned offsets_Q = Q_offset;


      device_gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 56, 9, 49, 1.0, V3mTo2nTWDivM(0,1), 56, _tmp0, 52, 1.0, Q, 56, offsets_V3mTo2nTWDivM, offsets__tmp0, offsets_Q, num_elements);
    }
    tmp_manager.free();
  }
  void kernel::nodalFlux::execute5() {
    assert(Q != nullptr);
    assert(V3mTo2nTWDivM(1,1) != nullptr);
    assert(fluxSolver != nullptr);
    assert(godunovState != nullptr);
    assert(kernel::nodalFlux::num_elements != 0);
    // get device temporary memory menager
    DeviceTemporaryMemoryMenager &tmp_manager = DeviceTemporaryMemoryMenager::get_instance();
    double *_tmp0;
    double *d_buffer0 = (double*)tmp_manager.get_mem(468 * num_elements * sizeof(double));
    _tmp0 = d_buffer0;
    {
      unsigned offsets_godunovState = godunovState_offset;
      unsigned offsets_fluxSolver = fluxSolver_offset;
      unsigned offsets__tmp0 = 468;


      device_gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 52, 9, 9, 1.0, godunovState, 52, fluxSolver, 9, 0.0, _tmp0, 52, offsets_godunovState, offsets_fluxSolver, offsets__tmp0, num_elements);
    }
    {
      unsigned offsets_V3mTo2nTWDivM = 0;
      unsigned offsets__tmp0 = 468;
      unsigned offsets_Q = Q_offset;


      device_gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 56, 9, 49, 1.0, V3mTo2nTWDivM(1,1), 56, _tmp0, 52, 1.0, Q, 56, offsets_V3mTo2nTWDivM, offsets__tmp0, offsets_Q, num_elements);
    }
    tmp_manager.free();
  }
  void kernel::nodalFlux::execute6() {
    assert(Q != nullptr);
    assert(V3mTo2nTWDivM(2,1) != nullptr);
    assert(fluxSolver != nullptr);
    assert(godunovState != nullptr);
    assert(kernel::nodalFlux::num_elements != 0);
    // get device temporary memory menager
    DeviceTemporaryMemoryMenager &tmp_manager = DeviceTemporaryMemoryMenager::get_instance();
    double *_tmp0;
    double *d_buffer0 = (double*)tmp_manager.get_mem(468 * num_elements * sizeof(double));
    _tmp0 = d_buffer0;
    {
      unsigned offsets_godunovState = godunovState_offset;
      unsigned offsets_fluxSolver = fluxSolver_offset;
      unsigned offsets__tmp0 = 468;


      device_gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 52, 9, 9, 1.0, godunovState, 52, fluxSolver, 9, 0.0, _tmp0, 52, offsets_godunovState, offsets_fluxSolver, offsets__tmp0, num_elements);
    }
    {
      unsigned offsets_V3mTo2nTWDivM = 0;
      unsigned offsets__tmp0 = 468;
      unsigned offsets_Q = Q_offset;


      device_gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 56, 9, 49, 1.0, V3mTo2nTWDivM(2,1), 56, _tmp0, 52, 1.0, Q, 56, offsets_V3mTo2nTWDivM, offsets__tmp0, offsets_Q, num_elements);
    }
    tmp_manager.free();
  }
  void kernel::nodalFlux::execute7() {
    assert(Q != nullptr);
    assert(V3mTo2nTWDivM(3,1) != nullptr);
    assert(fluxSolver != nullptr);
    assert(godunovState != nullptr);
    assert(kernel::nodalFlux::num_elements != 0);
    // get device temporary memory menager
    DeviceTemporaryMemoryMenager &tmp_manager = DeviceTemporaryMemoryMenager::get_instance();
    double *_tmp0;
    double *d_buffer0 = (double*)tmp_manager.get_mem(468 * num_elements * sizeof(double));
    _tmp0 = d_buffer0;
    {
      unsigned offsets_godunovState = godunovState_offset;
      unsigned offsets_fluxSolver = fluxSolver_offset;
      unsigned offsets__tmp0 = 468;


      device_gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 52, 9, 9, 1.0, godunovState, 52, fluxSolver, 9, 0.0, _tmp0, 52, offsets_godunovState, offsets_fluxSolver, offsets__tmp0, num_elements);
    }
    {
      unsigned offsets_V3mTo2nTWDivM = 0;
      unsigned offsets__tmp0 = 468;
      unsigned offsets_Q = Q_offset;


      device_gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 56, 9, 49, 1.0, V3mTo2nTWDivM(3,1), 56, _tmp0, 52, 1.0, Q, 56, offsets_V3mTo2nTWDivM, offsets__tmp0, offsets_Q, num_elements);
    }
    tmp_manager.free();
  }
  void kernel::nodalFlux::execute8() {
    assert(Q != nullptr);
    assert(V3mTo2nTWDivM(0,2) != nullptr);
    assert(fluxSolver != nullptr);
    assert(godunovState != nullptr);
    assert(kernel::nodalFlux::num_elements != 0);
    // get device temporary memory menager
    DeviceTemporaryMemoryMenager &tmp_manager = DeviceTemporaryMemoryMenager::get_instance();
    double *_tmp0;
    double *d_buffer0 = (double*)tmp_manager.get_mem(468 * num_elements * sizeof(double));
    _tmp0 = d_buffer0;
    {
      unsigned offsets_godunovState = godunovState_offset;
      unsigned offsets_fluxSolver = fluxSolver_offset;
      unsigned offsets__tmp0 = 468;


      device_gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 52, 9, 9, 1.0, godunovState, 52, fluxSolver, 9, 0.0, _tmp0, 52, offsets_godunovState, offsets_fluxSolver, offsets__tmp0, num_elements);
    }
    {
      unsigned offsets_V3mTo2nTWDivM = 0;
      unsigned offsets__tmp0 = 468;
      unsigned offsets_Q = Q_offset;


      device_gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 56, 9, 49, 1.0, V3mTo2nTWDivM(0,2), 56, _tmp0, 52, 1.0, Q, 56, offsets_V3mTo2nTWDivM, offsets__tmp0, offsets_Q, num_elements);
    }
    tmp_manager.free();
  }
  void kernel::nodalFlux::execute9() {
    assert(Q != nullptr);
    assert(V3mTo2nTWDivM(1,2) != nullptr);
    assert(fluxSolver != nullptr);
    assert(godunovState != nullptr);
    assert(kernel::nodalFlux::num_elements != 0);
    // get device temporary memory menager
    DeviceTemporaryMemoryMenager &tmp_manager = DeviceTemporaryMemoryMenager::get_instance();
    double *_tmp0;
    double *d_buffer0 = (double*)tmp_manager.get_mem(468 * num_elements * sizeof(double));
    _tmp0 = d_buffer0;
    {
      unsigned offsets_godunovState = godunovState_offset;
      unsigned offsets_fluxSolver = fluxSolver_offset;
      unsigned offsets__tmp0 = 468;


      device_gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 52, 9, 9, 1.0, godunovState, 52, fluxSolver, 9, 0.0, _tmp0, 52, offsets_godunovState, offsets_fluxSolver, offsets__tmp0, num_elements);
    }
    {
      unsigned offsets_V3mTo2nTWDivM = 0;
      unsigned offsets__tmp0 = 468;
      unsigned offsets_Q = Q_offset;


      device_gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 56, 9, 49, 1.0, V3mTo2nTWDivM(1,2), 56, _tmp0, 52, 1.0, Q, 56, offsets_V3mTo2nTWDivM, offsets__tmp0, offsets_Q, num_elements);
    }
    tmp_manager.free();
  }
  void kernel::nodalFlux::execute10() {
    assert(Q != nullptr);
    assert(V3mTo2nTWDivM(2,2) != nullptr);
    assert(fluxSolver != nullptr);
    assert(godunovState != nullptr);
    assert(kernel::nodalFlux::num_elements != 0);
    // get device temporary memory menager
    DeviceTemporaryMemoryMenager &tmp_manager = DeviceTemporaryMemoryMenager::get_instance();
    double *_tmp0;
    double *d_buffer0 = (double*)tmp_manager.get_mem(468 * num_elements * sizeof(double));
    _tmp0 = d_buffer0;
    {
      unsigned offsets_godunovState = godunovState_offset;
      unsigned offsets_fluxSolver = fluxSolver_offset;
      unsigned offsets__tmp0 = 468;


      device_gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 52, 9, 9, 1.0, godunovState, 52, fluxSolver, 9, 0.0, _tmp0, 52, offsets_godunovState, offsets_fluxSolver, offsets__tmp0, num_elements);
    }
    {
      unsigned offsets_V3mTo2nTWDivM = 0;
      unsigned offsets__tmp0 = 468;
      unsigned offsets_Q = Q_offset;


      device_gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 56, 9, 49, 1.0, V3mTo2nTWDivM(2,2), 56, _tmp0, 52, 1.0, Q, 56, offsets_V3mTo2nTWDivM, offsets__tmp0, offsets_Q, num_elements);
    }
    tmp_manager.free();
  }
  void kernel::nodalFlux::execute11() {
    assert(Q != nullptr);
    assert(V3mTo2nTWDivM(3,2) != nullptr);
    assert(fluxSolver != nullptr);
    assert(godunovState != nullptr);
    assert(kernel::nodalFlux::num_elements != 0);
    // get device temporary memory menager
    DeviceTemporaryMemoryMenager &tmp_manager = DeviceTemporaryMemoryMenager::get_instance();
    double *_tmp0;
    double *d_buffer0 = (double*)tmp_manager.get_mem(468 * num_elements * sizeof(double));
    _tmp0 = d_buffer0;
    {
      unsigned offsets_godunovState = godunovState_offset;
      unsigned offsets_fluxSolver = fluxSolver_offset;
      unsigned offsets__tmp0 = 468;


      device_gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 52, 9, 9, 1.0, godunovState, 52, fluxSolver, 9, 0.0, _tmp0, 52, offsets_godunovState, offsets_fluxSolver, offsets__tmp0, num_elements);
    }
    {
      unsigned offsets_V3mTo2nTWDivM = 0;
      unsigned offsets__tmp0 = 468;
      unsigned offsets_Q = Q_offset;


      device_gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 56, 9, 49, 1.0, V3mTo2nTWDivM(3,2), 56, _tmp0, 52, 1.0, Q, 56, offsets_V3mTo2nTWDivM, offsets__tmp0, offsets_Q, num_elements);
    }
    tmp_manager.free();
  }
  void kernel::nodalFlux::execute12() {
    assert(Q != nullptr);
    assert(V3mTo2nTWDivM(0,3) != nullptr);
    assert(fluxSolver != nullptr);
    assert(godunovState != nullptr);
    assert(kernel::nodalFlux::num_elements != 0);
    // get device temporary memory menager
    DeviceTemporaryMemoryMenager &tmp_manager = DeviceTemporaryMemoryMenager::get_instance();
    double *_tmp0;
    double *d_buffer0 = (double*)tmp_manager.get_mem(468 * num_elements * sizeof(double));
    _tmp0 = d_buffer0;
    {
      unsigned offsets_godunovState = godunovState_offset;
      unsigned offsets_fluxSolver = fluxSolver_offset;
      unsigned offsets__tmp0 = 468;


      device_gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 52, 9, 9, 1.0, godunovState, 52, fluxSolver, 9, 0.0, _tmp0, 52, offsets_godunovState, offsets_fluxSolver, offsets__tmp0, num_elements);
    }
    {
      unsigned offsets_V3mTo2nTWDivM = 0;
      unsigned offsets__tmp0 = 468;
      unsigned offsets_Q = Q_offset;


      device_gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 56, 9, 49, 1.0, V3mTo2nTWDivM(0,3), 56, _tmp0, 52, 1.0, Q, 56, offsets_V3mTo2nTWDivM, offsets__tmp0, offsets_Q, num_elements);
    }
    tmp_manager.free();
  }
  void kernel::nodalFlux::execute13() {
    assert(Q != nullptr);
    assert(V3mTo2nTWDivM(1,3) != nullptr);
    assert(fluxSolver != nullptr);
    assert(godunovState != nullptr);
    assert(kernel::nodalFlux::num_elements != 0);
    // get device temporary memory menager
    DeviceTemporaryMemoryMenager &tmp_manager = DeviceTemporaryMemoryMenager::get_instance();
    double *_tmp0;
    double *d_buffer0 = (double*)tmp_manager.get_mem(468 * num_elements * sizeof(double));
    _tmp0 = d_buffer0;
    {
      unsigned offsets_godunovState = godunovState_offset;
      unsigned offsets_fluxSolver = fluxSolver_offset;
      unsigned offsets__tmp0 = 468;


      device_gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 52, 9, 9, 1.0, godunovState, 52, fluxSolver, 9, 0.0, _tmp0, 52, offsets_godunovState, offsets_fluxSolver, offsets__tmp0, num_elements);
    }
    {
      unsigned offsets_V3mTo2nTWDivM = 0;
      unsigned offsets__tmp0 = 468;
      unsigned offsets_Q = Q_offset;


      device_gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 56, 9, 49, 1.0, V3mTo2nTWDivM(1,3), 56, _tmp0, 52, 1.0, Q, 56, offsets_V3mTo2nTWDivM, offsets__tmp0, offsets_Q, num_elements);
    }
    tmp_manager.free();
  }
  void kernel::nodalFlux::execute14() {
    assert(Q != nullptr);
    assert(V3mTo2nTWDivM(2,3) != nullptr);
    assert(fluxSolver != nullptr);
    assert(godunovState != nullptr);
    assert(kernel::nodalFlux::num_elements != 0);
    // get device temporary memory menager
    DeviceTemporaryMemoryMenager &tmp_manager = DeviceTemporaryMemoryMenager::get_instance();
    double *_tmp0;
    double *d_buffer0 = (double*)tmp_manager.get_mem(468 * num_elements * sizeof(double));
    _tmp0 = d_buffer0;
    {
      unsigned offsets_godunovState = godunovState_offset;
      unsigned offsets_fluxSolver = fluxSolver_offset;
      unsigned offsets__tmp0 = 468;


      device_gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 52, 9, 9, 1.0, godunovState, 52, fluxSolver, 9, 0.0, _tmp0, 52, offsets_godunovState, offsets_fluxSolver, offsets__tmp0, num_elements);
    }
    {
      unsigned offsets_V3mTo2nTWDivM = 0;
      unsigned offsets__tmp0 = 468;
      unsigned offsets_Q = Q_offset;


      device_gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 56, 9, 49, 1.0, V3mTo2nTWDivM(2,3), 56, _tmp0, 52, 1.0, Q, 56, offsets_V3mTo2nTWDivM, offsets__tmp0, offsets_Q, num_elements);
    }
    tmp_manager.free();
  }
  void kernel::nodalFlux::execute15() {
    assert(Q != nullptr);
    assert(V3mTo2nTWDivM(3,3) != nullptr);
    assert(fluxSolver != nullptr);
    assert(godunovState != nullptr);
    assert(kernel::nodalFlux::num_elements != 0);
    // get device temporary memory menager
    DeviceTemporaryMemoryMenager &tmp_manager = DeviceTemporaryMemoryMenager::get_instance();
    double *_tmp0;
    double *d_buffer0 = (double*)tmp_manager.get_mem(468 * num_elements * sizeof(double));
    _tmp0 = d_buffer0;
    {
      unsigned offsets_godunovState = godunovState_offset;
      unsigned offsets_fluxSolver = fluxSolver_offset;
      unsigned offsets__tmp0 = 468;


      device_gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 52, 9, 9, 1.0, godunovState, 52, fluxSolver, 9, 0.0, _tmp0, 52, offsets_godunovState, offsets_fluxSolver, offsets__tmp0, num_elements);
    }
    {
      unsigned offsets_V3mTo2nTWDivM = 0;
      unsigned offsets__tmp0 = 468;
      unsigned offsets_Q = Q_offset;


      device_gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 56, 9, 49, 1.0, V3mTo2nTWDivM(3,3), 56, _tmp0, 52, 1.0, Q, 56, offsets_V3mTo2nTWDivM, offsets__tmp0, offsets_Q, num_elements);
    }
    tmp_manager.free();
  }
}
#endif