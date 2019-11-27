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
    unsigned *offsets_I = I_indices;
    unsigned *offsets_star = star_indices(0);
    unsigned offsets__tmp0 = 324;
    
    
    device_gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 36, 9, 9, 1.0, I, 56, star(0), 9, 0.0, _tmp0, 36, offsets_I, offsets_star, offsets__tmp0, num_elements);
    }
    {
    unsigned offsets_kDivM = 0;
    unsigned offsets__tmp0 = 324;
    unsigned *offsets_Q = Q_indices;
    
    
    device_gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 56, 9, 35, 1.0, kDivM(0), 56, _tmp0, 36, 1.0, Q, 56, offsets_kDivM, offsets__tmp0, offsets_Q, num_elements);
    }
    _tmp2 = d_buffer0;
    {
    unsigned *offsets_I = I_indices;
    unsigned *offsets_star = star_indices(1);
    unsigned offsets__tmp2 = 324;
    
    
    device_gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 36, 9, 9, 1.0, I, 56, star(1), 9, 0.0, _tmp2, 36, offsets_I, offsets_star, offsets__tmp2, num_elements);
    }
    {
    unsigned offsets_kDivM = 0;
    unsigned offsets__tmp2 = 324;
    unsigned *offsets_Q = Q_indices;
    
    
    device_gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 56, 9, 35, 1.0, kDivM(1), 56, _tmp2, 36, 1.0, Q, 56, offsets_kDivM, offsets__tmp2, offsets_Q, num_elements);
    }
    _tmp4 = d_buffer0;
    {
    unsigned *offsets_I = I_indices;
    unsigned *offsets_star = star_indices(2);
    unsigned offsets__tmp4 = 324;
    
    
    device_gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 36, 9, 9, 1.0, I, 56, star(2), 9, 0.0, _tmp4, 36, offsets_I, offsets_star, offsets__tmp4, num_elements);
    }
    {
    unsigned offsets_kDivM = 0;
    unsigned offsets__tmp4 = 324;
    unsigned *offsets_Q = Q_indices;
    
    
    device_gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 56, 9, 35, 1.0, kDivM(2), 56, _tmp4, 36, 1.0, Q, 56, offsets_kDivM, offsets__tmp4, offsets_Q, num_elements);
    }
    tmp_manager.free();
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
    unsigned *offsets_I = I_indices;
    unsigned offsets__tmp0 = 216;
    
    
    device_gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 24, 9, 56, 1.0, fMrT(0), 24, I, 56, 0.0, _tmp0, 24, offsets_fMrT, offsets_I, offsets__tmp0, num_elements);
    }
    _tmp1 = d_buffer1;
    {
    unsigned offsets__tmp0 = 216;
    unsigned *offsets_AplusT = AplusT_indices;
    unsigned offsets__tmp1 = 216;
    
    
    device_gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 24, 9, 9, 1.0, _tmp0, 24, AplusT, 9, 0.0, _tmp1, 24, offsets__tmp0, offsets_AplusT, offsets__tmp1, num_elements);
    }
    {
    unsigned offsets_rDivM = 0;
    unsigned offsets__tmp1 = 216;
    unsigned *offsets_Q = Q_indices;
    
    
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
    unsigned *offsets_I = I_indices;
    unsigned offsets__tmp0 = 216;
    
    
    device_gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 24, 9, 56, 1.0, fMrT(1), 24, I, 56, 0.0, _tmp0, 24, offsets_fMrT, offsets_I, offsets__tmp0, num_elements);
    }
    _tmp1 = d_buffer1;
    {
    unsigned offsets__tmp0 = 216;
    unsigned *offsets_AplusT = AplusT_indices;
    unsigned offsets__tmp1 = 216;
    
    
    device_gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 24, 9, 9, 1.0, _tmp0, 24, AplusT, 9, 0.0, _tmp1, 24, offsets__tmp0, offsets_AplusT, offsets__tmp1, num_elements);
    }
    {
    unsigned offsets_rDivM = 0;
    unsigned offsets__tmp1 = 216;
    unsigned *offsets_Q = Q_indices;
    
    
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
    unsigned *offsets_I = I_indices;
    unsigned offsets__tmp0 = 216;
    
    
    device_gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 24, 9, 56, 1.0, fMrT(2), 24, I, 56, 0.0, _tmp0, 24, offsets_fMrT, offsets_I, offsets__tmp0, num_elements);
    }
    _tmp1 = d_buffer1;
    {
    unsigned offsets__tmp0 = 216;
    unsigned *offsets_AplusT = AplusT_indices;
    unsigned offsets__tmp1 = 216;
    
    
    device_gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 24, 9, 9, 1.0, _tmp0, 24, AplusT, 9, 0.0, _tmp1, 24, offsets__tmp0, offsets_AplusT, offsets__tmp1, num_elements);
    }
    {
    unsigned offsets_rDivM = 0;
    unsigned offsets__tmp1 = 216;
    unsigned *offsets_Q = Q_indices;
    
    
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
    unsigned *offsets_I = I_indices;
    unsigned offsets__tmp0 = 216;
    
    
    device_gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 24, 9, 56, 1.0, fMrT(3), 24, I, 56, 0.0, _tmp0, 24, offsets_fMrT, offsets_I, offsets__tmp0, num_elements);
    }
    _tmp1 = d_buffer1;
    {
    unsigned offsets__tmp0 = 216;
    unsigned *offsets_AplusT = AplusT_indices;
    unsigned offsets__tmp1 = 216;
    
    
    device_gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 24, 9, 9, 1.0, _tmp0, 24, AplusT, 9, 0.0, _tmp1, 24, offsets__tmp0, offsets_AplusT, offsets__tmp1, num_elements);
    }
    {
    unsigned offsets_rDivM = 0;
    unsigned offsets__tmp1 = 216;
    unsigned *offsets_Q = Q_indices;
    
    
    device_gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 56, 9, 21, 1.0, rDivM(3), 56, _tmp1, 24, 1.0, Q, 56, offsets_rDivM, offsets__tmp1, offsets_Q, num_elements);
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
    device_copy_add_scale(56, 9, power, dQ(0) + 0, 56, 0.0, I + 0, 56, dQ_indices(0), I_indices, num_elements);
  }
  void kernel::derivativeTaylorExpansion::execute1() {
    assert(!std::isnan(power));
    assert(I != nullptr);
    assert(dQ(1) != nullptr);
    assert(kernel::derivativeTaylorExpansion::num_elements != 0);
    device_copy_add_scale(35, 9, power, dQ(1) + 0, 36, 1.0, I + 0, 56, dQ_indices(1), I_indices, num_elements);
  }
  void kernel::derivativeTaylorExpansion::execute2() {
    assert(!std::isnan(power));
    assert(I != nullptr);
    assert(dQ(2) != nullptr);
    assert(kernel::derivativeTaylorExpansion::num_elements != 0);
    device_copy_add_scale(20, 9, power, dQ(2) + 0, 20, 1.0, I + 0, 56, dQ_indices(2), I_indices, num_elements);
  }
  void kernel::derivativeTaylorExpansion::execute3() {
    assert(!std::isnan(power));
    assert(I != nullptr);
    assert(dQ(3) != nullptr);
    assert(kernel::derivativeTaylorExpansion::num_elements != 0);
    device_copy_add_scale(10, 9, power, dQ(3) + 0, 12, 1.0, I + 0, 56, dQ_indices(3), I_indices, num_elements);
  }
  void kernel::derivativeTaylorExpansion::execute4() {
    assert(!std::isnan(power));
    assert(I != nullptr);
    assert(dQ(4) != nullptr);
    assert(kernel::derivativeTaylorExpansion::num_elements != 0);
    device_copy_add_scale(4, 9, power, dQ(4) + 0, 4, 1.0, I + 0, 56, dQ_indices(4), I_indices, num_elements);
  }
  void kernel::derivativeTaylorExpansion::execute5() {
    assert(!std::isnan(power));
    assert(I != nullptr);
    assert(dQ(5) != nullptr);
    assert(kernel::derivativeTaylorExpansion::num_elements != 0);
    device_copy_add_scale(1, 9, power, dQ(5) + 0, 4, 1.0, I + 0, 56, dQ_indices(5), I_indices, num_elements);
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
    unsigned *offsets_dQ = dQ_indices(0);
    unsigned offsets__tmp0 = 324;
    
    
    device_gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 36, 9, 53, 1.0, kDivMT(0), 36, dQ(0) + 1, 56, 0.0, _tmp0, 36, offsets_kDivMT, offsets_dQ, offsets__tmp0, num_elements);
    }
    {
    unsigned offsets__tmp0 = 324;
    unsigned *offsets_star = star_indices(0);
    unsigned *offsets_dQ = dQ_indices(1);
    
    
    device_gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 36, 9, 9, 1.0, _tmp0, 36, star(0), 9, 0.0, dQ(1), 36, offsets__tmp0, offsets_star, offsets_dQ, num_elements);
    }
    _tmp2 = d_buffer0;
    {
    unsigned offsets_kDivMT = 0;
    unsigned *offsets_dQ = dQ_indices(0);
    unsigned offsets__tmp2 = 324;
    
    
    device_gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 36, 9, 54, 1.0, kDivMT(1), 36, dQ(0) + 1, 56, 0.0, _tmp2, 36, offsets_kDivMT, offsets_dQ, offsets__tmp2, num_elements);
    }
    {
    unsigned offsets__tmp2 = 324;
    unsigned *offsets_star = star_indices(1);
    unsigned *offsets_dQ = dQ_indices(1);
    
    
    device_gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 36, 9, 9, 1.0, _tmp2, 36, star(1), 9, 1.0, dQ(1), 36, offsets__tmp2, offsets_star, offsets_dQ, num_elements);
    }
    _tmp4 = d_buffer0;
    {
    unsigned offsets_kDivMT = 0;
    unsigned *offsets_dQ = dQ_indices(0);
    unsigned offsets__tmp4 = 324;
    
    
    device_gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 36, 9, 55, 1.0, kDivMT(2), 36, dQ(0) + 1, 56, 0.0, _tmp4, 36, offsets_kDivMT, offsets_dQ, offsets__tmp4, num_elements);
    }
    {
    unsigned offsets__tmp4 = 324;
    unsigned *offsets_star = star_indices(2);
    unsigned *offsets_dQ = dQ_indices(1);
    
    
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
    unsigned *offsets_dQ = dQ_indices(1);
    unsigned offsets__tmp0 = 180;
    
    
    device_gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 20, 9, 32, 1.0, kDivMT(0), 36, dQ(1) + 1, 36, 0.0, _tmp0, 20, offsets_kDivMT, offsets_dQ, offsets__tmp0, num_elements);
    }
    {
    unsigned offsets__tmp0 = 180;
    unsigned *offsets_star = star_indices(0);
    unsigned *offsets_dQ = dQ_indices(2);
    
    
    device_gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 20, 9, 9, 1.0, _tmp0, 20, star(0), 9, 0.0, dQ(2), 20, offsets__tmp0, offsets_star, offsets_dQ, num_elements);
    }
    _tmp2 = d_buffer0;
    {
    unsigned offsets_kDivMT = 0;
    unsigned *offsets_dQ = dQ_indices(1);
    unsigned offsets__tmp2 = 180;
    
    
    device_gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 20, 9, 33, 1.0, kDivMT(1), 36, dQ(1) + 1, 36, 0.0, _tmp2, 20, offsets_kDivMT, offsets_dQ, offsets__tmp2, num_elements);
    }
    {
    unsigned offsets__tmp2 = 180;
    unsigned *offsets_star = star_indices(1);
    unsigned *offsets_dQ = dQ_indices(2);
    
    
    device_gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 20, 9, 9, 1.0, _tmp2, 20, star(1), 9, 1.0, dQ(2), 20, offsets__tmp2, offsets_star, offsets_dQ, num_elements);
    }
    _tmp4 = d_buffer0;
    {
    unsigned offsets_kDivMT = 0;
    unsigned *offsets_dQ = dQ_indices(1);
    unsigned offsets__tmp4 = 180;
    
    
    device_gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 20, 9, 34, 1.0, kDivMT(2), 36, dQ(1) + 1, 36, 0.0, _tmp4, 20, offsets_kDivMT, offsets_dQ, offsets__tmp4, num_elements);
    }
    {
    unsigned offsets__tmp4 = 180;
    unsigned *offsets_star = star_indices(2);
    unsigned *offsets_dQ = dQ_indices(2);
    
    
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
    unsigned *offsets_dQ = dQ_indices(2);
    unsigned offsets__tmp0 = 108;
    
    
    device_gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 12, 9, 17, 1.0, kDivMT(0), 36, dQ(2) + 1, 20, 0.0, _tmp0, 12, offsets_kDivMT, offsets_dQ, offsets__tmp0, num_elements);
    }
    {
    unsigned offsets__tmp0 = 108;
    unsigned *offsets_star = star_indices(0);
    unsigned *offsets_dQ = dQ_indices(3);
    
    
    device_gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 12, 9, 9, 1.0, _tmp0, 12, star(0), 9, 0.0, dQ(3), 12, offsets__tmp0, offsets_star, offsets_dQ, num_elements);
    }
    _tmp2 = d_buffer0;
    {
    unsigned offsets_kDivMT = 0;
    unsigned *offsets_dQ = dQ_indices(2);
    unsigned offsets__tmp2 = 108;
    
    
    device_gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 12, 9, 18, 1.0, kDivMT(1), 36, dQ(2) + 1, 20, 0.0, _tmp2, 12, offsets_kDivMT, offsets_dQ, offsets__tmp2, num_elements);
    }
    {
    unsigned offsets__tmp2 = 108;
    unsigned *offsets_star = star_indices(1);
    unsigned *offsets_dQ = dQ_indices(3);
    
    
    device_gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 12, 9, 9, 1.0, _tmp2, 12, star(1), 9, 1.0, dQ(3), 12, offsets__tmp2, offsets_star, offsets_dQ, num_elements);
    }
    _tmp4 = d_buffer0;
    {
    unsigned offsets_kDivMT = 0;
    unsigned *offsets_dQ = dQ_indices(2);
    unsigned offsets__tmp4 = 108;
    
    
    device_gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 12, 9, 19, 1.0, kDivMT(2), 36, dQ(2) + 1, 20, 0.0, _tmp4, 12, offsets_kDivMT, offsets_dQ, offsets__tmp4, num_elements);
    }
    {
    unsigned offsets__tmp4 = 108;
    unsigned *offsets_star = star_indices(2);
    unsigned *offsets_dQ = dQ_indices(3);
    
    
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
    unsigned *offsets_dQ = dQ_indices(3);
    unsigned offsets__tmp0 = 36;
    
    
    device_gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 4, 9, 7, 1.0, kDivMT(0), 36, dQ(3) + 1, 12, 0.0, _tmp0, 4, offsets_kDivMT, offsets_dQ, offsets__tmp0, num_elements);
    }
    {
    unsigned offsets__tmp0 = 36;
    unsigned *offsets_star = star_indices(0);
    unsigned *offsets_dQ = dQ_indices(4);
    
    
    device_gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 4, 9, 9, 1.0, _tmp0, 4, star(0), 9, 0.0, dQ(4), 4, offsets__tmp0, offsets_star, offsets_dQ, num_elements);
    }
    _tmp2 = d_buffer0;
    {
    unsigned offsets_kDivMT = 0;
    unsigned *offsets_dQ = dQ_indices(3);
    unsigned offsets__tmp2 = 36;
    
    
    device_gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 4, 9, 8, 1.0, kDivMT(1), 36, dQ(3) + 1, 12, 0.0, _tmp2, 4, offsets_kDivMT, offsets_dQ, offsets__tmp2, num_elements);
    }
    {
    unsigned offsets__tmp2 = 36;
    unsigned *offsets_star = star_indices(1);
    unsigned *offsets_dQ = dQ_indices(4);
    
    
    device_gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 4, 9, 9, 1.0, _tmp2, 4, star(1), 9, 1.0, dQ(4), 4, offsets__tmp2, offsets_star, offsets_dQ, num_elements);
    }
    _tmp4 = d_buffer0;
    {
    unsigned offsets_kDivMT = 0;
    unsigned *offsets_dQ = dQ_indices(3);
    unsigned offsets__tmp4 = 36;
    
    
    device_gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 4, 9, 9, 1.0, kDivMT(2), 36, dQ(3) + 1, 12, 0.0, _tmp4, 4, offsets_kDivMT, offsets_dQ, offsets__tmp4, num_elements);
    }
    {
    unsigned offsets__tmp4 = 36;
    unsigned *offsets_star = star_indices(2);
    unsigned *offsets_dQ = dQ_indices(4);
    
    
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
    unsigned *offsets_dQ = dQ_indices(4);
    unsigned offsets__tmp0 = 36;
    
    
    device_gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 4, 9, 1, 1.0, kDivMT(0), 36, dQ(4) + 1, 4, 0.0, _tmp0, 4, offsets_kDivMT, offsets_dQ, offsets__tmp0, num_elements);
    }
    {
    unsigned offsets__tmp0 = 36;
    unsigned *offsets_star = star_indices(0);
    unsigned *offsets_dQ = dQ_indices(5);
    
    
    device_gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 4, 9, 9, 1.0, _tmp0, 4, star(0), 9, 0.0, dQ(5), 4, offsets__tmp0, offsets_star, offsets_dQ, num_elements);
    }
    _tmp2 = d_buffer0;
    {
    unsigned offsets_kDivMT = 0;
    unsigned *offsets_dQ = dQ_indices(4);
    unsigned offsets__tmp2 = 36;
    
    
    device_gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 4, 9, 2, 1.0, kDivMT(1), 36, dQ(4) + 1, 4, 0.0, _tmp2, 4, offsets_kDivMT, offsets_dQ, offsets__tmp2, num_elements);
    }
    {
    unsigned offsets__tmp2 = 36;
    unsigned *offsets_star = star_indices(1);
    unsigned *offsets_dQ = dQ_indices(5);
    
    
    device_gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 4, 9, 9, 1.0, _tmp2, 4, star(1), 9, 1.0, dQ(5), 4, offsets__tmp2, offsets_star, offsets_dQ, num_elements);
    }
    _tmp4 = d_buffer0;
    {
    unsigned offsets_kDivMT = 0;
    unsigned *offsets_dQ = dQ_indices(4);
    unsigned offsets__tmp4 = 36;
    
    
    device_gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 4, 9, 3, 1.0, kDivMT(2), 36, dQ(4) + 1, 4, 0.0, _tmp4, 4, offsets_kDivMT, offsets_dQ, offsets__tmp4, num_elements);
    }
    {
    unsigned offsets__tmp4 = 36;
    unsigned *offsets_star = star_indices(2);
    unsigned *offsets_dQ = dQ_indices(5);
    
    
    device_gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 4, 9, 9, 1.0, _tmp4, 4, star(2), 9, 1.0, dQ(5), 4, offsets__tmp4, offsets_star, offsets_dQ, num_elements);
    }
    tmp_manager.free();
  }
}
#endif