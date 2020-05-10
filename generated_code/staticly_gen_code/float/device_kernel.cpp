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
    assert(NumElements != 0);
    device::DeviceInstance& device = device::DeviceInstance::getInstance();
    float *_tmp0, *_tmp2, *_tmp4;
    float* _buffer0 = (float*)device.api->getStackMemory(sizeof(float) * NumElements * 432);
    _tmp0 = _buffer0;

    sgemm_NT_NT_m35_64_n9_9_k9_pps_91656fc(const_cast<const float **>(I), ExtraOffset_I, const_cast<const float **>(star(0)), ExtraOffset_star(0), _tmp0, 0, NumElements);

    sgemm_NT_NT_m53_64_n9_48_k35_nsp_49df72e(const_cast<const float *>(kDivM(0)), 0, const_cast<const float *>(_tmp0), 0, Q, ExtraOffset_Q, NumElements);
    _tmp2 = _buffer0;

    sgemm_NT_NT_m35_64_n9_9_k9_pps_91656fc(const_cast<const float **>(I), ExtraOffset_I, const_cast<const float **>(star(1)), ExtraOffset_star(1), _tmp2, 0, NumElements);

    sgemm_NT_NT_m54_64_n9_48_k35_nsp_13f0270(const_cast<const float *>(kDivM(1)), 0, const_cast<const float *>(_tmp2), 0, Q, ExtraOffset_Q, NumElements);
    _tmp4 = _buffer0;

    sgemm_NT_NT_m35_64_n9_9_k9_pps_91656fc(const_cast<const float **>(I), ExtraOffset_I, const_cast<const float **>(star(2)), ExtraOffset_star(2), _tmp4, 0, NumElements);

    sgemm_NT_NT_m55_64_n9_48_k35_nsp_6daa8d7(const_cast<const float *>(kDivM(2)), 0, const_cast<const float *>(_tmp4), 0, Q, ExtraOffset_Q, NumElements);
    device.api->popStackMemory();
  }
  constexpr unsigned long const kernel::rotateGodunovStateLocal::NonZeroFlops;
  constexpr unsigned long const kernel::rotateGodunovStateLocal::HardwareFlops;
  void kernel::rotateGodunovStateLocal::execute() {
    assert(QgodLocal != nullptr);
    assert(Tinv != nullptr);
    assert(godunovMatrix != nullptr);
    assert(NumElements != 0);
    device::DeviceInstance& device = device::DeviceInstance::getInstance();

    sgemm_T_NT_m9_9_n9_9_k9_ppp_db79476(const_cast<const float **>(Tinv), ExtraOffset_Tinv, const_cast<const float **>(QgodLocal), ExtraOffset_QgodLocal, godunovMatrix, ExtraOffset_godunovMatrix, NumElements);
  }
  constexpr unsigned long const kernel::rotateGodunovStateNeighbor::NonZeroFlops;
  constexpr unsigned long const kernel::rotateGodunovStateNeighbor::HardwareFlops;
  void kernel::rotateGodunovStateNeighbor::execute() {
    assert(QgodNeighbor != nullptr);
    assert(Tinv != nullptr);
    assert(godunovMatrix != nullptr);
    assert(NumElements != 0);
    device::DeviceInstance& device = device::DeviceInstance::getInstance();

    sgemm_T_NT_m9_9_n9_9_k9_ppp_db79476(const_cast<const float **>(Tinv), ExtraOffset_Tinv, const_cast<const float **>(QgodNeighbor), ExtraOffset_QgodNeighbor, godunovMatrix, ExtraOffset_godunovMatrix, NumElements);
  }
  constexpr unsigned long const kernel::rotateFluxMatrix::NonZeroFlops;
  constexpr unsigned long const kernel::rotateFluxMatrix::HardwareFlops;
  void kernel::rotateFluxMatrix::execute() {
    assert(!std::isnan(fluxScale));
    assert(T != nullptr);
    assert(fluxSolver != nullptr);
    assert(star(0) != nullptr);
    assert(NumElements != 0);
    device::DeviceInstance& device = device::DeviceInstance::getInstance();

    sgemm_NT_T_m9_9_n9_9_k9_ppp_0c970c0(fluxScale, const_cast<const float **>(star(0)), ExtraOffset_star(0), const_cast<const float **>(T), ExtraOffset_T, fluxSolver, ExtraOffset_fluxSolver, NumElements);
  }
  constexpr unsigned long const kernel::plConvertToNodal::NonZeroFlops;
  constexpr unsigned long const kernel::plConvertToNodal::HardwareFlops;
  void kernel::plConvertToNodal::execute() {
    assert(QStress != nullptr);
    assert(QStressNodal != nullptr);
    assert(v != nullptr);
    assert(NumElements != 0);
    device::DeviceInstance& device = device::DeviceInstance::getInstance();

    sgemm_NT_NT_m56_64_n6_64_k56_npp_968d18f(const_cast<const float *>(v), 0, const_cast<const float **>(QStress), ExtraOffset_QStress, QStressNodal, ExtraOffset_QStressNodal, NumElements);
  }
  constexpr unsigned long const kernel::plConvertToModal::NonZeroFlops;
  constexpr unsigned long const kernel::plConvertToModal::HardwareFlops;
  void kernel::plConvertToModal::execute() {
    assert(QStress != nullptr);
    assert(QStressNodal != nullptr);
    assert(vInv != nullptr);
    assert(NumElements != 0);
    device::DeviceInstance& device = device::DeviceInstance::getInstance();

    sgemm_NT_NT_m56_64_n6_64_k56_npp_3250429(const_cast<const float *>(vInv), 0, const_cast<const float **>(QStressNodal), ExtraOffset_QStressNodal, QStress, ExtraOffset_QStress, NumElements);
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
    assert(NumElements != 0);
    device::DeviceInstance& device = device::DeviceInstance::getInstance();
    float *_tmp0, *_tmp1;
    float* _buffer0 = (float*)device.api->getStackMemory(sizeof(float) * NumElements * 288);
    float* _buffer1 = (float*)device.api->getStackMemory(sizeof(float) * NumElements * 288);
    _tmp0 = _buffer0;

    sgemm_NT_NT_m21_32_n9_64_k56_nps_3c76782(const_cast<const float *>(fMrT(0)), 0, const_cast<const float **>(I), ExtraOffset_I, _tmp0, 0, NumElements);
    _tmp1 = _buffer1;

    sgemm_NT_NT_m21_32_n9_9_k9_sps_97c1866(const_cast<const float *>(_tmp0), 0, const_cast<const float **>(AplusT), ExtraOffset_AplusT, _tmp1, 0, NumElements);

    sgemm_NT_NT_m56_64_n9_32_k21_nsp_0e604b1(const_cast<const float *>(rDivM(0)), 0, const_cast<const float *>(_tmp1), 0, Q, ExtraOffset_Q, NumElements);
    device.api->popStackMemory();
    device.api->popStackMemory();
  }
  void kernel::localFlux::execute1() {
    assert(AplusT != nullptr);
    assert(I != nullptr);
    assert(Q != nullptr);
    assert(fMrT(1) != nullptr);
    assert(rDivM(1) != nullptr);
    assert(NumElements != 0);
    device::DeviceInstance& device = device::DeviceInstance::getInstance();
    float *_tmp0, *_tmp1;
    float* _buffer0 = (float*)device.api->getStackMemory(sizeof(float) * NumElements * 288);
    float* _buffer1 = (float*)device.api->getStackMemory(sizeof(float) * NumElements * 288);
    _tmp0 = _buffer0;

    sgemm_NT_NT_m21_32_n9_64_k56_nps_3c76782(const_cast<const float *>(fMrT(1)), 0, const_cast<const float **>(I), ExtraOffset_I, _tmp0, 0, NumElements);
    _tmp1 = _buffer1;

    sgemm_NT_NT_m21_32_n9_9_k9_sps_97c1866(const_cast<const float *>(_tmp0), 0, const_cast<const float **>(AplusT), ExtraOffset_AplusT, _tmp1, 0, NumElements);

    sgemm_NT_NT_m56_64_n9_32_k21_nsp_0e604b1(const_cast<const float *>(rDivM(1)), 0, const_cast<const float *>(_tmp1), 0, Q, ExtraOffset_Q, NumElements);
    device.api->popStackMemory();
    device.api->popStackMemory();
  }
  void kernel::localFlux::execute2() {
    assert(AplusT != nullptr);
    assert(I != nullptr);
    assert(Q != nullptr);
    assert(fMrT(2) != nullptr);
    assert(rDivM(2) != nullptr);
    assert(NumElements != 0);
    device::DeviceInstance& device = device::DeviceInstance::getInstance();
    float *_tmp0, *_tmp1;
    float* _buffer0 = (float*)device.api->getStackMemory(sizeof(float) * NumElements * 288);
    float* _buffer1 = (float*)device.api->getStackMemory(sizeof(float) * NumElements * 288);
    _tmp0 = _buffer0;

    sgemm_NT_NT_m21_32_n9_64_k56_nps_3c76782(const_cast<const float *>(fMrT(2)), 0, const_cast<const float **>(I), ExtraOffset_I, _tmp0, 0, NumElements);
    _tmp1 = _buffer1;

    sgemm_NT_NT_m21_32_n9_9_k9_sps_97c1866(const_cast<const float *>(_tmp0), 0, const_cast<const float **>(AplusT), ExtraOffset_AplusT, _tmp1, 0, NumElements);

    sgemm_NT_NT_m56_64_n9_32_k21_nsp_0e604b1(const_cast<const float *>(rDivM(2)), 0, const_cast<const float *>(_tmp1), 0, Q, ExtraOffset_Q, NumElements);
    device.api->popStackMemory();
    device.api->popStackMemory();
  }
  void kernel::localFlux::execute3() {
    assert(AplusT != nullptr);
    assert(I != nullptr);
    assert(Q != nullptr);
    assert(fMrT(3) != nullptr);
    assert(rDivM(3) != nullptr);
    assert(NumElements != 0);
    device::DeviceInstance& device = device::DeviceInstance::getInstance();
    float *_tmp0, *_tmp1;
    float* _buffer0 = (float*)device.api->getStackMemory(sizeof(float) * NumElements * 288);
    float* _buffer1 = (float*)device.api->getStackMemory(sizeof(float) * NumElements * 288);
    _tmp0 = _buffer0;

    sgemm_NT_NT_m21_32_n9_64_k56_nps_3c76782(const_cast<const float *>(fMrT(3)), 0, const_cast<const float **>(I), ExtraOffset_I, _tmp0, 0, NumElements);
    _tmp1 = _buffer1;

    sgemm_NT_NT_m21_32_n9_9_k9_sps_97c1866(const_cast<const float *>(_tmp0), 0, const_cast<const float **>(AplusT), ExtraOffset_AplusT, _tmp1, 0, NumElements);

    sgemm_NT_NT_m56_64_n9_32_k21_nsp_0e604b1(const_cast<const float *>(rDivM(3)), 0, const_cast<const float *>(_tmp1), 0, Q, ExtraOffset_Q, NumElements);
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
    assert(NumElements != 0);
    device::DeviceInstance& device = device::DeviceInstance::getInstance();
    float *_tmp0, *_tmp1, *_tmp2;
    float* _buffer0 = (float*)device.api->getStackMemory(sizeof(float) * NumElements * 288);
    float* _buffer1 = (float*)device.api->getStackMemory(sizeof(float) * NumElements * 288);
    _tmp0 = _buffer0;

    sgemm_NT_NT_m21_32_n9_64_k56_nps_3c76782(const_cast<const float *>(rT(0)), 0, const_cast<const float **>(I), ExtraOffset_I, _tmp0, 0, NumElements);
    _tmp1 = _buffer1;

    sgemm_NT_NT_m21_32_n9_32_k21_nss_3b271a8(const_cast<const float *>(fP(0)), 0, const_cast<const float *>(_tmp0), 0, _tmp1, 0, NumElements);
    _tmp2 = _buffer0;

    sgemm_NT_NT_m21_32_n9_9_k9_sps_97c1866(const_cast<const float *>(_tmp1), 0, const_cast<const float **>(AminusT), ExtraOffset_AminusT, _tmp2, 0, NumElements);

    sgemm_NT_NT_m56_64_n9_32_k21_nsp_0e604b1(const_cast<const float *>(rDivM(0)), 0, const_cast<const float *>(_tmp2), 0, Q, ExtraOffset_Q, NumElements);
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
    assert(NumElements != 0);
    device::DeviceInstance& device = device::DeviceInstance::getInstance();
    float *_tmp0, *_tmp1, *_tmp2;
    float* _buffer0 = (float*)device.api->getStackMemory(sizeof(float) * NumElements * 288);
    float* _buffer1 = (float*)device.api->getStackMemory(sizeof(float) * NumElements * 288);
    _tmp0 = _buffer0;

    sgemm_NT_NT_m21_32_n9_64_k56_nps_3c76782(const_cast<const float *>(rT(0)), 0, const_cast<const float **>(I), ExtraOffset_I, _tmp0, 0, NumElements);
    _tmp1 = _buffer1;

    sgemm_NT_NT_m21_32_n9_32_k21_nss_3b271a8(const_cast<const float *>(fP(1)), 0, const_cast<const float *>(_tmp0), 0, _tmp1, 0, NumElements);
    _tmp2 = _buffer0;

    sgemm_NT_NT_m21_32_n9_9_k9_sps_97c1866(const_cast<const float *>(_tmp1), 0, const_cast<const float **>(AminusT), ExtraOffset_AminusT, _tmp2, 0, NumElements);

    sgemm_NT_NT_m56_64_n9_32_k21_nsp_0e604b1(const_cast<const float *>(rDivM(0)), 0, const_cast<const float *>(_tmp2), 0, Q, ExtraOffset_Q, NumElements);
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
    assert(NumElements != 0);
    device::DeviceInstance& device = device::DeviceInstance::getInstance();
    float *_tmp0, *_tmp1, *_tmp2;
    float* _buffer0 = (float*)device.api->getStackMemory(sizeof(float) * NumElements * 288);
    float* _buffer1 = (float*)device.api->getStackMemory(sizeof(float) * NumElements * 288);
    _tmp0 = _buffer0;

    sgemm_NT_NT_m21_32_n9_64_k56_nps_3c76782(const_cast<const float *>(rT(0)), 0, const_cast<const float **>(I), ExtraOffset_I, _tmp0, 0, NumElements);
    _tmp1 = _buffer1;

    sgemm_NT_NT_m21_32_n9_32_k21_nss_3b271a8(const_cast<const float *>(fP(2)), 0, const_cast<const float *>(_tmp0), 0, _tmp1, 0, NumElements);
    _tmp2 = _buffer0;

    sgemm_NT_NT_m21_32_n9_9_k9_sps_97c1866(const_cast<const float *>(_tmp1), 0, const_cast<const float **>(AminusT), ExtraOffset_AminusT, _tmp2, 0, NumElements);

    sgemm_NT_NT_m56_64_n9_32_k21_nsp_0e604b1(const_cast<const float *>(rDivM(0)), 0, const_cast<const float *>(_tmp2), 0, Q, ExtraOffset_Q, NumElements);
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
    assert(NumElements != 0);
    device::DeviceInstance& device = device::DeviceInstance::getInstance();
    float *_tmp0, *_tmp1, *_tmp2;
    float* _buffer0 = (float*)device.api->getStackMemory(sizeof(float) * NumElements * 288);
    float* _buffer1 = (float*)device.api->getStackMemory(sizeof(float) * NumElements * 288);
    _tmp0 = _buffer0;

    sgemm_NT_NT_m21_32_n9_64_k56_nps_3c76782(const_cast<const float *>(rT(1)), 0, const_cast<const float **>(I), ExtraOffset_I, _tmp0, 0, NumElements);
    _tmp1 = _buffer1;

    sgemm_NT_NT_m21_32_n9_32_k21_nss_3b271a8(const_cast<const float *>(fP(0)), 0, const_cast<const float *>(_tmp0), 0, _tmp1, 0, NumElements);
    _tmp2 = _buffer0;

    sgemm_NT_NT_m21_32_n9_9_k9_sps_97c1866(const_cast<const float *>(_tmp1), 0, const_cast<const float **>(AminusT), ExtraOffset_AminusT, _tmp2, 0, NumElements);

    sgemm_NT_NT_m56_64_n9_32_k21_nsp_0e604b1(const_cast<const float *>(rDivM(0)), 0, const_cast<const float *>(_tmp2), 0, Q, ExtraOffset_Q, NumElements);
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
    assert(NumElements != 0);
    device::DeviceInstance& device = device::DeviceInstance::getInstance();
    float *_tmp0, *_tmp1, *_tmp2;
    float* _buffer0 = (float*)device.api->getStackMemory(sizeof(float) * NumElements * 288);
    float* _buffer1 = (float*)device.api->getStackMemory(sizeof(float) * NumElements * 288);
    _tmp0 = _buffer0;

    sgemm_NT_NT_m21_32_n9_64_k56_nps_3c76782(const_cast<const float *>(rT(1)), 0, const_cast<const float **>(I), ExtraOffset_I, _tmp0, 0, NumElements);
    _tmp1 = _buffer1;

    sgemm_NT_NT_m21_32_n9_32_k21_nss_3b271a8(const_cast<const float *>(fP(1)), 0, const_cast<const float *>(_tmp0), 0, _tmp1, 0, NumElements);
    _tmp2 = _buffer0;

    sgemm_NT_NT_m21_32_n9_9_k9_sps_97c1866(const_cast<const float *>(_tmp1), 0, const_cast<const float **>(AminusT), ExtraOffset_AminusT, _tmp2, 0, NumElements);

    sgemm_NT_NT_m56_64_n9_32_k21_nsp_0e604b1(const_cast<const float *>(rDivM(0)), 0, const_cast<const float *>(_tmp2), 0, Q, ExtraOffset_Q, NumElements);
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
    assert(NumElements != 0);
    device::DeviceInstance& device = device::DeviceInstance::getInstance();
    float *_tmp0, *_tmp1, *_tmp2;
    float* _buffer0 = (float*)device.api->getStackMemory(sizeof(float) * NumElements * 288);
    float* _buffer1 = (float*)device.api->getStackMemory(sizeof(float) * NumElements * 288);
    _tmp0 = _buffer0;

    sgemm_NT_NT_m21_32_n9_64_k56_nps_3c76782(const_cast<const float *>(rT(1)), 0, const_cast<const float **>(I), ExtraOffset_I, _tmp0, 0, NumElements);
    _tmp1 = _buffer1;

    sgemm_NT_NT_m21_32_n9_32_k21_nss_3b271a8(const_cast<const float *>(fP(2)), 0, const_cast<const float *>(_tmp0), 0, _tmp1, 0, NumElements);
    _tmp2 = _buffer0;

    sgemm_NT_NT_m21_32_n9_9_k9_sps_97c1866(const_cast<const float *>(_tmp1), 0, const_cast<const float **>(AminusT), ExtraOffset_AminusT, _tmp2, 0, NumElements);

    sgemm_NT_NT_m56_64_n9_32_k21_nsp_0e604b1(const_cast<const float *>(rDivM(0)), 0, const_cast<const float *>(_tmp2), 0, Q, ExtraOffset_Q, NumElements);
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
    assert(NumElements != 0);
    device::DeviceInstance& device = device::DeviceInstance::getInstance();
    float *_tmp0, *_tmp1, *_tmp2;
    float* _buffer0 = (float*)device.api->getStackMemory(sizeof(float) * NumElements * 288);
    float* _buffer1 = (float*)device.api->getStackMemory(sizeof(float) * NumElements * 288);
    _tmp0 = _buffer0;

    sgemm_NT_NT_m21_32_n9_64_k56_nps_3c76782(const_cast<const float *>(rT(2)), 0, const_cast<const float **>(I), ExtraOffset_I, _tmp0, 0, NumElements);
    _tmp1 = _buffer1;

    sgemm_NT_NT_m21_32_n9_32_k21_nss_3b271a8(const_cast<const float *>(fP(0)), 0, const_cast<const float *>(_tmp0), 0, _tmp1, 0, NumElements);
    _tmp2 = _buffer0;

    sgemm_NT_NT_m21_32_n9_9_k9_sps_97c1866(const_cast<const float *>(_tmp1), 0, const_cast<const float **>(AminusT), ExtraOffset_AminusT, _tmp2, 0, NumElements);

    sgemm_NT_NT_m56_64_n9_32_k21_nsp_0e604b1(const_cast<const float *>(rDivM(0)), 0, const_cast<const float *>(_tmp2), 0, Q, ExtraOffset_Q, NumElements);
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
    assert(NumElements != 0);
    device::DeviceInstance& device = device::DeviceInstance::getInstance();
    float *_tmp0, *_tmp1, *_tmp2;
    float* _buffer0 = (float*)device.api->getStackMemory(sizeof(float) * NumElements * 288);
    float* _buffer1 = (float*)device.api->getStackMemory(sizeof(float) * NumElements * 288);
    _tmp0 = _buffer0;

    sgemm_NT_NT_m21_32_n9_64_k56_nps_3c76782(const_cast<const float *>(rT(2)), 0, const_cast<const float **>(I), ExtraOffset_I, _tmp0, 0, NumElements);
    _tmp1 = _buffer1;

    sgemm_NT_NT_m21_32_n9_32_k21_nss_3b271a8(const_cast<const float *>(fP(1)), 0, const_cast<const float *>(_tmp0), 0, _tmp1, 0, NumElements);
    _tmp2 = _buffer0;

    sgemm_NT_NT_m21_32_n9_9_k9_sps_97c1866(const_cast<const float *>(_tmp1), 0, const_cast<const float **>(AminusT), ExtraOffset_AminusT, _tmp2, 0, NumElements);

    sgemm_NT_NT_m56_64_n9_32_k21_nsp_0e604b1(const_cast<const float *>(rDivM(0)), 0, const_cast<const float *>(_tmp2), 0, Q, ExtraOffset_Q, NumElements);
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
    assert(NumElements != 0);
    device::DeviceInstance& device = device::DeviceInstance::getInstance();
    float *_tmp0, *_tmp1, *_tmp2;
    float* _buffer0 = (float*)device.api->getStackMemory(sizeof(float) * NumElements * 288);
    float* _buffer1 = (float*)device.api->getStackMemory(sizeof(float) * NumElements * 288);
    _tmp0 = _buffer0;

    sgemm_NT_NT_m21_32_n9_64_k56_nps_3c76782(const_cast<const float *>(rT(2)), 0, const_cast<const float **>(I), ExtraOffset_I, _tmp0, 0, NumElements);
    _tmp1 = _buffer1;

    sgemm_NT_NT_m21_32_n9_32_k21_nss_3b271a8(const_cast<const float *>(fP(2)), 0, const_cast<const float *>(_tmp0), 0, _tmp1, 0, NumElements);
    _tmp2 = _buffer0;

    sgemm_NT_NT_m21_32_n9_9_k9_sps_97c1866(const_cast<const float *>(_tmp1), 0, const_cast<const float **>(AminusT), ExtraOffset_AminusT, _tmp2, 0, NumElements);

    sgemm_NT_NT_m56_64_n9_32_k21_nsp_0e604b1(const_cast<const float *>(rDivM(0)), 0, const_cast<const float *>(_tmp2), 0, Q, ExtraOffset_Q, NumElements);
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
    assert(NumElements != 0);
    device::DeviceInstance& device = device::DeviceInstance::getInstance();
    float *_tmp0, *_tmp1, *_tmp2;
    float* _buffer0 = (float*)device.api->getStackMemory(sizeof(float) * NumElements * 288);
    float* _buffer1 = (float*)device.api->getStackMemory(sizeof(float) * NumElements * 288);
    _tmp0 = _buffer0;

    sgemm_NT_NT_m21_32_n9_64_k56_nps_3c76782(const_cast<const float *>(rT(3)), 0, const_cast<const float **>(I), ExtraOffset_I, _tmp0, 0, NumElements);
    _tmp1 = _buffer1;

    sgemm_NT_NT_m21_32_n9_32_k21_nss_3b271a8(const_cast<const float *>(fP(0)), 0, const_cast<const float *>(_tmp0), 0, _tmp1, 0, NumElements);
    _tmp2 = _buffer0;

    sgemm_NT_NT_m21_32_n9_9_k9_sps_97c1866(const_cast<const float *>(_tmp1), 0, const_cast<const float **>(AminusT), ExtraOffset_AminusT, _tmp2, 0, NumElements);

    sgemm_NT_NT_m56_64_n9_32_k21_nsp_0e604b1(const_cast<const float *>(rDivM(0)), 0, const_cast<const float *>(_tmp2), 0, Q, ExtraOffset_Q, NumElements);
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
    assert(NumElements != 0);
    device::DeviceInstance& device = device::DeviceInstance::getInstance();
    float *_tmp0, *_tmp1, *_tmp2;
    float* _buffer0 = (float*)device.api->getStackMemory(sizeof(float) * NumElements * 288);
    float* _buffer1 = (float*)device.api->getStackMemory(sizeof(float) * NumElements * 288);
    _tmp0 = _buffer0;

    sgemm_NT_NT_m21_32_n9_64_k56_nps_3c76782(const_cast<const float *>(rT(3)), 0, const_cast<const float **>(I), ExtraOffset_I, _tmp0, 0, NumElements);
    _tmp1 = _buffer1;

    sgemm_NT_NT_m21_32_n9_32_k21_nss_3b271a8(const_cast<const float *>(fP(1)), 0, const_cast<const float *>(_tmp0), 0, _tmp1, 0, NumElements);
    _tmp2 = _buffer0;

    sgemm_NT_NT_m21_32_n9_9_k9_sps_97c1866(const_cast<const float *>(_tmp1), 0, const_cast<const float **>(AminusT), ExtraOffset_AminusT, _tmp2, 0, NumElements);

    sgemm_NT_NT_m56_64_n9_32_k21_nsp_0e604b1(const_cast<const float *>(rDivM(0)), 0, const_cast<const float *>(_tmp2), 0, Q, ExtraOffset_Q, NumElements);
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
    assert(NumElements != 0);
    device::DeviceInstance& device = device::DeviceInstance::getInstance();
    float *_tmp0, *_tmp1, *_tmp2;
    float* _buffer0 = (float*)device.api->getStackMemory(sizeof(float) * NumElements * 288);
    float* _buffer1 = (float*)device.api->getStackMemory(sizeof(float) * NumElements * 288);
    _tmp0 = _buffer0;

    sgemm_NT_NT_m21_32_n9_64_k56_nps_3c76782(const_cast<const float *>(rT(3)), 0, const_cast<const float **>(I), ExtraOffset_I, _tmp0, 0, NumElements);
    _tmp1 = _buffer1;

    sgemm_NT_NT_m21_32_n9_32_k21_nss_3b271a8(const_cast<const float *>(fP(2)), 0, const_cast<const float *>(_tmp0), 0, _tmp1, 0, NumElements);
    _tmp2 = _buffer0;

    sgemm_NT_NT_m21_32_n9_9_k9_sps_97c1866(const_cast<const float *>(_tmp1), 0, const_cast<const float **>(AminusT), ExtraOffset_AminusT, _tmp2, 0, NumElements);

    sgemm_NT_NT_m56_64_n9_32_k21_nsp_0e604b1(const_cast<const float *>(rDivM(0)), 0, const_cast<const float *>(_tmp2), 0, Q, ExtraOffset_Q, NumElements);
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
    assert(NumElements != 0);
    device::DeviceInstance& device = device::DeviceInstance::getInstance();
    float *_tmp0, *_tmp1, *_tmp2;
    float* _buffer0 = (float*)device.api->getStackMemory(sizeof(float) * NumElements * 288);
    float* _buffer1 = (float*)device.api->getStackMemory(sizeof(float) * NumElements * 288);
    _tmp0 = _buffer0;

    sgemm_NT_NT_m21_32_n9_64_k56_nps_3c76782(const_cast<const float *>(rT(0)), 0, const_cast<const float **>(I), ExtraOffset_I, _tmp0, 0, NumElements);
    _tmp1 = _buffer1;

    sgemm_NT_NT_m21_32_n9_32_k21_nss_3b271a8(const_cast<const float *>(fP(0)), 0, const_cast<const float *>(_tmp0), 0, _tmp1, 0, NumElements);
    _tmp2 = _buffer0;

    sgemm_NT_NT_m21_32_n9_9_k9_sps_97c1866(const_cast<const float *>(_tmp1), 0, const_cast<const float **>(AminusT), ExtraOffset_AminusT, _tmp2, 0, NumElements);

    sgemm_NT_NT_m56_64_n9_32_k21_nsp_0e604b1(const_cast<const float *>(rDivM(1)), 0, const_cast<const float *>(_tmp2), 0, Q, ExtraOffset_Q, NumElements);
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
    assert(NumElements != 0);
    device::DeviceInstance& device = device::DeviceInstance::getInstance();
    float *_tmp0, *_tmp1, *_tmp2;
    float* _buffer0 = (float*)device.api->getStackMemory(sizeof(float) * NumElements * 288);
    float* _buffer1 = (float*)device.api->getStackMemory(sizeof(float) * NumElements * 288);
    _tmp0 = _buffer0;

    sgemm_NT_NT_m21_32_n9_64_k56_nps_3c76782(const_cast<const float *>(rT(0)), 0, const_cast<const float **>(I), ExtraOffset_I, _tmp0, 0, NumElements);
    _tmp1 = _buffer1;

    sgemm_NT_NT_m21_32_n9_32_k21_nss_3b271a8(const_cast<const float *>(fP(1)), 0, const_cast<const float *>(_tmp0), 0, _tmp1, 0, NumElements);
    _tmp2 = _buffer0;

    sgemm_NT_NT_m21_32_n9_9_k9_sps_97c1866(const_cast<const float *>(_tmp1), 0, const_cast<const float **>(AminusT), ExtraOffset_AminusT, _tmp2, 0, NumElements);

    sgemm_NT_NT_m56_64_n9_32_k21_nsp_0e604b1(const_cast<const float *>(rDivM(1)), 0, const_cast<const float *>(_tmp2), 0, Q, ExtraOffset_Q, NumElements);
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
    assert(NumElements != 0);
    device::DeviceInstance& device = device::DeviceInstance::getInstance();
    float *_tmp0, *_tmp1, *_tmp2;
    float* _buffer0 = (float*)device.api->getStackMemory(sizeof(float) * NumElements * 288);
    float* _buffer1 = (float*)device.api->getStackMemory(sizeof(float) * NumElements * 288);
    _tmp0 = _buffer0;

    sgemm_NT_NT_m21_32_n9_64_k56_nps_3c76782(const_cast<const float *>(rT(0)), 0, const_cast<const float **>(I), ExtraOffset_I, _tmp0, 0, NumElements);
    _tmp1 = _buffer1;

    sgemm_NT_NT_m21_32_n9_32_k21_nss_3b271a8(const_cast<const float *>(fP(2)), 0, const_cast<const float *>(_tmp0), 0, _tmp1, 0, NumElements);
    _tmp2 = _buffer0;

    sgemm_NT_NT_m21_32_n9_9_k9_sps_97c1866(const_cast<const float *>(_tmp1), 0, const_cast<const float **>(AminusT), ExtraOffset_AminusT, _tmp2, 0, NumElements);

    sgemm_NT_NT_m56_64_n9_32_k21_nsp_0e604b1(const_cast<const float *>(rDivM(1)), 0, const_cast<const float *>(_tmp2), 0, Q, ExtraOffset_Q, NumElements);
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
    assert(NumElements != 0);
    device::DeviceInstance& device = device::DeviceInstance::getInstance();
    float *_tmp0, *_tmp1, *_tmp2;
    float* _buffer0 = (float*)device.api->getStackMemory(sizeof(float) * NumElements * 288);
    float* _buffer1 = (float*)device.api->getStackMemory(sizeof(float) * NumElements * 288);
    _tmp0 = _buffer0;

    sgemm_NT_NT_m21_32_n9_64_k56_nps_3c76782(const_cast<const float *>(rT(1)), 0, const_cast<const float **>(I), ExtraOffset_I, _tmp0, 0, NumElements);
    _tmp1 = _buffer1;

    sgemm_NT_NT_m21_32_n9_32_k21_nss_3b271a8(const_cast<const float *>(fP(0)), 0, const_cast<const float *>(_tmp0), 0, _tmp1, 0, NumElements);
    _tmp2 = _buffer0;

    sgemm_NT_NT_m21_32_n9_9_k9_sps_97c1866(const_cast<const float *>(_tmp1), 0, const_cast<const float **>(AminusT), ExtraOffset_AminusT, _tmp2, 0, NumElements);

    sgemm_NT_NT_m56_64_n9_32_k21_nsp_0e604b1(const_cast<const float *>(rDivM(1)), 0, const_cast<const float *>(_tmp2), 0, Q, ExtraOffset_Q, NumElements);
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
    assert(NumElements != 0);
    device::DeviceInstance& device = device::DeviceInstance::getInstance();
    float *_tmp0, *_tmp1, *_tmp2;
    float* _buffer0 = (float*)device.api->getStackMemory(sizeof(float) * NumElements * 288);
    float* _buffer1 = (float*)device.api->getStackMemory(sizeof(float) * NumElements * 288);
    _tmp0 = _buffer0;

    sgemm_NT_NT_m21_32_n9_64_k56_nps_3c76782(const_cast<const float *>(rT(1)), 0, const_cast<const float **>(I), ExtraOffset_I, _tmp0, 0, NumElements);
    _tmp1 = _buffer1;

    sgemm_NT_NT_m21_32_n9_32_k21_nss_3b271a8(const_cast<const float *>(fP(1)), 0, const_cast<const float *>(_tmp0), 0, _tmp1, 0, NumElements);
    _tmp2 = _buffer0;

    sgemm_NT_NT_m21_32_n9_9_k9_sps_97c1866(const_cast<const float *>(_tmp1), 0, const_cast<const float **>(AminusT), ExtraOffset_AminusT, _tmp2, 0, NumElements);

    sgemm_NT_NT_m56_64_n9_32_k21_nsp_0e604b1(const_cast<const float *>(rDivM(1)), 0, const_cast<const float *>(_tmp2), 0, Q, ExtraOffset_Q, NumElements);
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
    assert(NumElements != 0);
    device::DeviceInstance& device = device::DeviceInstance::getInstance();
    float *_tmp0, *_tmp1, *_tmp2;
    float* _buffer0 = (float*)device.api->getStackMemory(sizeof(float) * NumElements * 288);
    float* _buffer1 = (float*)device.api->getStackMemory(sizeof(float) * NumElements * 288);
    _tmp0 = _buffer0;

    sgemm_NT_NT_m21_32_n9_64_k56_nps_3c76782(const_cast<const float *>(rT(1)), 0, const_cast<const float **>(I), ExtraOffset_I, _tmp0, 0, NumElements);
    _tmp1 = _buffer1;

    sgemm_NT_NT_m21_32_n9_32_k21_nss_3b271a8(const_cast<const float *>(fP(2)), 0, const_cast<const float *>(_tmp0), 0, _tmp1, 0, NumElements);
    _tmp2 = _buffer0;

    sgemm_NT_NT_m21_32_n9_9_k9_sps_97c1866(const_cast<const float *>(_tmp1), 0, const_cast<const float **>(AminusT), ExtraOffset_AminusT, _tmp2, 0, NumElements);

    sgemm_NT_NT_m56_64_n9_32_k21_nsp_0e604b1(const_cast<const float *>(rDivM(1)), 0, const_cast<const float *>(_tmp2), 0, Q, ExtraOffset_Q, NumElements);
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
    assert(NumElements != 0);
    device::DeviceInstance& device = device::DeviceInstance::getInstance();
    float *_tmp0, *_tmp1, *_tmp2;
    float* _buffer0 = (float*)device.api->getStackMemory(sizeof(float) * NumElements * 288);
    float* _buffer1 = (float*)device.api->getStackMemory(sizeof(float) * NumElements * 288);
    _tmp0 = _buffer0;

    sgemm_NT_NT_m21_32_n9_64_k56_nps_3c76782(const_cast<const float *>(rT(2)), 0, const_cast<const float **>(I), ExtraOffset_I, _tmp0, 0, NumElements);
    _tmp1 = _buffer1;

    sgemm_NT_NT_m21_32_n9_32_k21_nss_3b271a8(const_cast<const float *>(fP(0)), 0, const_cast<const float *>(_tmp0), 0, _tmp1, 0, NumElements);
    _tmp2 = _buffer0;

    sgemm_NT_NT_m21_32_n9_9_k9_sps_97c1866(const_cast<const float *>(_tmp1), 0, const_cast<const float **>(AminusT), ExtraOffset_AminusT, _tmp2, 0, NumElements);

    sgemm_NT_NT_m56_64_n9_32_k21_nsp_0e604b1(const_cast<const float *>(rDivM(1)), 0, const_cast<const float *>(_tmp2), 0, Q, ExtraOffset_Q, NumElements);
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
    assert(NumElements != 0);
    device::DeviceInstance& device = device::DeviceInstance::getInstance();
    float *_tmp0, *_tmp1, *_tmp2;
    float* _buffer0 = (float*)device.api->getStackMemory(sizeof(float) * NumElements * 288);
    float* _buffer1 = (float*)device.api->getStackMemory(sizeof(float) * NumElements * 288);
    _tmp0 = _buffer0;

    sgemm_NT_NT_m21_32_n9_64_k56_nps_3c76782(const_cast<const float *>(rT(2)), 0, const_cast<const float **>(I), ExtraOffset_I, _tmp0, 0, NumElements);
    _tmp1 = _buffer1;

    sgemm_NT_NT_m21_32_n9_32_k21_nss_3b271a8(const_cast<const float *>(fP(1)), 0, const_cast<const float *>(_tmp0), 0, _tmp1, 0, NumElements);
    _tmp2 = _buffer0;

    sgemm_NT_NT_m21_32_n9_9_k9_sps_97c1866(const_cast<const float *>(_tmp1), 0, const_cast<const float **>(AminusT), ExtraOffset_AminusT, _tmp2, 0, NumElements);

    sgemm_NT_NT_m56_64_n9_32_k21_nsp_0e604b1(const_cast<const float *>(rDivM(1)), 0, const_cast<const float *>(_tmp2), 0, Q, ExtraOffset_Q, NumElements);
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
    assert(NumElements != 0);
    device::DeviceInstance& device = device::DeviceInstance::getInstance();
    float *_tmp0, *_tmp1, *_tmp2;
    float* _buffer0 = (float*)device.api->getStackMemory(sizeof(float) * NumElements * 288);
    float* _buffer1 = (float*)device.api->getStackMemory(sizeof(float) * NumElements * 288);
    _tmp0 = _buffer0;

    sgemm_NT_NT_m21_32_n9_64_k56_nps_3c76782(const_cast<const float *>(rT(2)), 0, const_cast<const float **>(I), ExtraOffset_I, _tmp0, 0, NumElements);
    _tmp1 = _buffer1;

    sgemm_NT_NT_m21_32_n9_32_k21_nss_3b271a8(const_cast<const float *>(fP(2)), 0, const_cast<const float *>(_tmp0), 0, _tmp1, 0, NumElements);
    _tmp2 = _buffer0;

    sgemm_NT_NT_m21_32_n9_9_k9_sps_97c1866(const_cast<const float *>(_tmp1), 0, const_cast<const float **>(AminusT), ExtraOffset_AminusT, _tmp2, 0, NumElements);

    sgemm_NT_NT_m56_64_n9_32_k21_nsp_0e604b1(const_cast<const float *>(rDivM(1)), 0, const_cast<const float *>(_tmp2), 0, Q, ExtraOffset_Q, NumElements);
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
    assert(NumElements != 0);
    device::DeviceInstance& device = device::DeviceInstance::getInstance();
    float *_tmp0, *_tmp1, *_tmp2;
    float* _buffer0 = (float*)device.api->getStackMemory(sizeof(float) * NumElements * 288);
    float* _buffer1 = (float*)device.api->getStackMemory(sizeof(float) * NumElements * 288);
    _tmp0 = _buffer0;

    sgemm_NT_NT_m21_32_n9_64_k56_nps_3c76782(const_cast<const float *>(rT(3)), 0, const_cast<const float **>(I), ExtraOffset_I, _tmp0, 0, NumElements);
    _tmp1 = _buffer1;

    sgemm_NT_NT_m21_32_n9_32_k21_nss_3b271a8(const_cast<const float *>(fP(0)), 0, const_cast<const float *>(_tmp0), 0, _tmp1, 0, NumElements);
    _tmp2 = _buffer0;

    sgemm_NT_NT_m21_32_n9_9_k9_sps_97c1866(const_cast<const float *>(_tmp1), 0, const_cast<const float **>(AminusT), ExtraOffset_AminusT, _tmp2, 0, NumElements);

    sgemm_NT_NT_m56_64_n9_32_k21_nsp_0e604b1(const_cast<const float *>(rDivM(1)), 0, const_cast<const float *>(_tmp2), 0, Q, ExtraOffset_Q, NumElements);
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
    assert(NumElements != 0);
    device::DeviceInstance& device = device::DeviceInstance::getInstance();
    float *_tmp0, *_tmp1, *_tmp2;
    float* _buffer0 = (float*)device.api->getStackMemory(sizeof(float) * NumElements * 288);
    float* _buffer1 = (float*)device.api->getStackMemory(sizeof(float) * NumElements * 288);
    _tmp0 = _buffer0;

    sgemm_NT_NT_m21_32_n9_64_k56_nps_3c76782(const_cast<const float *>(rT(3)), 0, const_cast<const float **>(I), ExtraOffset_I, _tmp0, 0, NumElements);
    _tmp1 = _buffer1;

    sgemm_NT_NT_m21_32_n9_32_k21_nss_3b271a8(const_cast<const float *>(fP(1)), 0, const_cast<const float *>(_tmp0), 0, _tmp1, 0, NumElements);
    _tmp2 = _buffer0;

    sgemm_NT_NT_m21_32_n9_9_k9_sps_97c1866(const_cast<const float *>(_tmp1), 0, const_cast<const float **>(AminusT), ExtraOffset_AminusT, _tmp2, 0, NumElements);

    sgemm_NT_NT_m56_64_n9_32_k21_nsp_0e604b1(const_cast<const float *>(rDivM(1)), 0, const_cast<const float *>(_tmp2), 0, Q, ExtraOffset_Q, NumElements);
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
    assert(NumElements != 0);
    device::DeviceInstance& device = device::DeviceInstance::getInstance();
    float *_tmp0, *_tmp1, *_tmp2;
    float* _buffer0 = (float*)device.api->getStackMemory(sizeof(float) * NumElements * 288);
    float* _buffer1 = (float*)device.api->getStackMemory(sizeof(float) * NumElements * 288);
    _tmp0 = _buffer0;

    sgemm_NT_NT_m21_32_n9_64_k56_nps_3c76782(const_cast<const float *>(rT(3)), 0, const_cast<const float **>(I), ExtraOffset_I, _tmp0, 0, NumElements);
    _tmp1 = _buffer1;

    sgemm_NT_NT_m21_32_n9_32_k21_nss_3b271a8(const_cast<const float *>(fP(2)), 0, const_cast<const float *>(_tmp0), 0, _tmp1, 0, NumElements);
    _tmp2 = _buffer0;

    sgemm_NT_NT_m21_32_n9_9_k9_sps_97c1866(const_cast<const float *>(_tmp1), 0, const_cast<const float **>(AminusT), ExtraOffset_AminusT, _tmp2, 0, NumElements);

    sgemm_NT_NT_m56_64_n9_32_k21_nsp_0e604b1(const_cast<const float *>(rDivM(1)), 0, const_cast<const float *>(_tmp2), 0, Q, ExtraOffset_Q, NumElements);
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
    assert(NumElements != 0);
    device::DeviceInstance& device = device::DeviceInstance::getInstance();
    float *_tmp0, *_tmp1, *_tmp2;
    float* _buffer0 = (float*)device.api->getStackMemory(sizeof(float) * NumElements * 288);
    float* _buffer1 = (float*)device.api->getStackMemory(sizeof(float) * NumElements * 288);
    _tmp0 = _buffer0;

    sgemm_NT_NT_m21_32_n9_64_k56_nps_3c76782(const_cast<const float *>(rT(0)), 0, const_cast<const float **>(I), ExtraOffset_I, _tmp0, 0, NumElements);
    _tmp1 = _buffer1;

    sgemm_NT_NT_m21_32_n9_32_k21_nss_3b271a8(const_cast<const float *>(fP(0)), 0, const_cast<const float *>(_tmp0), 0, _tmp1, 0, NumElements);
    _tmp2 = _buffer0;

    sgemm_NT_NT_m21_32_n9_9_k9_sps_97c1866(const_cast<const float *>(_tmp1), 0, const_cast<const float **>(AminusT), ExtraOffset_AminusT, _tmp2, 0, NumElements);

    sgemm_NT_NT_m56_64_n9_32_k21_nsp_0e604b1(const_cast<const float *>(rDivM(2)), 0, const_cast<const float *>(_tmp2), 0, Q, ExtraOffset_Q, NumElements);
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
    assert(NumElements != 0);
    device::DeviceInstance& device = device::DeviceInstance::getInstance();
    float *_tmp0, *_tmp1, *_tmp2;
    float* _buffer0 = (float*)device.api->getStackMemory(sizeof(float) * NumElements * 288);
    float* _buffer1 = (float*)device.api->getStackMemory(sizeof(float) * NumElements * 288);
    _tmp0 = _buffer0;

    sgemm_NT_NT_m21_32_n9_64_k56_nps_3c76782(const_cast<const float *>(rT(0)), 0, const_cast<const float **>(I), ExtraOffset_I, _tmp0, 0, NumElements);
    _tmp1 = _buffer1;

    sgemm_NT_NT_m21_32_n9_32_k21_nss_3b271a8(const_cast<const float *>(fP(1)), 0, const_cast<const float *>(_tmp0), 0, _tmp1, 0, NumElements);
    _tmp2 = _buffer0;

    sgemm_NT_NT_m21_32_n9_9_k9_sps_97c1866(const_cast<const float *>(_tmp1), 0, const_cast<const float **>(AminusT), ExtraOffset_AminusT, _tmp2, 0, NumElements);

    sgemm_NT_NT_m56_64_n9_32_k21_nsp_0e604b1(const_cast<const float *>(rDivM(2)), 0, const_cast<const float *>(_tmp2), 0, Q, ExtraOffset_Q, NumElements);
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
    assert(NumElements != 0);
    device::DeviceInstance& device = device::DeviceInstance::getInstance();
    float *_tmp0, *_tmp1, *_tmp2;
    float* _buffer0 = (float*)device.api->getStackMemory(sizeof(float) * NumElements * 288);
    float* _buffer1 = (float*)device.api->getStackMemory(sizeof(float) * NumElements * 288);
    _tmp0 = _buffer0;

    sgemm_NT_NT_m21_32_n9_64_k56_nps_3c76782(const_cast<const float *>(rT(0)), 0, const_cast<const float **>(I), ExtraOffset_I, _tmp0, 0, NumElements);
    _tmp1 = _buffer1;

    sgemm_NT_NT_m21_32_n9_32_k21_nss_3b271a8(const_cast<const float *>(fP(2)), 0, const_cast<const float *>(_tmp0), 0, _tmp1, 0, NumElements);
    _tmp2 = _buffer0;

    sgemm_NT_NT_m21_32_n9_9_k9_sps_97c1866(const_cast<const float *>(_tmp1), 0, const_cast<const float **>(AminusT), ExtraOffset_AminusT, _tmp2, 0, NumElements);

    sgemm_NT_NT_m56_64_n9_32_k21_nsp_0e604b1(const_cast<const float *>(rDivM(2)), 0, const_cast<const float *>(_tmp2), 0, Q, ExtraOffset_Q, NumElements);
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
    assert(NumElements != 0);
    device::DeviceInstance& device = device::DeviceInstance::getInstance();
    float *_tmp0, *_tmp1, *_tmp2;
    float* _buffer0 = (float*)device.api->getStackMemory(sizeof(float) * NumElements * 288);
    float* _buffer1 = (float*)device.api->getStackMemory(sizeof(float) * NumElements * 288);
    _tmp0 = _buffer0;

    sgemm_NT_NT_m21_32_n9_64_k56_nps_3c76782(const_cast<const float *>(rT(1)), 0, const_cast<const float **>(I), ExtraOffset_I, _tmp0, 0, NumElements);
    _tmp1 = _buffer1;

    sgemm_NT_NT_m21_32_n9_32_k21_nss_3b271a8(const_cast<const float *>(fP(0)), 0, const_cast<const float *>(_tmp0), 0, _tmp1, 0, NumElements);
    _tmp2 = _buffer0;

    sgemm_NT_NT_m21_32_n9_9_k9_sps_97c1866(const_cast<const float *>(_tmp1), 0, const_cast<const float **>(AminusT), ExtraOffset_AminusT, _tmp2, 0, NumElements);

    sgemm_NT_NT_m56_64_n9_32_k21_nsp_0e604b1(const_cast<const float *>(rDivM(2)), 0, const_cast<const float *>(_tmp2), 0, Q, ExtraOffset_Q, NumElements);
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
    assert(NumElements != 0);
    device::DeviceInstance& device = device::DeviceInstance::getInstance();
    float *_tmp0, *_tmp1, *_tmp2;
    float* _buffer0 = (float*)device.api->getStackMemory(sizeof(float) * NumElements * 288);
    float* _buffer1 = (float*)device.api->getStackMemory(sizeof(float) * NumElements * 288);
    _tmp0 = _buffer0;

    sgemm_NT_NT_m21_32_n9_64_k56_nps_3c76782(const_cast<const float *>(rT(1)), 0, const_cast<const float **>(I), ExtraOffset_I, _tmp0, 0, NumElements);
    _tmp1 = _buffer1;

    sgemm_NT_NT_m21_32_n9_32_k21_nss_3b271a8(const_cast<const float *>(fP(1)), 0, const_cast<const float *>(_tmp0), 0, _tmp1, 0, NumElements);
    _tmp2 = _buffer0;

    sgemm_NT_NT_m21_32_n9_9_k9_sps_97c1866(const_cast<const float *>(_tmp1), 0, const_cast<const float **>(AminusT), ExtraOffset_AminusT, _tmp2, 0, NumElements);

    sgemm_NT_NT_m56_64_n9_32_k21_nsp_0e604b1(const_cast<const float *>(rDivM(2)), 0, const_cast<const float *>(_tmp2), 0, Q, ExtraOffset_Q, NumElements);
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
    assert(NumElements != 0);
    device::DeviceInstance& device = device::DeviceInstance::getInstance();
    float *_tmp0, *_tmp1, *_tmp2;
    float* _buffer0 = (float*)device.api->getStackMemory(sizeof(float) * NumElements * 288);
    float* _buffer1 = (float*)device.api->getStackMemory(sizeof(float) * NumElements * 288);
    _tmp0 = _buffer0;

    sgemm_NT_NT_m21_32_n9_64_k56_nps_3c76782(const_cast<const float *>(rT(1)), 0, const_cast<const float **>(I), ExtraOffset_I, _tmp0, 0, NumElements);
    _tmp1 = _buffer1;

    sgemm_NT_NT_m21_32_n9_32_k21_nss_3b271a8(const_cast<const float *>(fP(2)), 0, const_cast<const float *>(_tmp0), 0, _tmp1, 0, NumElements);
    _tmp2 = _buffer0;

    sgemm_NT_NT_m21_32_n9_9_k9_sps_97c1866(const_cast<const float *>(_tmp1), 0, const_cast<const float **>(AminusT), ExtraOffset_AminusT, _tmp2, 0, NumElements);

    sgemm_NT_NT_m56_64_n9_32_k21_nsp_0e604b1(const_cast<const float *>(rDivM(2)), 0, const_cast<const float *>(_tmp2), 0, Q, ExtraOffset_Q, NumElements);
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
    assert(NumElements != 0);
    device::DeviceInstance& device = device::DeviceInstance::getInstance();
    float *_tmp0, *_tmp1, *_tmp2;
    float* _buffer0 = (float*)device.api->getStackMemory(sizeof(float) * NumElements * 288);
    float* _buffer1 = (float*)device.api->getStackMemory(sizeof(float) * NumElements * 288);
    _tmp0 = _buffer0;

    sgemm_NT_NT_m21_32_n9_64_k56_nps_3c76782(const_cast<const float *>(rT(2)), 0, const_cast<const float **>(I), ExtraOffset_I, _tmp0, 0, NumElements);
    _tmp1 = _buffer1;

    sgemm_NT_NT_m21_32_n9_32_k21_nss_3b271a8(const_cast<const float *>(fP(0)), 0, const_cast<const float *>(_tmp0), 0, _tmp1, 0, NumElements);
    _tmp2 = _buffer0;

    sgemm_NT_NT_m21_32_n9_9_k9_sps_97c1866(const_cast<const float *>(_tmp1), 0, const_cast<const float **>(AminusT), ExtraOffset_AminusT, _tmp2, 0, NumElements);

    sgemm_NT_NT_m56_64_n9_32_k21_nsp_0e604b1(const_cast<const float *>(rDivM(2)), 0, const_cast<const float *>(_tmp2), 0, Q, ExtraOffset_Q, NumElements);
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
    assert(NumElements != 0);
    device::DeviceInstance& device = device::DeviceInstance::getInstance();
    float *_tmp0, *_tmp1, *_tmp2;
    float* _buffer0 = (float*)device.api->getStackMemory(sizeof(float) * NumElements * 288);
    float* _buffer1 = (float*)device.api->getStackMemory(sizeof(float) * NumElements * 288);
    _tmp0 = _buffer0;

    sgemm_NT_NT_m21_32_n9_64_k56_nps_3c76782(const_cast<const float *>(rT(2)), 0, const_cast<const float **>(I), ExtraOffset_I, _tmp0, 0, NumElements);
    _tmp1 = _buffer1;

    sgemm_NT_NT_m21_32_n9_32_k21_nss_3b271a8(const_cast<const float *>(fP(1)), 0, const_cast<const float *>(_tmp0), 0, _tmp1, 0, NumElements);
    _tmp2 = _buffer0;

    sgemm_NT_NT_m21_32_n9_9_k9_sps_97c1866(const_cast<const float *>(_tmp1), 0, const_cast<const float **>(AminusT), ExtraOffset_AminusT, _tmp2, 0, NumElements);

    sgemm_NT_NT_m56_64_n9_32_k21_nsp_0e604b1(const_cast<const float *>(rDivM(2)), 0, const_cast<const float *>(_tmp2), 0, Q, ExtraOffset_Q, NumElements);
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
    assert(NumElements != 0);
    device::DeviceInstance& device = device::DeviceInstance::getInstance();
    float *_tmp0, *_tmp1, *_tmp2;
    float* _buffer0 = (float*)device.api->getStackMemory(sizeof(float) * NumElements * 288);
    float* _buffer1 = (float*)device.api->getStackMemory(sizeof(float) * NumElements * 288);
    _tmp0 = _buffer0;

    sgemm_NT_NT_m21_32_n9_64_k56_nps_3c76782(const_cast<const float *>(rT(2)), 0, const_cast<const float **>(I), ExtraOffset_I, _tmp0, 0, NumElements);
    _tmp1 = _buffer1;

    sgemm_NT_NT_m21_32_n9_32_k21_nss_3b271a8(const_cast<const float *>(fP(2)), 0, const_cast<const float *>(_tmp0), 0, _tmp1, 0, NumElements);
    _tmp2 = _buffer0;

    sgemm_NT_NT_m21_32_n9_9_k9_sps_97c1866(const_cast<const float *>(_tmp1), 0, const_cast<const float **>(AminusT), ExtraOffset_AminusT, _tmp2, 0, NumElements);

    sgemm_NT_NT_m56_64_n9_32_k21_nsp_0e604b1(const_cast<const float *>(rDivM(2)), 0, const_cast<const float *>(_tmp2), 0, Q, ExtraOffset_Q, NumElements);
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
    assert(NumElements != 0);
    device::DeviceInstance& device = device::DeviceInstance::getInstance();
    float *_tmp0, *_tmp1, *_tmp2;
    float* _buffer0 = (float*)device.api->getStackMemory(sizeof(float) * NumElements * 288);
    float* _buffer1 = (float*)device.api->getStackMemory(sizeof(float) * NumElements * 288);
    _tmp0 = _buffer0;

    sgemm_NT_NT_m21_32_n9_64_k56_nps_3c76782(const_cast<const float *>(rT(3)), 0, const_cast<const float **>(I), ExtraOffset_I, _tmp0, 0, NumElements);
    _tmp1 = _buffer1;

    sgemm_NT_NT_m21_32_n9_32_k21_nss_3b271a8(const_cast<const float *>(fP(0)), 0, const_cast<const float *>(_tmp0), 0, _tmp1, 0, NumElements);
    _tmp2 = _buffer0;

    sgemm_NT_NT_m21_32_n9_9_k9_sps_97c1866(const_cast<const float *>(_tmp1), 0, const_cast<const float **>(AminusT), ExtraOffset_AminusT, _tmp2, 0, NumElements);

    sgemm_NT_NT_m56_64_n9_32_k21_nsp_0e604b1(const_cast<const float *>(rDivM(2)), 0, const_cast<const float *>(_tmp2), 0, Q, ExtraOffset_Q, NumElements);
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
    assert(NumElements != 0);
    device::DeviceInstance& device = device::DeviceInstance::getInstance();
    float *_tmp0, *_tmp1, *_tmp2;
    float* _buffer0 = (float*)device.api->getStackMemory(sizeof(float) * NumElements * 288);
    float* _buffer1 = (float*)device.api->getStackMemory(sizeof(float) * NumElements * 288);
    _tmp0 = _buffer0;

    sgemm_NT_NT_m21_32_n9_64_k56_nps_3c76782(const_cast<const float *>(rT(3)), 0, const_cast<const float **>(I), ExtraOffset_I, _tmp0, 0, NumElements);
    _tmp1 = _buffer1;

    sgemm_NT_NT_m21_32_n9_32_k21_nss_3b271a8(const_cast<const float *>(fP(1)), 0, const_cast<const float *>(_tmp0), 0, _tmp1, 0, NumElements);
    _tmp2 = _buffer0;

    sgemm_NT_NT_m21_32_n9_9_k9_sps_97c1866(const_cast<const float *>(_tmp1), 0, const_cast<const float **>(AminusT), ExtraOffset_AminusT, _tmp2, 0, NumElements);

    sgemm_NT_NT_m56_64_n9_32_k21_nsp_0e604b1(const_cast<const float *>(rDivM(2)), 0, const_cast<const float *>(_tmp2), 0, Q, ExtraOffset_Q, NumElements);
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
    assert(NumElements != 0);
    device::DeviceInstance& device = device::DeviceInstance::getInstance();
    float *_tmp0, *_tmp1, *_tmp2;
    float* _buffer0 = (float*)device.api->getStackMemory(sizeof(float) * NumElements * 288);
    float* _buffer1 = (float*)device.api->getStackMemory(sizeof(float) * NumElements * 288);
    _tmp0 = _buffer0;

    sgemm_NT_NT_m21_32_n9_64_k56_nps_3c76782(const_cast<const float *>(rT(3)), 0, const_cast<const float **>(I), ExtraOffset_I, _tmp0, 0, NumElements);
    _tmp1 = _buffer1;

    sgemm_NT_NT_m21_32_n9_32_k21_nss_3b271a8(const_cast<const float *>(fP(2)), 0, const_cast<const float *>(_tmp0), 0, _tmp1, 0, NumElements);
    _tmp2 = _buffer0;

    sgemm_NT_NT_m21_32_n9_9_k9_sps_97c1866(const_cast<const float *>(_tmp1), 0, const_cast<const float **>(AminusT), ExtraOffset_AminusT, _tmp2, 0, NumElements);

    sgemm_NT_NT_m56_64_n9_32_k21_nsp_0e604b1(const_cast<const float *>(rDivM(2)), 0, const_cast<const float *>(_tmp2), 0, Q, ExtraOffset_Q, NumElements);
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
    assert(NumElements != 0);
    device::DeviceInstance& device = device::DeviceInstance::getInstance();
    float *_tmp0, *_tmp1, *_tmp2;
    float* _buffer0 = (float*)device.api->getStackMemory(sizeof(float) * NumElements * 288);
    float* _buffer1 = (float*)device.api->getStackMemory(sizeof(float) * NumElements * 288);
    _tmp0 = _buffer0;

    sgemm_NT_NT_m21_32_n9_64_k56_nps_3c76782(const_cast<const float *>(rT(0)), 0, const_cast<const float **>(I), ExtraOffset_I, _tmp0, 0, NumElements);
    _tmp1 = _buffer1;

    sgemm_NT_NT_m21_32_n9_32_k21_nss_3b271a8(const_cast<const float *>(fP(0)), 0, const_cast<const float *>(_tmp0), 0, _tmp1, 0, NumElements);
    _tmp2 = _buffer0;

    sgemm_NT_NT_m21_32_n9_9_k9_sps_97c1866(const_cast<const float *>(_tmp1), 0, const_cast<const float **>(AminusT), ExtraOffset_AminusT, _tmp2, 0, NumElements);

    sgemm_NT_NT_m56_64_n9_32_k21_nsp_0e604b1(const_cast<const float *>(rDivM(3)), 0, const_cast<const float *>(_tmp2), 0, Q, ExtraOffset_Q, NumElements);
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
    assert(NumElements != 0);
    device::DeviceInstance& device = device::DeviceInstance::getInstance();
    float *_tmp0, *_tmp1, *_tmp2;
    float* _buffer0 = (float*)device.api->getStackMemory(sizeof(float) * NumElements * 288);
    float* _buffer1 = (float*)device.api->getStackMemory(sizeof(float) * NumElements * 288);
    _tmp0 = _buffer0;

    sgemm_NT_NT_m21_32_n9_64_k56_nps_3c76782(const_cast<const float *>(rT(0)), 0, const_cast<const float **>(I), ExtraOffset_I, _tmp0, 0, NumElements);
    _tmp1 = _buffer1;

    sgemm_NT_NT_m21_32_n9_32_k21_nss_3b271a8(const_cast<const float *>(fP(1)), 0, const_cast<const float *>(_tmp0), 0, _tmp1, 0, NumElements);
    _tmp2 = _buffer0;

    sgemm_NT_NT_m21_32_n9_9_k9_sps_97c1866(const_cast<const float *>(_tmp1), 0, const_cast<const float **>(AminusT), ExtraOffset_AminusT, _tmp2, 0, NumElements);

    sgemm_NT_NT_m56_64_n9_32_k21_nsp_0e604b1(const_cast<const float *>(rDivM(3)), 0, const_cast<const float *>(_tmp2), 0, Q, ExtraOffset_Q, NumElements);
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
    assert(NumElements != 0);
    device::DeviceInstance& device = device::DeviceInstance::getInstance();
    float *_tmp0, *_tmp1, *_tmp2;
    float* _buffer0 = (float*)device.api->getStackMemory(sizeof(float) * NumElements * 288);
    float* _buffer1 = (float*)device.api->getStackMemory(sizeof(float) * NumElements * 288);
    _tmp0 = _buffer0;

    sgemm_NT_NT_m21_32_n9_64_k56_nps_3c76782(const_cast<const float *>(rT(0)), 0, const_cast<const float **>(I), ExtraOffset_I, _tmp0, 0, NumElements);
    _tmp1 = _buffer1;

    sgemm_NT_NT_m21_32_n9_32_k21_nss_3b271a8(const_cast<const float *>(fP(2)), 0, const_cast<const float *>(_tmp0), 0, _tmp1, 0, NumElements);
    _tmp2 = _buffer0;

    sgemm_NT_NT_m21_32_n9_9_k9_sps_97c1866(const_cast<const float *>(_tmp1), 0, const_cast<const float **>(AminusT), ExtraOffset_AminusT, _tmp2, 0, NumElements);

    sgemm_NT_NT_m56_64_n9_32_k21_nsp_0e604b1(const_cast<const float *>(rDivM(3)), 0, const_cast<const float *>(_tmp2), 0, Q, ExtraOffset_Q, NumElements);
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
    assert(NumElements != 0);
    device::DeviceInstance& device = device::DeviceInstance::getInstance();
    float *_tmp0, *_tmp1, *_tmp2;
    float* _buffer0 = (float*)device.api->getStackMemory(sizeof(float) * NumElements * 288);
    float* _buffer1 = (float*)device.api->getStackMemory(sizeof(float) * NumElements * 288);
    _tmp0 = _buffer0;

    sgemm_NT_NT_m21_32_n9_64_k56_nps_3c76782(const_cast<const float *>(rT(1)), 0, const_cast<const float **>(I), ExtraOffset_I, _tmp0, 0, NumElements);
    _tmp1 = _buffer1;

    sgemm_NT_NT_m21_32_n9_32_k21_nss_3b271a8(const_cast<const float *>(fP(0)), 0, const_cast<const float *>(_tmp0), 0, _tmp1, 0, NumElements);
    _tmp2 = _buffer0;

    sgemm_NT_NT_m21_32_n9_9_k9_sps_97c1866(const_cast<const float *>(_tmp1), 0, const_cast<const float **>(AminusT), ExtraOffset_AminusT, _tmp2, 0, NumElements);

    sgemm_NT_NT_m56_64_n9_32_k21_nsp_0e604b1(const_cast<const float *>(rDivM(3)), 0, const_cast<const float *>(_tmp2), 0, Q, ExtraOffset_Q, NumElements);
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
    assert(NumElements != 0);
    device::DeviceInstance& device = device::DeviceInstance::getInstance();
    float *_tmp0, *_tmp1, *_tmp2;
    float* _buffer0 = (float*)device.api->getStackMemory(sizeof(float) * NumElements * 288);
    float* _buffer1 = (float*)device.api->getStackMemory(sizeof(float) * NumElements * 288);
    _tmp0 = _buffer0;

    sgemm_NT_NT_m21_32_n9_64_k56_nps_3c76782(const_cast<const float *>(rT(1)), 0, const_cast<const float **>(I), ExtraOffset_I, _tmp0, 0, NumElements);
    _tmp1 = _buffer1;

    sgemm_NT_NT_m21_32_n9_32_k21_nss_3b271a8(const_cast<const float *>(fP(1)), 0, const_cast<const float *>(_tmp0), 0, _tmp1, 0, NumElements);
    _tmp2 = _buffer0;

    sgemm_NT_NT_m21_32_n9_9_k9_sps_97c1866(const_cast<const float *>(_tmp1), 0, const_cast<const float **>(AminusT), ExtraOffset_AminusT, _tmp2, 0, NumElements);

    sgemm_NT_NT_m56_64_n9_32_k21_nsp_0e604b1(const_cast<const float *>(rDivM(3)), 0, const_cast<const float *>(_tmp2), 0, Q, ExtraOffset_Q, NumElements);
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
    assert(NumElements != 0);
    device::DeviceInstance& device = device::DeviceInstance::getInstance();
    float *_tmp0, *_tmp1, *_tmp2;
    float* _buffer0 = (float*)device.api->getStackMemory(sizeof(float) * NumElements * 288);
    float* _buffer1 = (float*)device.api->getStackMemory(sizeof(float) * NumElements * 288);
    _tmp0 = _buffer0;

    sgemm_NT_NT_m21_32_n9_64_k56_nps_3c76782(const_cast<const float *>(rT(1)), 0, const_cast<const float **>(I), ExtraOffset_I, _tmp0, 0, NumElements);
    _tmp1 = _buffer1;

    sgemm_NT_NT_m21_32_n9_32_k21_nss_3b271a8(const_cast<const float *>(fP(2)), 0, const_cast<const float *>(_tmp0), 0, _tmp1, 0, NumElements);
    _tmp2 = _buffer0;

    sgemm_NT_NT_m21_32_n9_9_k9_sps_97c1866(const_cast<const float *>(_tmp1), 0, const_cast<const float **>(AminusT), ExtraOffset_AminusT, _tmp2, 0, NumElements);

    sgemm_NT_NT_m56_64_n9_32_k21_nsp_0e604b1(const_cast<const float *>(rDivM(3)), 0, const_cast<const float *>(_tmp2), 0, Q, ExtraOffset_Q, NumElements);
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
    assert(NumElements != 0);
    device::DeviceInstance& device = device::DeviceInstance::getInstance();
    float *_tmp0, *_tmp1, *_tmp2;
    float* _buffer0 = (float*)device.api->getStackMemory(sizeof(float) * NumElements * 288);
    float* _buffer1 = (float*)device.api->getStackMemory(sizeof(float) * NumElements * 288);
    _tmp0 = _buffer0;

    sgemm_NT_NT_m21_32_n9_64_k56_nps_3c76782(const_cast<const float *>(rT(2)), 0, const_cast<const float **>(I), ExtraOffset_I, _tmp0, 0, NumElements);
    _tmp1 = _buffer1;

    sgemm_NT_NT_m21_32_n9_32_k21_nss_3b271a8(const_cast<const float *>(fP(0)), 0, const_cast<const float *>(_tmp0), 0, _tmp1, 0, NumElements);
    _tmp2 = _buffer0;

    sgemm_NT_NT_m21_32_n9_9_k9_sps_97c1866(const_cast<const float *>(_tmp1), 0, const_cast<const float **>(AminusT), ExtraOffset_AminusT, _tmp2, 0, NumElements);

    sgemm_NT_NT_m56_64_n9_32_k21_nsp_0e604b1(const_cast<const float *>(rDivM(3)), 0, const_cast<const float *>(_tmp2), 0, Q, ExtraOffset_Q, NumElements);
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
    assert(NumElements != 0);
    device::DeviceInstance& device = device::DeviceInstance::getInstance();
    float *_tmp0, *_tmp1, *_tmp2;
    float* _buffer0 = (float*)device.api->getStackMemory(sizeof(float) * NumElements * 288);
    float* _buffer1 = (float*)device.api->getStackMemory(sizeof(float) * NumElements * 288);
    _tmp0 = _buffer0;

    sgemm_NT_NT_m21_32_n9_64_k56_nps_3c76782(const_cast<const float *>(rT(2)), 0, const_cast<const float **>(I), ExtraOffset_I, _tmp0, 0, NumElements);
    _tmp1 = _buffer1;

    sgemm_NT_NT_m21_32_n9_32_k21_nss_3b271a8(const_cast<const float *>(fP(1)), 0, const_cast<const float *>(_tmp0), 0, _tmp1, 0, NumElements);
    _tmp2 = _buffer0;

    sgemm_NT_NT_m21_32_n9_9_k9_sps_97c1866(const_cast<const float *>(_tmp1), 0, const_cast<const float **>(AminusT), ExtraOffset_AminusT, _tmp2, 0, NumElements);

    sgemm_NT_NT_m56_64_n9_32_k21_nsp_0e604b1(const_cast<const float *>(rDivM(3)), 0, const_cast<const float *>(_tmp2), 0, Q, ExtraOffset_Q, NumElements);
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
    assert(NumElements != 0);
    device::DeviceInstance& device = device::DeviceInstance::getInstance();
    float *_tmp0, *_tmp1, *_tmp2;
    float* _buffer0 = (float*)device.api->getStackMemory(sizeof(float) * NumElements * 288);
    float* _buffer1 = (float*)device.api->getStackMemory(sizeof(float) * NumElements * 288);
    _tmp0 = _buffer0;

    sgemm_NT_NT_m21_32_n9_64_k56_nps_3c76782(const_cast<const float *>(rT(2)), 0, const_cast<const float **>(I), ExtraOffset_I, _tmp0, 0, NumElements);
    _tmp1 = _buffer1;

    sgemm_NT_NT_m21_32_n9_32_k21_nss_3b271a8(const_cast<const float *>(fP(2)), 0, const_cast<const float *>(_tmp0), 0, _tmp1, 0, NumElements);
    _tmp2 = _buffer0;

    sgemm_NT_NT_m21_32_n9_9_k9_sps_97c1866(const_cast<const float *>(_tmp1), 0, const_cast<const float **>(AminusT), ExtraOffset_AminusT, _tmp2, 0, NumElements);

    sgemm_NT_NT_m56_64_n9_32_k21_nsp_0e604b1(const_cast<const float *>(rDivM(3)), 0, const_cast<const float *>(_tmp2), 0, Q, ExtraOffset_Q, NumElements);
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
    assert(NumElements != 0);
    device::DeviceInstance& device = device::DeviceInstance::getInstance();
    float *_tmp0, *_tmp1, *_tmp2;
    float* _buffer0 = (float*)device.api->getStackMemory(sizeof(float) * NumElements * 288);
    float* _buffer1 = (float*)device.api->getStackMemory(sizeof(float) * NumElements * 288);
    _tmp0 = _buffer0;

    sgemm_NT_NT_m21_32_n9_64_k56_nps_3c76782(const_cast<const float *>(rT(3)), 0, const_cast<const float **>(I), ExtraOffset_I, _tmp0, 0, NumElements);
    _tmp1 = _buffer1;

    sgemm_NT_NT_m21_32_n9_32_k21_nss_3b271a8(const_cast<const float *>(fP(0)), 0, const_cast<const float *>(_tmp0), 0, _tmp1, 0, NumElements);
    _tmp2 = _buffer0;

    sgemm_NT_NT_m21_32_n9_9_k9_sps_97c1866(const_cast<const float *>(_tmp1), 0, const_cast<const float **>(AminusT), ExtraOffset_AminusT, _tmp2, 0, NumElements);

    sgemm_NT_NT_m56_64_n9_32_k21_nsp_0e604b1(const_cast<const float *>(rDivM(3)), 0, const_cast<const float *>(_tmp2), 0, Q, ExtraOffset_Q, NumElements);
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
    assert(NumElements != 0);
    device::DeviceInstance& device = device::DeviceInstance::getInstance();
    float *_tmp0, *_tmp1, *_tmp2;
    float* _buffer0 = (float*)device.api->getStackMemory(sizeof(float) * NumElements * 288);
    float* _buffer1 = (float*)device.api->getStackMemory(sizeof(float) * NumElements * 288);
    _tmp0 = _buffer0;

    sgemm_NT_NT_m21_32_n9_64_k56_nps_3c76782(const_cast<const float *>(rT(3)), 0, const_cast<const float **>(I), ExtraOffset_I, _tmp0, 0, NumElements);
    _tmp1 = _buffer1;

    sgemm_NT_NT_m21_32_n9_32_k21_nss_3b271a8(const_cast<const float *>(fP(1)), 0, const_cast<const float *>(_tmp0), 0, _tmp1, 0, NumElements);
    _tmp2 = _buffer0;

    sgemm_NT_NT_m21_32_n9_9_k9_sps_97c1866(const_cast<const float *>(_tmp1), 0, const_cast<const float **>(AminusT), ExtraOffset_AminusT, _tmp2, 0, NumElements);

    sgemm_NT_NT_m56_64_n9_32_k21_nsp_0e604b1(const_cast<const float *>(rDivM(3)), 0, const_cast<const float *>(_tmp2), 0, Q, ExtraOffset_Q, NumElements);
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
    assert(NumElements != 0);
    device::DeviceInstance& device = device::DeviceInstance::getInstance();
    float *_tmp0, *_tmp1, *_tmp2;
    float* _buffer0 = (float*)device.api->getStackMemory(sizeof(float) * NumElements * 288);
    float* _buffer1 = (float*)device.api->getStackMemory(sizeof(float) * NumElements * 288);
    _tmp0 = _buffer0;

    sgemm_NT_NT_m21_32_n9_64_k56_nps_3c76782(const_cast<const float *>(rT(3)), 0, const_cast<const float **>(I), ExtraOffset_I, _tmp0, 0, NumElements);
    _tmp1 = _buffer1;

    sgemm_NT_NT_m21_32_n9_32_k21_nss_3b271a8(const_cast<const float *>(fP(2)), 0, const_cast<const float *>(_tmp0), 0, _tmp1, 0, NumElements);
    _tmp2 = _buffer0;

    sgemm_NT_NT_m21_32_n9_9_k9_sps_97c1866(const_cast<const float *>(_tmp1), 0, const_cast<const float **>(AminusT), ExtraOffset_AminusT, _tmp2, 0, NumElements);

    sgemm_NT_NT_m56_64_n9_32_k21_nsp_0e604b1(const_cast<const float *>(rDivM(3)), 0, const_cast<const float *>(_tmp2), 0, Q, ExtraOffset_Q, NumElements);
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
    assert(NumElements != 0);
    device::DeviceInstance& device = device::DeviceInstance::getInstance();
    scopyAddScale_m56_64_n9_64_pp_1869608(power, dQ(0), ExtraOffset_dQ(0), I, ExtraOffset_I, NumElements);
  }
  void kernel::derivativeTaylorExpansion::execute1() {
    assert(!std::isnan(power));
    assert(I != nullptr);
    assert(dQ(1) != nullptr);
    assert(NumElements != 0);
    device::DeviceInstance& device = device::DeviceInstance::getInstance();
    scopyAddScale_m35_64_n9_64_pp_8aadec7(power, dQ(1), ExtraOffset_dQ(1), I, ExtraOffset_I, NumElements);
  }
  void kernel::derivativeTaylorExpansion::execute2() {
    assert(!std::isnan(power));
    assert(I != nullptr);
    assert(dQ(2) != nullptr);
    assert(NumElements != 0);
    device::DeviceInstance& device = device::DeviceInstance::getInstance();
    scopyAddScale_m20_64_n9_64_pp_f68fa48(power, dQ(2), ExtraOffset_dQ(2), I, ExtraOffset_I, NumElements);
  }
  void kernel::derivativeTaylorExpansion::execute3() {
    assert(!std::isnan(power));
    assert(I != nullptr);
    assert(dQ(3) != nullptr);
    assert(NumElements != 0);
    device::DeviceInstance& device = device::DeviceInstance::getInstance();
    scopyAddScale_m10_64_n9_64_pp_5f5bfb7(power, dQ(3), ExtraOffset_dQ(3), I, ExtraOffset_I, NumElements);
  }
  void kernel::derivativeTaylorExpansion::execute4() {
    assert(!std::isnan(power));
    assert(I != nullptr);
    assert(dQ(4) != nullptr);
    assert(NumElements != 0);
    device::DeviceInstance& device = device::DeviceInstance::getInstance();
    scopyAddScale_m4_64_n9_64_pp_8b6c4ec(power, dQ(4), ExtraOffset_dQ(4), I, ExtraOffset_I, NumElements);
  }
  void kernel::derivativeTaylorExpansion::execute5() {
    assert(!std::isnan(power));
    assert(I != nullptr);
    assert(dQ(5) != nullptr);
    assert(NumElements != 0);
    device::DeviceInstance& device = device::DeviceInstance::getInstance();
    scopyAddScale_m1_64_n9_64_pp_ddfc1ea(power, dQ(5), ExtraOffset_dQ(5), I, ExtraOffset_I, NumElements);
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
    assert(NumElements != 0);
    device::DeviceInstance& device = device::DeviceInstance::getInstance();
    float *_tmp0, *_tmp2, *_tmp4;
    float* _buffer0 = (float*)device.api->getStackMemory(sizeof(float) * NumElements * 432);
    _tmp0 = _buffer0;

    sgemm_NT_NT_m35_48_n9_64_k53_nps_7539a56(const_cast<const float *>(kDivMT(0)), 0, const_cast<const float **>(dQ(0)), ExtraOffset_dQ(0), _tmp0, 0, NumElements);

    sgemm_NT_NT_m35_48_n9_9_k9_spp_bcf4f26(const_cast<const float *>(_tmp0), 0, const_cast<const float **>(star(0)), ExtraOffset_star(0), dQ(1), ExtraOffset_dQ(1), NumElements);
    _tmp2 = _buffer0;

    sgemm_NT_NT_m35_48_n9_64_k54_nps_bfcf3fd(const_cast<const float *>(kDivMT(1)), 0, const_cast<const float **>(dQ(0)), ExtraOffset_dQ(0), _tmp2, 0, NumElements);

    sgemm_NT_NT_m35_48_n9_9_k9_spp_d135d69(const_cast<const float *>(_tmp2), 0, const_cast<const float **>(star(1)), ExtraOffset_star(1), dQ(1), ExtraOffset_dQ(1), NumElements);
    _tmp4 = _buffer0;

    sgemm_NT_NT_m35_48_n9_64_k55_nps_72bd187(const_cast<const float *>(kDivMT(2)), 0, const_cast<const float **>(dQ(0)), ExtraOffset_dQ(0), _tmp4, 0, NumElements);

    sgemm_NT_NT_m35_48_n9_9_k9_spp_d135d69(const_cast<const float *>(_tmp4), 0, const_cast<const float **>(star(2)), ExtraOffset_star(2), dQ(1), ExtraOffset_dQ(1), NumElements);
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
    assert(NumElements != 0);
    device::DeviceInstance& device = device::DeviceInstance::getInstance();
    float *_tmp0, *_tmp2, *_tmp4;
    float* _buffer0 = (float*)device.api->getStackMemory(sizeof(float) * NumElements * 288);
    _tmp0 = _buffer0;

    sgemm_NT_NT_m20_48_n9_48_k32_nps_c879df6(const_cast<const float *>(kDivMT(0)), 0, const_cast<const float **>(dQ(1)), ExtraOffset_dQ(1), _tmp0, 0, NumElements);

    sgemm_NT_NT_m20_32_n9_9_k9_spp_7cafb26(const_cast<const float *>(_tmp0), 0, const_cast<const float **>(star(0)), ExtraOffset_star(0), dQ(2), ExtraOffset_dQ(2), NumElements);
    _tmp2 = _buffer0;

    sgemm_NT_NT_m20_48_n9_48_k33_nps_34e9903(const_cast<const float *>(kDivMT(1)), 0, const_cast<const float **>(dQ(1)), ExtraOffset_dQ(1), _tmp2, 0, NumElements);

    sgemm_NT_NT_m20_32_n9_9_k9_spp_51f73de(const_cast<const float *>(_tmp2), 0, const_cast<const float **>(star(1)), ExtraOffset_star(1), dQ(2), ExtraOffset_dQ(2), NumElements);
    _tmp4 = _buffer0;

    sgemm_NT_NT_m20_48_n9_48_k34_nps_c44f1e3(const_cast<const float *>(kDivMT(2)), 0, const_cast<const float **>(dQ(1)), ExtraOffset_dQ(1), _tmp4, 0, NumElements);

    sgemm_NT_NT_m20_32_n9_9_k9_spp_51f73de(const_cast<const float *>(_tmp4), 0, const_cast<const float **>(star(2)), ExtraOffset_star(2), dQ(2), ExtraOffset_dQ(2), NumElements);
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
    assert(NumElements != 0);
    device::DeviceInstance& device = device::DeviceInstance::getInstance();
    float *_tmp0, *_tmp2, *_tmp4;
    float* _buffer0 = (float*)device.api->getStackMemory(sizeof(float) * NumElements * 144);
    _tmp0 = _buffer0;

    sgemm_NT_NT_m10_48_n9_32_k17_nps_eeac906(const_cast<const float *>(kDivMT(0)), 0, const_cast<const float **>(dQ(2)), ExtraOffset_dQ(2), _tmp0, 0, NumElements);

    sgemm_NT_NT_m10_16_n9_9_k9_spp_afbfa38(const_cast<const float *>(_tmp0), 0, const_cast<const float **>(star(0)), ExtraOffset_star(0), dQ(3), ExtraOffset_dQ(3), NumElements);
    _tmp2 = _buffer0;

    sgemm_NT_NT_m10_48_n9_32_k18_nps_2125e35(const_cast<const float *>(kDivMT(1)), 0, const_cast<const float **>(dQ(2)), ExtraOffset_dQ(2), _tmp2, 0, NumElements);

    sgemm_NT_NT_m10_16_n9_9_k9_spp_fe81cdb(const_cast<const float *>(_tmp2), 0, const_cast<const float **>(star(1)), ExtraOffset_star(1), dQ(3), ExtraOffset_dQ(3), NumElements);
    _tmp4 = _buffer0;

    sgemm_NT_NT_m10_48_n9_32_k19_nps_63bce85(const_cast<const float *>(kDivMT(2)), 0, const_cast<const float **>(dQ(2)), ExtraOffset_dQ(2), _tmp4, 0, NumElements);

    sgemm_NT_NT_m10_16_n9_9_k9_spp_fe81cdb(const_cast<const float *>(_tmp4), 0, const_cast<const float **>(star(2)), ExtraOffset_star(2), dQ(3), ExtraOffset_dQ(3), NumElements);
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
    assert(NumElements != 0);
    device::DeviceInstance& device = device::DeviceInstance::getInstance();
    float *_tmp0, *_tmp2, *_tmp4;
    float* _buffer0 = (float*)device.api->getStackMemory(sizeof(float) * NumElements * 144);
    _tmp0 = _buffer0;

    sgemm_NT_NT_m4_48_n9_16_k7_nps_dff4b53(const_cast<const float *>(kDivMT(0)), 0, const_cast<const float **>(dQ(3)), ExtraOffset_dQ(3), _tmp0, 0, NumElements);

    sgemm_NT_NT_m4_16_n9_9_k9_spp_b787603(const_cast<const float *>(_tmp0), 0, const_cast<const float **>(star(0)), ExtraOffset_star(0), dQ(4), ExtraOffset_dQ(4), NumElements);
    _tmp2 = _buffer0;

    sgemm_NT_NT_m4_48_n9_16_k8_nps_7c6bd3e(const_cast<const float *>(kDivMT(1)), 0, const_cast<const float **>(dQ(3)), ExtraOffset_dQ(3), _tmp2, 0, NumElements);

    sgemm_NT_NT_m4_16_n9_9_k9_spp_fbdc4fe(const_cast<const float *>(_tmp2), 0, const_cast<const float **>(star(1)), ExtraOffset_star(1), dQ(4), ExtraOffset_dQ(4), NumElements);
    _tmp4 = _buffer0;

    sgemm_NT_NT_m4_48_n9_16_k9_nps_83bc4d7(const_cast<const float *>(kDivMT(2)), 0, const_cast<const float **>(dQ(3)), ExtraOffset_dQ(3), _tmp4, 0, NumElements);

    sgemm_NT_NT_m4_16_n9_9_k9_spp_fbdc4fe(const_cast<const float *>(_tmp4), 0, const_cast<const float **>(star(2)), ExtraOffset_star(2), dQ(4), ExtraOffset_dQ(4), NumElements);
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
    assert(NumElements != 0);
    device::DeviceInstance& device = device::DeviceInstance::getInstance();
    float *_tmp0, *_tmp2, *_tmp4;
    float* _buffer0 = (float*)device.api->getStackMemory(sizeof(float) * NumElements * 144);
    _tmp0 = _buffer0;

    sgemm_NT_NT_m1_48_n9_16_k1_nps_7a8cc39(const_cast<const float *>(kDivMT(0)), 0, const_cast<const float **>(dQ(4)), ExtraOffset_dQ(4), _tmp0, 0, NumElements);

    sgemm_NT_NT_m1_16_n9_9_k9_spp_8da6d94(const_cast<const float *>(_tmp0), 0, const_cast<const float **>(star(0)), ExtraOffset_star(0), dQ(5), ExtraOffset_dQ(5), NumElements);
    _tmp2 = _buffer0;

    sgemm_NT_NT_m1_48_n9_16_k2_nps_04350e1(const_cast<const float *>(kDivMT(1)), 0, const_cast<const float **>(dQ(4)), ExtraOffset_dQ(4), _tmp2, 0, NumElements);

    sgemm_NT_NT_m1_16_n9_9_k9_spp_fe6bf6f(const_cast<const float *>(_tmp2), 0, const_cast<const float **>(star(1)), ExtraOffset_star(1), dQ(5), ExtraOffset_dQ(5), NumElements);
    _tmp4 = _buffer0;

    sgemm_NT_NT_m1_48_n9_16_k3_nps_c033cdc(const_cast<const float *>(kDivMT(2)), 0, const_cast<const float **>(dQ(4)), ExtraOffset_dQ(4), _tmp4, 0, NumElements);

    sgemm_NT_NT_m1_16_n9_9_k9_spp_fe6bf6f(const_cast<const float *>(_tmp4), 0, const_cast<const float **>(star(2)), ExtraOffset_star(2), dQ(5), ExtraOffset_dQ(5), NumElements);
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
    assert(NumElements != 0);
    device::DeviceInstance& device = device::DeviceInstance::getInstance();
    float *_tmp0;
    float* _buffer0 = (float*)device.api->getStackMemory(sizeof(float) * NumElements * 576);
    _tmp0 = _buffer0;

    sgemm_NT_NT_m49_64_n9_64_k56_nps_29706f6(const_cast<const float *>(V3mTo2n(0,0)), 0, const_cast<const float **>(Q), ExtraOffset_Q, _tmp0, 0, NumElements);

    sgemm_NT_NT_m49_64_n9_9_k9_spp_43ea6a7(const_cast<const float *>(_tmp0), 0, const_cast<const float **>(godunovMatrix), ExtraOffset_godunovMatrix, godunovState, ExtraOffset_godunovState, NumElements);
    device.api->popStackMemory();
  }
  void kernel::godunovState::execute1() {
    assert(Q != nullptr);
    assert(V3mTo2n(1,0) != nullptr);
    assert(godunovMatrix != nullptr);
    assert(godunovState != nullptr);
    assert(NumElements != 0);
    device::DeviceInstance& device = device::DeviceInstance::getInstance();
    float *_tmp0;
    float* _buffer0 = (float*)device.api->getStackMemory(sizeof(float) * NumElements * 576);
    _tmp0 = _buffer0;

    sgemm_NT_NT_m49_64_n9_64_k56_nps_29706f6(const_cast<const float *>(V3mTo2n(1,0)), 0, const_cast<const float **>(Q), ExtraOffset_Q, _tmp0, 0, NumElements);

    sgemm_NT_NT_m49_64_n9_9_k9_spp_43ea6a7(const_cast<const float *>(_tmp0), 0, const_cast<const float **>(godunovMatrix), ExtraOffset_godunovMatrix, godunovState, ExtraOffset_godunovState, NumElements);
    device.api->popStackMemory();
  }
  void kernel::godunovState::execute2() {
    assert(Q != nullptr);
    assert(V3mTo2n(2,0) != nullptr);
    assert(godunovMatrix != nullptr);
    assert(godunovState != nullptr);
    assert(NumElements != 0);
    device::DeviceInstance& device = device::DeviceInstance::getInstance();
    float *_tmp0;
    float* _buffer0 = (float*)device.api->getStackMemory(sizeof(float) * NumElements * 576);
    _tmp0 = _buffer0;

    sgemm_NT_NT_m49_64_n9_64_k56_nps_29706f6(const_cast<const float *>(V3mTo2n(2,0)), 0, const_cast<const float **>(Q), ExtraOffset_Q, _tmp0, 0, NumElements);

    sgemm_NT_NT_m49_64_n9_9_k9_spp_43ea6a7(const_cast<const float *>(_tmp0), 0, const_cast<const float **>(godunovMatrix), ExtraOffset_godunovMatrix, godunovState, ExtraOffset_godunovState, NumElements);
    device.api->popStackMemory();
  }
  void kernel::godunovState::execute3() {
    assert(Q != nullptr);
    assert(V3mTo2n(3,0) != nullptr);
    assert(godunovMatrix != nullptr);
    assert(godunovState != nullptr);
    assert(NumElements != 0);
    device::DeviceInstance& device = device::DeviceInstance::getInstance();
    float *_tmp0;
    float* _buffer0 = (float*)device.api->getStackMemory(sizeof(float) * NumElements * 576);
    _tmp0 = _buffer0;

    sgemm_NT_NT_m49_64_n9_64_k56_nps_29706f6(const_cast<const float *>(V3mTo2n(3,0)), 0, const_cast<const float **>(Q), ExtraOffset_Q, _tmp0, 0, NumElements);

    sgemm_NT_NT_m49_64_n9_9_k9_spp_43ea6a7(const_cast<const float *>(_tmp0), 0, const_cast<const float **>(godunovMatrix), ExtraOffset_godunovMatrix, godunovState, ExtraOffset_godunovState, NumElements);
    device.api->popStackMemory();
  }
  void kernel::godunovState::execute4() {
    assert(Q != nullptr);
    assert(V3mTo2n(0,1) != nullptr);
    assert(godunovMatrix != nullptr);
    assert(godunovState != nullptr);
    assert(NumElements != 0);
    device::DeviceInstance& device = device::DeviceInstance::getInstance();
    float *_tmp0;
    float* _buffer0 = (float*)device.api->getStackMemory(sizeof(float) * NumElements * 576);
    _tmp0 = _buffer0;

    sgemm_NT_NT_m49_64_n9_64_k56_nps_29706f6(const_cast<const float *>(V3mTo2n(0,1)), 0, const_cast<const float **>(Q), ExtraOffset_Q, _tmp0, 0, NumElements);

    sgemm_NT_NT_m49_64_n9_9_k9_spp_b4ee274(const_cast<const float *>(_tmp0), 0, const_cast<const float **>(godunovMatrix), ExtraOffset_godunovMatrix, godunovState, ExtraOffset_godunovState, NumElements);
    device.api->popStackMemory();
  }
  void kernel::godunovState::execute5() {
    assert(Q != nullptr);
    assert(V3mTo2n(1,1) != nullptr);
    assert(godunovMatrix != nullptr);
    assert(godunovState != nullptr);
    assert(NumElements != 0);
    device::DeviceInstance& device = device::DeviceInstance::getInstance();
    float *_tmp0;
    float* _buffer0 = (float*)device.api->getStackMemory(sizeof(float) * NumElements * 576);
    _tmp0 = _buffer0;

    sgemm_NT_NT_m49_64_n9_64_k56_nps_29706f6(const_cast<const float *>(V3mTo2n(1,1)), 0, const_cast<const float **>(Q), ExtraOffset_Q, _tmp0, 0, NumElements);

    sgemm_NT_NT_m49_64_n9_9_k9_spp_b4ee274(const_cast<const float *>(_tmp0), 0, const_cast<const float **>(godunovMatrix), ExtraOffset_godunovMatrix, godunovState, ExtraOffset_godunovState, NumElements);
    device.api->popStackMemory();
  }
  void kernel::godunovState::execute6() {
    assert(Q != nullptr);
    assert(V3mTo2n(2,1) != nullptr);
    assert(godunovMatrix != nullptr);
    assert(godunovState != nullptr);
    assert(NumElements != 0);
    device::DeviceInstance& device = device::DeviceInstance::getInstance();
    float *_tmp0;
    float* _buffer0 = (float*)device.api->getStackMemory(sizeof(float) * NumElements * 576);
    _tmp0 = _buffer0;

    sgemm_NT_NT_m49_64_n9_64_k56_nps_29706f6(const_cast<const float *>(V3mTo2n(2,1)), 0, const_cast<const float **>(Q), ExtraOffset_Q, _tmp0, 0, NumElements);

    sgemm_NT_NT_m49_64_n9_9_k9_spp_b4ee274(const_cast<const float *>(_tmp0), 0, const_cast<const float **>(godunovMatrix), ExtraOffset_godunovMatrix, godunovState, ExtraOffset_godunovState, NumElements);
    device.api->popStackMemory();
  }
  void kernel::godunovState::execute7() {
    assert(Q != nullptr);
    assert(V3mTo2n(3,1) != nullptr);
    assert(godunovMatrix != nullptr);
    assert(godunovState != nullptr);
    assert(NumElements != 0);
    device::DeviceInstance& device = device::DeviceInstance::getInstance();
    float *_tmp0;
    float* _buffer0 = (float*)device.api->getStackMemory(sizeof(float) * NumElements * 576);
    _tmp0 = _buffer0;

    sgemm_NT_NT_m49_64_n9_64_k56_nps_29706f6(const_cast<const float *>(V3mTo2n(3,1)), 0, const_cast<const float **>(Q), ExtraOffset_Q, _tmp0, 0, NumElements);

    sgemm_NT_NT_m49_64_n9_9_k9_spp_b4ee274(const_cast<const float *>(_tmp0), 0, const_cast<const float **>(godunovMatrix), ExtraOffset_godunovMatrix, godunovState, ExtraOffset_godunovState, NumElements);
    device.api->popStackMemory();
  }
  void kernel::godunovState::execute8() {
    assert(Q != nullptr);
    assert(V3mTo2n(0,2) != nullptr);
    assert(godunovMatrix != nullptr);
    assert(godunovState != nullptr);
    assert(NumElements != 0);
    device::DeviceInstance& device = device::DeviceInstance::getInstance();
    float *_tmp0;
    float* _buffer0 = (float*)device.api->getStackMemory(sizeof(float) * NumElements * 576);
    _tmp0 = _buffer0;

    sgemm_NT_NT_m49_64_n9_64_k56_nps_29706f6(const_cast<const float *>(V3mTo2n(0,2)), 0, const_cast<const float **>(Q), ExtraOffset_Q, _tmp0, 0, NumElements);

    sgemm_NT_NT_m49_64_n9_9_k9_spp_b4ee274(const_cast<const float *>(_tmp0), 0, const_cast<const float **>(godunovMatrix), ExtraOffset_godunovMatrix, godunovState, ExtraOffset_godunovState, NumElements);
    device.api->popStackMemory();
  }
  void kernel::godunovState::execute9() {
    assert(Q != nullptr);
    assert(V3mTo2n(1,2) != nullptr);
    assert(godunovMatrix != nullptr);
    assert(godunovState != nullptr);
    assert(NumElements != 0);
    device::DeviceInstance& device = device::DeviceInstance::getInstance();
    float *_tmp0;
    float* _buffer0 = (float*)device.api->getStackMemory(sizeof(float) * NumElements * 576);
    _tmp0 = _buffer0;

    sgemm_NT_NT_m49_64_n9_64_k56_nps_29706f6(const_cast<const float *>(V3mTo2n(1,2)), 0, const_cast<const float **>(Q), ExtraOffset_Q, _tmp0, 0, NumElements);

    sgemm_NT_NT_m49_64_n9_9_k9_spp_b4ee274(const_cast<const float *>(_tmp0), 0, const_cast<const float **>(godunovMatrix), ExtraOffset_godunovMatrix, godunovState, ExtraOffset_godunovState, NumElements);
    device.api->popStackMemory();
  }
  void kernel::godunovState::execute10() {
    assert(Q != nullptr);
    assert(V3mTo2n(2,2) != nullptr);
    assert(godunovMatrix != nullptr);
    assert(godunovState != nullptr);
    assert(NumElements != 0);
    device::DeviceInstance& device = device::DeviceInstance::getInstance();
    float *_tmp0;
    float* _buffer0 = (float*)device.api->getStackMemory(sizeof(float) * NumElements * 576);
    _tmp0 = _buffer0;

    sgemm_NT_NT_m49_64_n9_64_k56_nps_29706f6(const_cast<const float *>(V3mTo2n(2,2)), 0, const_cast<const float **>(Q), ExtraOffset_Q, _tmp0, 0, NumElements);

    sgemm_NT_NT_m49_64_n9_9_k9_spp_b4ee274(const_cast<const float *>(_tmp0), 0, const_cast<const float **>(godunovMatrix), ExtraOffset_godunovMatrix, godunovState, ExtraOffset_godunovState, NumElements);
    device.api->popStackMemory();
  }
  void kernel::godunovState::execute11() {
    assert(Q != nullptr);
    assert(V3mTo2n(3,2) != nullptr);
    assert(godunovMatrix != nullptr);
    assert(godunovState != nullptr);
    assert(NumElements != 0);
    device::DeviceInstance& device = device::DeviceInstance::getInstance();
    float *_tmp0;
    float* _buffer0 = (float*)device.api->getStackMemory(sizeof(float) * NumElements * 576);
    _tmp0 = _buffer0;

    sgemm_NT_NT_m49_64_n9_64_k56_nps_29706f6(const_cast<const float *>(V3mTo2n(3,2)), 0, const_cast<const float **>(Q), ExtraOffset_Q, _tmp0, 0, NumElements);

    sgemm_NT_NT_m49_64_n9_9_k9_spp_b4ee274(const_cast<const float *>(_tmp0), 0, const_cast<const float **>(godunovMatrix), ExtraOffset_godunovMatrix, godunovState, ExtraOffset_godunovState, NumElements);
    device.api->popStackMemory();
  }
  void kernel::godunovState::execute12() {
    assert(Q != nullptr);
    assert(V3mTo2n(0,3) != nullptr);
    assert(godunovMatrix != nullptr);
    assert(godunovState != nullptr);
    assert(NumElements != 0);
    device::DeviceInstance& device = device::DeviceInstance::getInstance();
    float *_tmp0;
    float* _buffer0 = (float*)device.api->getStackMemory(sizeof(float) * NumElements * 576);
    _tmp0 = _buffer0;

    sgemm_NT_NT_m49_64_n9_64_k56_nps_29706f6(const_cast<const float *>(V3mTo2n(0,3)), 0, const_cast<const float **>(Q), ExtraOffset_Q, _tmp0, 0, NumElements);

    sgemm_NT_NT_m49_64_n9_9_k9_spp_b4ee274(const_cast<const float *>(_tmp0), 0, const_cast<const float **>(godunovMatrix), ExtraOffset_godunovMatrix, godunovState, ExtraOffset_godunovState, NumElements);
    device.api->popStackMemory();
  }
  void kernel::godunovState::execute13() {
    assert(Q != nullptr);
    assert(V3mTo2n(1,3) != nullptr);
    assert(godunovMatrix != nullptr);
    assert(godunovState != nullptr);
    assert(NumElements != 0);
    device::DeviceInstance& device = device::DeviceInstance::getInstance();
    float *_tmp0;
    float* _buffer0 = (float*)device.api->getStackMemory(sizeof(float) * NumElements * 576);
    _tmp0 = _buffer0;

    sgemm_NT_NT_m49_64_n9_64_k56_nps_29706f6(const_cast<const float *>(V3mTo2n(1,3)), 0, const_cast<const float **>(Q), ExtraOffset_Q, _tmp0, 0, NumElements);

    sgemm_NT_NT_m49_64_n9_9_k9_spp_b4ee274(const_cast<const float *>(_tmp0), 0, const_cast<const float **>(godunovMatrix), ExtraOffset_godunovMatrix, godunovState, ExtraOffset_godunovState, NumElements);
    device.api->popStackMemory();
  }
  void kernel::godunovState::execute14() {
    assert(Q != nullptr);
    assert(V3mTo2n(2,3) != nullptr);
    assert(godunovMatrix != nullptr);
    assert(godunovState != nullptr);
    assert(NumElements != 0);
    device::DeviceInstance& device = device::DeviceInstance::getInstance();
    float *_tmp0;
    float* _buffer0 = (float*)device.api->getStackMemory(sizeof(float) * NumElements * 576);
    _tmp0 = _buffer0;

    sgemm_NT_NT_m49_64_n9_64_k56_nps_29706f6(const_cast<const float *>(V3mTo2n(2,3)), 0, const_cast<const float **>(Q), ExtraOffset_Q, _tmp0, 0, NumElements);

    sgemm_NT_NT_m49_64_n9_9_k9_spp_b4ee274(const_cast<const float *>(_tmp0), 0, const_cast<const float **>(godunovMatrix), ExtraOffset_godunovMatrix, godunovState, ExtraOffset_godunovState, NumElements);
    device.api->popStackMemory();
  }
  void kernel::godunovState::execute15() {
    assert(Q != nullptr);
    assert(V3mTo2n(3,3) != nullptr);
    assert(godunovMatrix != nullptr);
    assert(godunovState != nullptr);
    assert(NumElements != 0);
    device::DeviceInstance& device = device::DeviceInstance::getInstance();
    float *_tmp0;
    float* _buffer0 = (float*)device.api->getStackMemory(sizeof(float) * NumElements * 576);
    _tmp0 = _buffer0;

    sgemm_NT_NT_m49_64_n9_64_k56_nps_29706f6(const_cast<const float *>(V3mTo2n(3,3)), 0, const_cast<const float **>(Q), ExtraOffset_Q, _tmp0, 0, NumElements);

    sgemm_NT_NT_m49_64_n9_9_k9_spp_b4ee274(const_cast<const float *>(_tmp0), 0, const_cast<const float **>(godunovMatrix), ExtraOffset_godunovMatrix, godunovState, ExtraOffset_godunovState, NumElements);
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
    assert(NumElements != 0);
    device::DeviceInstance& device = device::DeviceInstance::getInstance();
    float *_tmp0;
    float* _buffer0 = (float*)device.api->getStackMemory(sizeof(float) * NumElements * 576);
    _tmp0 = _buffer0;

    sgemm_NT_NT_m49_64_n9_9_k9_pps_4046c83(const_cast<const float **>(godunovState), ExtraOffset_godunovState, const_cast<const float **>(fluxSolver), ExtraOffset_fluxSolver, _tmp0, 0, NumElements);

    sgemm_NT_NT_m56_64_n9_64_k49_nsp_01c641c(const_cast<const float *>(V3mTo2nTWDivM(0,0)), 0, const_cast<const float *>(_tmp0), 0, Q, ExtraOffset_Q, NumElements);
    device.api->popStackMemory();
  }
  void kernel::nodalFlux::execute1() {
    assert(Q != nullptr);
    assert(V3mTo2nTWDivM(1,0) != nullptr);
    assert(fluxSolver != nullptr);
    assert(godunovState != nullptr);
    assert(NumElements != 0);
    device::DeviceInstance& device = device::DeviceInstance::getInstance();
    float *_tmp0;
    float* _buffer0 = (float*)device.api->getStackMemory(sizeof(float) * NumElements * 576);
    _tmp0 = _buffer0;

    sgemm_NT_NT_m49_64_n9_9_k9_pps_4046c83(const_cast<const float **>(godunovState), ExtraOffset_godunovState, const_cast<const float **>(fluxSolver), ExtraOffset_fluxSolver, _tmp0, 0, NumElements);

    sgemm_NT_NT_m56_64_n9_64_k49_nsp_01c641c(const_cast<const float *>(V3mTo2nTWDivM(1,0)), 0, const_cast<const float *>(_tmp0), 0, Q, ExtraOffset_Q, NumElements);
    device.api->popStackMemory();
  }
  void kernel::nodalFlux::execute2() {
    assert(Q != nullptr);
    assert(V3mTo2nTWDivM(2,0) != nullptr);
    assert(fluxSolver != nullptr);
    assert(godunovState != nullptr);
    assert(NumElements != 0);
    device::DeviceInstance& device = device::DeviceInstance::getInstance();
    float *_tmp0;
    float* _buffer0 = (float*)device.api->getStackMemory(sizeof(float) * NumElements * 576);
    _tmp0 = _buffer0;

    sgemm_NT_NT_m49_64_n9_9_k9_pps_4046c83(const_cast<const float **>(godunovState), ExtraOffset_godunovState, const_cast<const float **>(fluxSolver), ExtraOffset_fluxSolver, _tmp0, 0, NumElements);

    sgemm_NT_NT_m56_64_n9_64_k49_nsp_01c641c(const_cast<const float *>(V3mTo2nTWDivM(2,0)), 0, const_cast<const float *>(_tmp0), 0, Q, ExtraOffset_Q, NumElements);
    device.api->popStackMemory();
  }
  void kernel::nodalFlux::execute3() {
    assert(Q != nullptr);
    assert(V3mTo2nTWDivM(3,0) != nullptr);
    assert(fluxSolver != nullptr);
    assert(godunovState != nullptr);
    assert(NumElements != 0);
    device::DeviceInstance& device = device::DeviceInstance::getInstance();
    float *_tmp0;
    float* _buffer0 = (float*)device.api->getStackMemory(sizeof(float) * NumElements * 576);
    _tmp0 = _buffer0;

    sgemm_NT_NT_m49_64_n9_9_k9_pps_4046c83(const_cast<const float **>(godunovState), ExtraOffset_godunovState, const_cast<const float **>(fluxSolver), ExtraOffset_fluxSolver, _tmp0, 0, NumElements);

    sgemm_NT_NT_m56_64_n9_64_k49_nsp_01c641c(const_cast<const float *>(V3mTo2nTWDivM(3,0)), 0, const_cast<const float *>(_tmp0), 0, Q, ExtraOffset_Q, NumElements);
    device.api->popStackMemory();
  }
  void kernel::nodalFlux::execute4() {
    assert(Q != nullptr);
    assert(V3mTo2nTWDivM(0,1) != nullptr);
    assert(fluxSolver != nullptr);
    assert(godunovState != nullptr);
    assert(NumElements != 0);
    device::DeviceInstance& device = device::DeviceInstance::getInstance();
    float *_tmp0;
    float* _buffer0 = (float*)device.api->getStackMemory(sizeof(float) * NumElements * 576);
    _tmp0 = _buffer0;

    sgemm_NT_NT_m49_64_n9_9_k9_pps_4046c83(const_cast<const float **>(godunovState), ExtraOffset_godunovState, const_cast<const float **>(fluxSolver), ExtraOffset_fluxSolver, _tmp0, 0, NumElements);

    sgemm_NT_NT_m56_64_n9_64_k49_nsp_01c641c(const_cast<const float *>(V3mTo2nTWDivM(0,1)), 0, const_cast<const float *>(_tmp0), 0, Q, ExtraOffset_Q, NumElements);
    device.api->popStackMemory();
  }
  void kernel::nodalFlux::execute5() {
    assert(Q != nullptr);
    assert(V3mTo2nTWDivM(1,1) != nullptr);
    assert(fluxSolver != nullptr);
    assert(godunovState != nullptr);
    assert(NumElements != 0);
    device::DeviceInstance& device = device::DeviceInstance::getInstance();
    float *_tmp0;
    float* _buffer0 = (float*)device.api->getStackMemory(sizeof(float) * NumElements * 576);
    _tmp0 = _buffer0;

    sgemm_NT_NT_m49_64_n9_9_k9_pps_4046c83(const_cast<const float **>(godunovState), ExtraOffset_godunovState, const_cast<const float **>(fluxSolver), ExtraOffset_fluxSolver, _tmp0, 0, NumElements);

    sgemm_NT_NT_m56_64_n9_64_k49_nsp_01c641c(const_cast<const float *>(V3mTo2nTWDivM(1,1)), 0, const_cast<const float *>(_tmp0), 0, Q, ExtraOffset_Q, NumElements);
    device.api->popStackMemory();
  }
  void kernel::nodalFlux::execute6() {
    assert(Q != nullptr);
    assert(V3mTo2nTWDivM(2,1) != nullptr);
    assert(fluxSolver != nullptr);
    assert(godunovState != nullptr);
    assert(NumElements != 0);
    device::DeviceInstance& device = device::DeviceInstance::getInstance();
    float *_tmp0;
    float* _buffer0 = (float*)device.api->getStackMemory(sizeof(float) * NumElements * 576);
    _tmp0 = _buffer0;

    sgemm_NT_NT_m49_64_n9_9_k9_pps_4046c83(const_cast<const float **>(godunovState), ExtraOffset_godunovState, const_cast<const float **>(fluxSolver), ExtraOffset_fluxSolver, _tmp0, 0, NumElements);

    sgemm_NT_NT_m56_64_n9_64_k49_nsp_01c641c(const_cast<const float *>(V3mTo2nTWDivM(2,1)), 0, const_cast<const float *>(_tmp0), 0, Q, ExtraOffset_Q, NumElements);
    device.api->popStackMemory();
  }
  void kernel::nodalFlux::execute7() {
    assert(Q != nullptr);
    assert(V3mTo2nTWDivM(3,1) != nullptr);
    assert(fluxSolver != nullptr);
    assert(godunovState != nullptr);
    assert(NumElements != 0);
    device::DeviceInstance& device = device::DeviceInstance::getInstance();
    float *_tmp0;
    float* _buffer0 = (float*)device.api->getStackMemory(sizeof(float) * NumElements * 576);
    _tmp0 = _buffer0;

    sgemm_NT_NT_m49_64_n9_9_k9_pps_4046c83(const_cast<const float **>(godunovState), ExtraOffset_godunovState, const_cast<const float **>(fluxSolver), ExtraOffset_fluxSolver, _tmp0, 0, NumElements);

    sgemm_NT_NT_m56_64_n9_64_k49_nsp_01c641c(const_cast<const float *>(V3mTo2nTWDivM(3,1)), 0, const_cast<const float *>(_tmp0), 0, Q, ExtraOffset_Q, NumElements);
    device.api->popStackMemory();
  }
  void kernel::nodalFlux::execute8() {
    assert(Q != nullptr);
    assert(V3mTo2nTWDivM(0,2) != nullptr);
    assert(fluxSolver != nullptr);
    assert(godunovState != nullptr);
    assert(NumElements != 0);
    device::DeviceInstance& device = device::DeviceInstance::getInstance();
    float *_tmp0;
    float* _buffer0 = (float*)device.api->getStackMemory(sizeof(float) * NumElements * 576);
    _tmp0 = _buffer0;

    sgemm_NT_NT_m49_64_n9_9_k9_pps_4046c83(const_cast<const float **>(godunovState), ExtraOffset_godunovState, const_cast<const float **>(fluxSolver), ExtraOffset_fluxSolver, _tmp0, 0, NumElements);

    sgemm_NT_NT_m56_64_n9_64_k49_nsp_01c641c(const_cast<const float *>(V3mTo2nTWDivM(0,2)), 0, const_cast<const float *>(_tmp0), 0, Q, ExtraOffset_Q, NumElements);
    device.api->popStackMemory();
  }
  void kernel::nodalFlux::execute9() {
    assert(Q != nullptr);
    assert(V3mTo2nTWDivM(1,2) != nullptr);
    assert(fluxSolver != nullptr);
    assert(godunovState != nullptr);
    assert(NumElements != 0);
    device::DeviceInstance& device = device::DeviceInstance::getInstance();
    float *_tmp0;
    float* _buffer0 = (float*)device.api->getStackMemory(sizeof(float) * NumElements * 576);
    _tmp0 = _buffer0;

    sgemm_NT_NT_m49_64_n9_9_k9_pps_4046c83(const_cast<const float **>(godunovState), ExtraOffset_godunovState, const_cast<const float **>(fluxSolver), ExtraOffset_fluxSolver, _tmp0, 0, NumElements);

    sgemm_NT_NT_m56_64_n9_64_k49_nsp_01c641c(const_cast<const float *>(V3mTo2nTWDivM(1,2)), 0, const_cast<const float *>(_tmp0), 0, Q, ExtraOffset_Q, NumElements);
    device.api->popStackMemory();
  }
  void kernel::nodalFlux::execute10() {
    assert(Q != nullptr);
    assert(V3mTo2nTWDivM(2,2) != nullptr);
    assert(fluxSolver != nullptr);
    assert(godunovState != nullptr);
    assert(NumElements != 0);
    device::DeviceInstance& device = device::DeviceInstance::getInstance();
    float *_tmp0;
    float* _buffer0 = (float*)device.api->getStackMemory(sizeof(float) * NumElements * 576);
    _tmp0 = _buffer0;

    sgemm_NT_NT_m49_64_n9_9_k9_pps_4046c83(const_cast<const float **>(godunovState), ExtraOffset_godunovState, const_cast<const float **>(fluxSolver), ExtraOffset_fluxSolver, _tmp0, 0, NumElements);

    sgemm_NT_NT_m56_64_n9_64_k49_nsp_01c641c(const_cast<const float *>(V3mTo2nTWDivM(2,2)), 0, const_cast<const float *>(_tmp0), 0, Q, ExtraOffset_Q, NumElements);
    device.api->popStackMemory();
  }
  void kernel::nodalFlux::execute11() {
    assert(Q != nullptr);
    assert(V3mTo2nTWDivM(3,2) != nullptr);
    assert(fluxSolver != nullptr);
    assert(godunovState != nullptr);
    assert(NumElements != 0);
    device::DeviceInstance& device = device::DeviceInstance::getInstance();
    float *_tmp0;
    float* _buffer0 = (float*)device.api->getStackMemory(sizeof(float) * NumElements * 576);
    _tmp0 = _buffer0;

    sgemm_NT_NT_m49_64_n9_9_k9_pps_4046c83(const_cast<const float **>(godunovState), ExtraOffset_godunovState, const_cast<const float **>(fluxSolver), ExtraOffset_fluxSolver, _tmp0, 0, NumElements);

    sgemm_NT_NT_m56_64_n9_64_k49_nsp_01c641c(const_cast<const float *>(V3mTo2nTWDivM(3,2)), 0, const_cast<const float *>(_tmp0), 0, Q, ExtraOffset_Q, NumElements);
    device.api->popStackMemory();
  }
  void kernel::nodalFlux::execute12() {
    assert(Q != nullptr);
    assert(V3mTo2nTWDivM(0,3) != nullptr);
    assert(fluxSolver != nullptr);
    assert(godunovState != nullptr);
    assert(NumElements != 0);
    device::DeviceInstance& device = device::DeviceInstance::getInstance();
    float *_tmp0;
    float* _buffer0 = (float*)device.api->getStackMemory(sizeof(float) * NumElements * 576);
    _tmp0 = _buffer0;

    sgemm_NT_NT_m49_64_n9_9_k9_pps_4046c83(const_cast<const float **>(godunovState), ExtraOffset_godunovState, const_cast<const float **>(fluxSolver), ExtraOffset_fluxSolver, _tmp0, 0, NumElements);

    sgemm_NT_NT_m56_64_n9_64_k49_nsp_01c641c(const_cast<const float *>(V3mTo2nTWDivM(0,3)), 0, const_cast<const float *>(_tmp0), 0, Q, ExtraOffset_Q, NumElements);
    device.api->popStackMemory();
  }
  void kernel::nodalFlux::execute13() {
    assert(Q != nullptr);
    assert(V3mTo2nTWDivM(1,3) != nullptr);
    assert(fluxSolver != nullptr);
    assert(godunovState != nullptr);
    assert(NumElements != 0);
    device::DeviceInstance& device = device::DeviceInstance::getInstance();
    float *_tmp0;
    float* _buffer0 = (float*)device.api->getStackMemory(sizeof(float) * NumElements * 576);
    _tmp0 = _buffer0;

    sgemm_NT_NT_m49_64_n9_9_k9_pps_4046c83(const_cast<const float **>(godunovState), ExtraOffset_godunovState, const_cast<const float **>(fluxSolver), ExtraOffset_fluxSolver, _tmp0, 0, NumElements);

    sgemm_NT_NT_m56_64_n9_64_k49_nsp_01c641c(const_cast<const float *>(V3mTo2nTWDivM(1,3)), 0, const_cast<const float *>(_tmp0), 0, Q, ExtraOffset_Q, NumElements);
    device.api->popStackMemory();
  }
  void kernel::nodalFlux::execute14() {
    assert(Q != nullptr);
    assert(V3mTo2nTWDivM(2,3) != nullptr);
    assert(fluxSolver != nullptr);
    assert(godunovState != nullptr);
    assert(NumElements != 0);
    device::DeviceInstance& device = device::DeviceInstance::getInstance();
    float *_tmp0;
    float* _buffer0 = (float*)device.api->getStackMemory(sizeof(float) * NumElements * 576);
    _tmp0 = _buffer0;

    sgemm_NT_NT_m49_64_n9_9_k9_pps_4046c83(const_cast<const float **>(godunovState), ExtraOffset_godunovState, const_cast<const float **>(fluxSolver), ExtraOffset_fluxSolver, _tmp0, 0, NumElements);

    sgemm_NT_NT_m56_64_n9_64_k49_nsp_01c641c(const_cast<const float *>(V3mTo2nTWDivM(2,3)), 0, const_cast<const float *>(_tmp0), 0, Q, ExtraOffset_Q, NumElements);
    device.api->popStackMemory();
  }
  void kernel::nodalFlux::execute15() {
    assert(Q != nullptr);
    assert(V3mTo2nTWDivM(3,3) != nullptr);
    assert(fluxSolver != nullptr);
    assert(godunovState != nullptr);
    assert(NumElements != 0);
    device::DeviceInstance& device = device::DeviceInstance::getInstance();
    float *_tmp0;
    float* _buffer0 = (float*)device.api->getStackMemory(sizeof(float) * NumElements * 576);
    _tmp0 = _buffer0;

    sgemm_NT_NT_m49_64_n9_9_k9_pps_4046c83(const_cast<const float **>(godunovState), ExtraOffset_godunovState, const_cast<const float **>(fluxSolver), ExtraOffset_fluxSolver, _tmp0, 0, NumElements);

    sgemm_NT_NT_m56_64_n9_64_k49_nsp_01c641c(const_cast<const float *>(V3mTo2nTWDivM(3,3)), 0, const_cast<const float *>(_tmp0), 0, Q, ExtraOffset_Q, NumElements);
    device.api->popStackMemory();
  }
}
