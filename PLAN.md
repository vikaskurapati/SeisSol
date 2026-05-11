# Triton Backend Integration Plan

**Goal:** Integrate Triton as a GPU kernel generator backend for SeisSol/YATeTo to benchmark against TensorForge.

**Approach:** 
- Follow TDD workflow: write tests first, then implement
- Use tinytc backend as reference pattern
- AOT compilation with LLVM/PTX (no Python runtime dependency)
- Start with basic GEMMs for proxy app, add fused GEMMs after benchmarking

**Key Constraints:**
- No GPU on dev machine (compilation checks only, performance on cluster)
- Avoid modifying ChainForge/GemmForge/TensorForge
- Write tests for all changes (YATeTo and SeisSol)

---

## Progress Summary

### ✅ Phase 1: YATeTo Triton GemmTool Class (COMPLETED)
**Status:** Done - Commit fc01dec

**What was done:**
- Created comprehensive test suite: `submodules/yateto/tests/backends/test_triton_backend.py` (23 tests)
- Implemented `Triton` class in `submodules/yateto/yateto/gemm_configuration.py`
- All tests passing (23/23)

**Key decisions:**
- Triton supports CUDA and HIP backends only (not SYCL/oneAPI)
- Flexible alpha/beta support (any values, unlike tinytc's restrictions)
- No sparse matrix support initially
- HIGHEST preference for operations it supports
- Follows CodeGenerator base class pattern

---

### ✅ Phase 2: YATeTo Common Infrastructure (COMPLETED)
**Status:** Done - Commit [hash]

**What was done:**
- Created test suite: `submodules/yateto/tests/codegen/test_triton_common.py` (18 tests)
- Implemented `submodules/yateto/yateto/codegen/triton_common.py` with:
  - TritonKernelArgument / TritonScalarKernelArgument classes
  - TritonWrapper class for C++ wrapper generation
  - compile_triton_kernel() for AOT compilation
  - make_triton_kernel_name() for kernel naming
- All tests passing (18/18)

**Key implementation details:**
- Wrappers use CUDA Driver API to load pre-compiled .so/.cubin files
- No source code embedded in C++ (unlike tinytc's runtime JIT)
- Supports arbitrary alpha/beta values
- Hash-based wrapper names for uniqueness

---

## Task Breakdown (14 total)

### YATeTo Infrastructure (Phase 1)
- [x] `yateto-triton-class` - Triton GemmTool class ✅ DONE
- [ ] `yateto-triton-common` - Common infrastructure 🔄 IN PROGRESS
- [ ] `yateto-triton-gemm` - GEMM generator (depends: triton-common)
- [ ] `yateto-triton-cache` - TritonWriter cache (depends: triton-common)
- [ ] `yateto-triton-factory` - Factory updates (depends: triton-cache)

### YATeTo Testing (Phase 2)
- [ ] `yateto-tests-triton` - Code generation tests (depends: triton-gemm, triton-cache)

### SeisSol Integration (Phase 3)
- [ ] `seissol-cmake-triton` - CMake integration (depends: yateto-tests-triton)
- [ ] `seissol-codegen-triton` - Code generator updates (depends: seissol-cmake-triton)
- [ ] `seissol-proxy-test` - Proxy app builds (depends: seissol-codegen-triton)

### SeisSol Testing (Phase 4)
- [ ] `seissol-tests-triton` - Integration tests (depends: seissol-proxy-test)

### Benchmarking (Phase 5)
- [ ] `benchmark-setup` - Create configs (depends: seissol-tests-triton)
- [ ] `benchmark-run` - Run on GPU cluster (depends: benchmark-setup) ⚠️ Requires Vista
- [ ] `benchmark-analysis` - Analyze results (depends: benchmark-run)

### Future Work (Phase 6)
- [ ] `fused-gemm-explore` - Fused GEMMs (depends: benchmark-analysis)

---

## Technical Details

### Architecture Comparison: tinytc vs Triton

**tinytc (existing):**
- Runtime JIT compilation via SYCL
- Compiles kernels at runtime using tinytc::parse_string()
- Wrapper contains source code as R"tinytc(...)tinytc" string literal
- Creates SYCL kernel bundle on-demand
- Restricts alpha=1.0, beta∈{0.0,1.0}

**Triton (new):**
- Build-time AOT compilation
- Compiles during code generation using triton.compile()
- Wrapper links to pre-compiled .cubin or .so files
- No Python runtime dependency
- Supports any alpha/beta values

### YATeTo Integration Pattern

**Files modified so far:**
1. `submodules/yateto/yateto/gemm_configuration.py` - Added Triton class (lines 277-331)
2. `submodules/yateto/tests/backends/test_triton_backend.py` - 23 unit tests

**Files to create:**
1. `submodules/yateto/yateto/codegen/triton_common.py` - Infrastructure
2. `submodules/yateto/yateto/codegen/gemms/triton.py` - GEMM generator
3. `submodules/yateto/yateto/codegen/fused_gemms/triton.py` - Fused GEMM generator (later)

**Files to modify:**
1. `submodules/yateto/yateto/codegen/cache.py` - Add TritonWriter
2. `submodules/yateto/yateto/gemm_configuration.py` - Update GeneratorCollection if needed

### Matrix Sizes from SeisSol (DG order 2-7 elastic)
- Volume/time: (10×9×10) to (120×9×120)
- Flux/face: (56×21×21), (84×28×28)
- ADER: (56×56×56), (84×84×84)

### Testing Framework
- YATeTo uses Python unittest (not pytest)
- Run tests: `cd submodules/yateto && python3 tests/backends/test_triton_backend.py -v`
- All changes must have corresponding tests

---

## Important Files Reference

### YATeTo
- `yateto/gemm_configuration.py` - Backend registry (Triton class: lines 277-331)
- `yateto/codegen/common.py` - TinytcKernelArgument/Wrapper classes (reference)
- `yateto/codegen/cache.py` - Code generation cache system
- `yateto/codegen/gemms/tinytc.py` - tinytc GEMM generator (reference)
- `yateto/codegen/fused_gemms/tinytc.py` - tinytc fused GEMM generator (reference)

### SeisSol
- To be determined during integration phase

---

## Git Branches
- YATeTo: `vikas/triton` (currently on this branch)
- SeisSol: main (will create branch when needed)

---

## Notes for Recovery
If this session crashes, next agent should:
1. Check git status and latest commit
2. Query SQL database for current todo status: `SELECT * FROM todos ORDER BY status, id;`
3. Read session files in `~/.copilot/session-state/.../files/`
4. Continue with current `in_progress` task or next `pending` task with no deps

**Current context:**
- Phase 3 (GEMM & Cache Integration) is complete.
- Working on Phase 4: `yateto-triton-csa` - generating CopyScaleAdd operations in Triton to fully replace TensorForge.

---

### ✅ Phase 3: YATeTo GEMM Generator & Testing (DONE)
**What's done so far:**
- ✅ Created `triton_common.py` infrastructure (Phases 1-2)
- ✅ Added `TritonWriter` class to `cache.py` for kernel compilation and linking
- ✅ Created `codegen/gemm/triton.py` with `tritonGemmGen()` function
- ✅ Created GEMM integration tests in `test_triton_gemm.py`
- ✅ Integrated into `GemmGen.generate()` method to handle Triton backend
- ✅ Tested cache integration and fallback mechanisms
- ✅ Added Triton to `DefaultGeneratorCollection` in `gemm_configuration.py`
- ✅ CMake Integration in SeisSol (`FindGemmTools.cmake` & `process_users_input.cmake`)

---

### 🔄 Phase 4: YATeTo CopyScaleAdd Generator + Pure Triton Codegen (IN PROGRESS)
**Context:** SeisSol requires auxiliary routines (Copy, Scale, Add) for the proxy app. Triton was initially only generating GEMMs, which caused compilation to fail when TensorForge was removed. We are implementing pure Triton code generation for both GEMM and CSA paths to achieve full autonomy from TensorForge/GemmForge during GPU code generation.

**What's done so far:**
- ✅ Created `yateto/codegen/copyscaleadd/triton.py` with `CopyScaleAddTriton` generator.
- ✅ Updated `yateto/codegen/copyscaleadd/factory.py` to route GPU CSA requests to Triton when active.
- ✅ Fixed `TritonWrapper` and `TritonWriter` argument mismatch issues.
- ✅ Improved Triton AOT kernel discovery logic in `triton_common.py` to detect modern Triton `@triton.jit` functions more robustly.
- ✅ Confirmed the active SeisSol codegen path uses `codegen/yateto/...`, not only `submodules/yateto/...`.
- ✅ Fixed Python import compatibility in `copyscaleadd/factory.py` (`import importlib.util`).
- ✅ Added Triton compile API compatibility in `triton_common.py`:
  - Tries `kernel_fn.compile(...)` when available.
  - Falls back to `triton.compile(...)` and `triton.compiler.compile(...)`.
  - Handles multiple compiled artifact shapes (`asm` dict, direct attrs, returned file path).
- ✅ Added Triton target compatibility for newer APIs that require `GPUTarget` objects:
  - Detects available `GPUTarget` classes from Triton backend modules.
  - Builds target candidates from backend/arch variants (e.g. CUDA `sm_90`/`90`/`sm90`).
  - Retries compilation across target candidates for `kernel_fn.compile`, `triton.compile`, and `triton.compiler.compile`.
- ✅ Added regression tests:
  - `tests/codegen/test_triton_gemm.py::test_gemm_gen_custom_kernel_name`
  - `tests/codegen/test_triton_common.py::test_compile_kernel_contains_api_compat_fallbacks` (now also checks `GPUTarget` support)

**Current Blockers:**
- ⚠️ Cluster validation still pending for full SeisSol proxy build with Triton-only device codegen.

**Next steps:**
1. Re-run YATeTo Triton tests on cluster Python.
2. Re-run SeisSol code generation + build for `DEVICE_CODEGEN=triton`.
3. Verify `seissol-proxy` runs, then proceed to benchmarking.

---

**Last updated:** 2026-05-11 (Added Triton `GPUTarget` compatibility for clusters where string targets are rejected)
