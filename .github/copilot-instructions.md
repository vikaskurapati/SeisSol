# SeisSol Development Instructions

This file helps AI assistants (GitHub Copilot, etc.) work effectively in the SeisSol codebase.

## Project Overview

SeisSol is a scientific software for numerical simulation of seismic wave phenomena and earthquake dynamics. It uses the discontinuous Galerkin method with ADER time discretization and supports both CPU and GPU execution.

**Core Technologies:** C++17, Python 3.9+, CMake, MPI, HDF5, OpenMP

## Build Commands

### Basic Build

```bash
# Configure with typical options (from SeisSol root)
mkdir -p build && cd build
cmake -DORDER=4 -DEQUATIONS=elastic -DPRECISION=double ..
make -j $(nproc)
```

### Key CMake Options

**Required for any build:**
- `-DORDER=<2|3|4|5|6|7|8>` - Polynomial convergence order (commonly 4 or 6)
- `-DEQUATIONS=<elastic|anisotropic|viscoelastic|viscoelastic2|poroelastic|acoustic>` - Equation set
- `-DPRECISION=<single|double>` - Floating point precision

**GPU builds:**
- `-DDEVICE_BACKEND=<cuda|hip|acpp|oneapi>` - Enable GPU support
- `-DDEVICE_ARCH=<sm_90|gfx90a|...>` - GPU architecture (see `cmake/process_users_input.cmake` for full list)

**Common development flags:**
- `-DTESTING=ON` - Enable test suite
- `-DCMAKE_BUILD_TYPE=<Release|Debug|RelWithDebInfo>` - Build type (default: Release)
- `-DHOST_ARCH=<auto|hsw|skx|...>` - CPU architecture (auto-detects by default)
- `-DGEMM_TOOLS_LIST=<libxsmm,pspamm,etc>` - Matrix multiplication libraries

**Optional features:**
- `-DASAGI=ON` - Enable ASAGI for material input
- `-DGRAPH_PARTITIONING_LIBS=<parmetis|parhip|ptscotch>` - Mesh partitioning library

Example GPU build for NVIDIA Hopper:
```bash
cmake -DORDER=4 -DPRECISION=single \
      -DDEVICE_BACKEND=cuda -DDEVICE_ARCH=sm_90 ..
```

### Important Build Notes

- **Always use `--recursive`** when cloning: `git clone --recursive https://github.com/SeisSol/SeisSol.git`
- **Update submodules** after pulling: `git submodule update --recursive --init`
- Code generation happens at build time via Python scripts in `codegen/`
- For clean rebuilds affecting generated code, remove the build directory entirely

## Test Commands

### Running Tests

```bash
# Enable and run full test suite
cd build
cmake .. -DTESTING=ON
make -j $(nproc)
ctest --output-on-failure
```

### Running Specific Tests

```bash
# Run tests matching a pattern
ctest -R TestModel --output-on-failure

# Run a specific test with verbose output
ctest -R TestGeometry -V

# List all available tests
ctest -N
```

### Test for Generated Code

```bash
# Validate generated kernels (CI "sanity check")
cmake .. -DTESTING=ON -DTESTING_GENERATED=ON
make -j $(nproc)
ctest --output-on-failure
```

If "yateto" tests fail, there may be an issue in the code generation system (`codegen/`).

## Lint and Format Commands

### Pre-commit Hooks

The project uses pre-commit for automated checks. **Install before making changes:**

```bash
pip install pre-commit
pre-commit install
```

This runs checks automatically on commit (clang-format, flake8, bandit, isort, black, etc.).

### Manual Formatting

**C++:**
```bash
# Format a single file
clang-format -i src/path/to/file.cpp

# Format modified files only (via pre-commit)
pre-commit run clang-format --files src/MyFile.cpp
```

**Python:**
```bash
# Format with black and isort (via pre-commit)
pre-commit run black --all-files
pre-commit run isort --all-files
```

### Static Analysis

**clang-tidy** (required for CI):
```bash
# Enable compile commands export
cmake .. -DCMAKE_EXPORT_COMPILE_COMMANDS=ON

# Run clang-tidy from repo root
.ci/tidy.sh ./ build/

# Run with auto-fix (use with caution)
.ci/tidy.sh ./ build/ -fix
```

Configuration: `.clang-tidy` (contains rules and exceptions)

**Python linting:**
```bash
# Run flake8 (config in .flake8)
python -m flake8

# Run bandit security checks
bandit -c .bandit -r codegen/
```

## Architecture Overview

### Directory Structure

```
SeisSol/
├── src/              # C++ core simulation code
│   ├── Geometry/         # Mesh handling and refinement
│   ├── Kernels/          # Computation kernels (CPU/GPU)
│   ├── Solver/           # Time stepping and solver logic
│   ├── DynamicRupture/   # Fault dynamics and friction laws
│   ├── Initializer/      # Setup and initialization
│   ├── IO/               # Input/output operations
│   ├── Equations/        # Equation-specific implementations
│   ├── Physics/          # Physical models
│   ├── Parallel/         # MPI and communication
│   └── SeisSol.cpp       # Main application entry
├── app/              # Main executables and proxy applications
│   ├── Main/             # Primary SeisSol executable
│   └── Proxy/            # Performance testing proxy
├── codegen/          # Python-based code generation system
│   ├── generate.py       # Main code generator
│   ├── kernels/          # Kernel generation logic
│   └── yateto/           # Matrix operation library (submodule)
├── tests/            # Unit tests (doctest framework)
├── cmake/            # CMake configuration and build system
├── preprocessing/    # Python tools for input preparation
├── postprocessing/   # Python tools for output analysis
├── docs/             # Sphinx documentation
└── submodules/       # External dependencies (git submodules)
```

### Key Architectural Concepts

**Code Generation:**
- SeisSol generates optimized matrix kernels **at build time** using Python (`codegen/generate.py`)
- Generated code depends on `ORDER`, `EQUATIONS`, `PRECISION`, `HOST_ARCH`, `DEVICE_BACKEND`
- Uses `yateto` library for tensor operations and kernel generation
- Generated files are placed in `build/src/generated_code/`

**Kernel Abstractions:**
- `Kernels/` contains both host (CPU) and device (GPU) implementations
- Separate subdirectories for different equation types (e.g., `LinearCK`, `LinearCKAnelastic`)
- GPU kernels support CUDA, HIP, and SYCL backends via abstraction layer

**Namespace Hierarchy:**
- Base namespace: `seissol::`
- Sub-namespaces follow directory structure: `seissol::geometry::`, `seissol::refinement::`, etc.
- Code consistently uses nested namespaces (C++17 style)

**Equation Sets:**
- Different physics models (elastic, anisotropic, viscoelastic, etc.) selected at **compile time**
- Equation-specific code in `src/Equations/`
- Number of variables, material parameters vary by equation type

**MPI Parallelization:**
- Domain decomposition using graph partitioning (ParMETIS, ParHIP, PT-SCOTCH)
- Asynchronous MPI communication with optional communication thread
- Node-local shared memory optimization

**Time Stepping:**
- ADER-DG method with local time stepping (LTS)
- `TimeCluster` abstraction for different time step groups
- Code in `src/Solver/` and time kernels in `src/Kernels/Time*.{h,cpp}`

## Code Conventions

### C++ Specifics

**Language Standard:** C++17 with `__restrict` keyword extension for array access optimization

**Formatting:** Enforced by `.clang-format` (clang-format 22.1.0):
- Run pre-commit hooks or manually format before committing
- Configuration includes project-specific exceptions

**Logging:**
```cpp
#include <utils/logger.h>

// Use streaming API (from utils submodule)
logInfo() << "Message with" << variable;
logWarning() << "Warning message";
logError() << "Error message";
```

**Testing:** Uses doctest framework
- Test files: `tests/<Module>/Test<Module>.cpp` include `.t.h` header files
- `.t.h` files contain actual test cases
- Example: `tests/Geometry/TestGeometry.cpp` includes `MeshReader.t.h`

**Code Organization:**
- Header guards use `#pragma once`
- Prefer nested namespace declarations: `namespace seissol::geometry { ... }`
- Follow existing patterns in the module you're modifying

**SPDX License Headers:**
All source files must include SPDX headers:
```cpp
// SPDX-FileCopyrightText: 2025 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
```

### Python Specifics

**Style:** PEP 8 enforced by flake8, black, and isort
- Configuration in `.flake8` (lists excepted rules)
- Black formatting with `--profile black` for isort
- Security checks via bandit

**Code Generation Scripts:**
- Located in `codegen/`
- Must pass flake8, black, isort, and bandit
- Critical for build process - test thoroughly

**Pre/Post-processing Tools:**
- In `preprocessing/` and `postprocessing/`
- More lenient formatting requirements (but encouraged to follow)

### Git Workflow

**Branching:**
- Trunk-based development on `master` branch
- Feature branches: `feature/my-feature`, `fix/bug-description`, `docs/update-guide`
- Always branch from up-to-date master with submodules

**Commit Messages:**
Encouraged pattern (not strictly enforced):
```
<type>: <short summary in imperative>

<optional body explaining motivation>

<optional footer: Fixes #123>
```

Types: `fix`, `feat`, `docs`, `refactor`, `test`, `perf`, `ci`, `build`

**Before Submitting PR:**
1. Run pre-commit hooks or manual formatting/linting
2. Build successfully with your configuration
3. Run test suite: `ctest --output-on-failure`
4. Run clang-tidy: `.ci/tidy.sh ./ build/`
5. Update documentation if adding features or parameters

**CI Requirements:**
- All builds must succeed (CPU and GPU configurations tested)
- All tests must pass
- clang-format, clang-tidy, flake8, bandit must pass
- Documentation must build without errors

## Common Development Tasks

### Adding a New CMake Option

1. Add option to `cmake/process_users_input.cmake`
2. Document in `docs/build-parameters.rst`
3. Update `README.md` if user-facing

### Modifying Code Generation

1. Edit scripts in `codegen/kernels/`
2. Clean build directory: `rm -rf build/`
3. Rebuild and verify generated code compiles
4. Run with `-DTESTING_GENERATED=ON` to validate

### Adding Tests

1. Create `Test<Module>.t.h` in `tests/<Module>/`
2. Include in `Test<Module>.cpp` (or create if new module)
3. Use doctest macros: `TEST_CASE`, `CHECK`, `REQUIRE`
4. Add test sources to `tests/CMakeLists.txt` if needed

### Working with Submodules

```bash
# Update a specific submodule
git submodule update --init --recursive submodules/yateto

# Update all submodules to latest
git submodule update --remote --recursive

# Commit submodule changes
git add submodules/yateto
git commit -m "refactor: update yateto submodule"
```

### Debugging Build Issues

1. Check submodules are initialized: `git submodule status`
2. Verify dependencies are installed (see `docs/build-dependencies.html`)
3. Clean build: `rm -rf build/ && mkdir build && cd build`
4. Use `CMAKE_BUILD_TYPE=Debug` for debug symbols
5. Check `CMakeCache.txt` for actual configuration values

## Documentation

**Location:** `docs/` (Sphinx/reStructuredText)

**Build locally:**
```bash
cd docs
pip install -r requirements.txt
make html
# Open _build/html/index.html
```

**When to update:**
- New build parameters → update `docs/parameters.par`
- New features → add to relevant `.rst` files
- Changed behavior → update user guide

**Full docs:** [seissol.readthedocs.io](https://seissol.readthedocs.io)

## Important Environment Variables

- `SEISSOL_COMMTHREAD=<0|1>` - Enable/disable communication thread (default: auto)
- `CMAKE_PREFIX_PATH` - Paths to dependencies (if not in standard locations)

## Additional Resources

- **Contributing Guide:** `CONTRIBUTING.md` - detailed PR process, code style, testing
- **Build Documentation:** `docs/build-overview.html` and `docs/build-parameters.html`
- **GitHub Discussions:** For questions and community support
- **Issues:** For bug reports with build config and reproduction steps

## Notes for AI Assistants

- **Always check submodules are up-to-date** when build errors occur
- **CMake options affect code generation** - changing ORDER/EQUATIONS requires rebuild
- **Multi-configuration project** - test changes with different ORDER/EQUATIONS if core code affected
- **CI is strict** - ensure formatting/linting passes before suggesting changes
- **Performance-critical code** - kernel modifications need careful consideration of generated code
- Pre/postprocessing scripts have **looser style requirements** than core code
