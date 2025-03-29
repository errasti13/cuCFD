# cuCFD: CUDA-accelerated Computational Fluid Dynamics Solver

A high-performance CFD solver using Fortran with GPU acceleration via OpenMP and CUDA.

## Features

- Finite Volume discretization schemes
- CUDA and OpenMP acceleration for massively parallel execution
- Multi-device support for scalability
- Built-in performance benchmarking (CPU vs GPU)
- Modular design for extensibility

## Project Structure

```
cuCFD/
├── src/
│   ├── benchmarks/             # Performance benchmark codes
│   │   ├── gpu_benchmark.f90   # GPU vs CPU benchmarking
│   │   └── gpu_utils.f90       # GPU utility module
│   ├── core/                   # Core solver components
│   │   ├── mesh.f90            # Mesh data structures
│   │   ├── field.f90           # Field data structures  
│   │   ├── boundary.f90        # Boundary conditions
│   │   └── parallel.f90        # Parallelization utilities
│   ├── numerics/               # Numerical methods
│   │   ├── gradient.f90        # Gradient calculation
│   │   ├── divergence.f90      # Divergence operators
│   │   ├── laplacian.f90       # Laplacian operators
│   │   ├── interpolation.f90   # Field interpolation
│   │   └── linear_solvers.f90  # Linear equation solvers
│   ├── models/                 # Physical models
│   │   ├── advection.f90       # Advection model
│   │   ├── diffusion.f90       # Diffusion model
│   │   ├── turbulence.f90      # Turbulence models
│   │   └── heat_transfer.f90   # Heat transfer models
│   ├── io/                     # Input/Output utilities
│   │   ├── reader.f90          # Mesh/data readers
│   │   └── writer.f90          # Result output
│   ├── tests/                  # Test cases
│   │   ├── taylor_green.f90    # Taylor-Green vortex
│   │   └── lid_driven_cavity.f90 # Lid-driven cavity flow
│   ├── utils/                  # Utility functions
│   │   ├── constants.f90       # Physical and numerical constants
│   │   └── error.f90           # Error handling
│   └── main.f90                # Main program
├── include/                    # Header files
├── build/                      # Build directory
├── bin/                        # Executable directory
├── examples/                   # Example cases
│   ├── 2d_cavity/              # 2D lid-driven cavity
│   ├── 3d_channel/             # 3D channel flow
│   └── heat_transfer/          # Heat transfer examples
├── tests/                      # Test suite
├── tools/                      # Helper tools
├── Makefile                    # Build system
└── README.md                   # This file
```

## Building

Different build configurations are available:

```bash
# Default build with CUDA acceleration
make

# Debug build with bounds checking
make debug

# Release build with optimizations
make release

# Clean build files
make clean

# Build for specific GPU architectures
make tesla  # For Tesla GPUs (cc75)
make volta  # For Volta GPUs (cc70)
make ampere # For Ampere GPUs (cc80)
```

## Performance Optimization

This CFD solver is optimized for high performance:
- GPU-accelerated core algorithms with OpenMP offloading
- Memory access patterns optimized for GPU coalescing
- Efficient parallel data structures
- Block-based computation for cache efficiency
- Mixed-precision computation where appropriate
- Multi-device support for larger problems

## Getting Started

```bash
# Clone the repository
git clone https://github.com/yourusername/cuCFD.git
cd cuCFD

# Build the project
make release

# Run a test case
./bin/cuCFD examples/2d_cavity/cavity.config
```

## Requirements

- NVIDIA HPC SDK (including nvfortran compiler)
- CUDA-capable GPU
- OpenMP support
- Make utility

## License

This project is open-source and available under the MIT License.
