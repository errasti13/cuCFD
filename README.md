# cuCFD - CUDA-accelerated Computational Fluid Dynamics Solver

A GPU-accelerated CFD solver written in modern Fortran with CUDA/OpenMP target offload.

## Features

- GPU acceleration using CUDA via OpenMP target offload for parallel computations
- Conditional CPU/GPU build capability
- Heat conduction solver with boundary conditions
- Solution I/O for saving and loading simulation results
- Paraview-compatible VTK export for visualization

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

### Requirements

- NVIDIA HPC SDK (nvfortran compiler) with CUDA
- NVIDIA GPU with compute capability 7.0 or higher (Volta/Turing/Ampere/Hopper)

### Compilation

Build with GPU support (default):
```bash
make
```

Build for CPU only:
```bash
make USE_GPU=no
```

Build for specific GPU architecture:
```bash
make tesla     # For Tesla V100 (cc70) or T4 (cc75)
make volta     # For Volta (cc70)
make turing    # For Turing (cc75)
make ampere    # For Ampere (cc80)
```

Build in debug mode:
```bash
make debug
```

## Running Simulations

Run the application with a configuration file:
```bash
./bin/cuCFD examples/2d_cavity/lid_driven_cavity.config
```

Example configuration files are provided in the `examples/` directory.

## Solution I/O

cuCFD can save simulation results to binary files and convert them to Paraview-compatible formats.

### Writing Solution Files

Solution writing is controlled via configuration files:
```
output_dir: results/my_simulation
output_frequency: 10  # Write every 10 time steps
write_solution: true  # Enable solution writing
```

### Converting Solution Files for Visualization

Use the solution converter to transform binary solution files to VTK format:
```bash
./bin/solution_converter results/my_simulation/solution_000100.bin -o temperature.vtk
```

Options:
- `-o, --output FILE` - Specify output file name
- `-d, --dir DIRECTORY` - Specify output directory
- `-f, --format FORMAT` - Specify output format (vtk or vtu, default: vtk)

## Examples

Heat conduction example:
```bash
./bin/cuCFD examples/heat_conduction/heat_solution_config.txt
```

This will simulate heat transfer and save solution files that can be converted to VTK format for visualization.

## License

This project is open source and licensed under the MIT License.
