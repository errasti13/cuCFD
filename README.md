# OpenMP GPU vs CPU Performance Benchmark

A Fortran program to benchmark and compare CPU and GPU performance using OpenMP offloading.

## Prerequisites

- NVIDIA HPC SDK (including nvfortran compiler)
- CUDA-capable GPU
- OpenMP support
- Make utility

## Project Structure

```
testCudaFortran/
├── src/
│   ├── gpu_utils.f90      # GPU utility module
│   └── testCudaFortran.f90 # Main program
├── Makefile
└── README.md
```

## Building

Different build configurations are available:

```bash
# Default build (Tesla GPU)
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

## Program Features

- Compares CPU and GPU performance for array operations
- Supports multiple array sizes
- Performs warmup iterations for accurate GPU timing
- Validates results between CPU and GPU computations
- Uses OpenMP for both CPU parallelization and GPU offloading
- Configurable block sizes for GPU optimization

## Configuration

Key parameters can be modified in `src/testCudaFortran.f90`:

- `NUM_SIZES`: Number of different array sizes to test
- `ITERATIONS`: Number of timing iterations
- `WARMUP_ITER`: Number of GPU warmup iterations
- `BLOCK_SIZE`: GPU thread block size
- `sizes`: Array of problem sizes to test

## Output Format

The program outputs:
- System information (CPU threads, GPU availability)
- Performance comparison table with:
  - Array size
  - CPU execution time
  - GPU execution time
  - Speedup ratio
  - Result validation

## Example Output

```
Number of CPU threads: 16
GPU execution enabled: T

Size        CPU Time(s)    GPU Time(s)    Speedup    Match?
--------------------------------------------------------
 10000000      0.006043      0.011628       0.52       T
 50000000      0.029200      0.050346       0.58       T
100000000      0.066245      0.121602       0.54       T
200000000      0.119392      0.251542       0.47       T
```

## License

This project is open-source and available under the MIT License.
