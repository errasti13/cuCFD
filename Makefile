# Compiler and flags
FC = nvfortran

# Configurable options
USE_GPU ?= yes
GPU_ARCH ?= cc75 # Default GPU architecture

# Base Flags (Common)
BASE_FFLAGS = -Mpreprocess
OPTFLAGS = -O3 -fast -Mvect
DEBUG_FLAGS = -g -O0 -Minfo=all -Mbounds -traceback -Mchkptr -Mchkstk -Mdclchk
RELEASE_FLAGS = $(OPTFLAGS) -Minfo=accel # Keep accel info for release builds for now

# Module availability flags for preprocessor
MODULE_FLAGS = -DHAVE_DIFFUSION -DHAVE_INCOMPRESSIBLE_FLOW

# Conditional Flags
ifeq ($(USE_GPU), yes)
    GPU_FFLAGS := $(BASE_FFLAGS) $(MODULE_FLAGS) -mp=gpu -cuda -gpu=$(GPU_ARCH)
    GPU_LDFLAGS := -gpu=$(GPU_ARCH)
    FINAL_FFLAGS := $(GPU_FFLAGS)
    FINAL_LDFLAGS := $(GPU_LDFLAGS)
    BUILD_MODE_MSG := "Building for GPU (arch: $(GPU_ARCH))"
else
    CPU_FFLAGS := $(BASE_FFLAGS) $(MODULE_FLAGS) # Add CPU-specific flags like -mp=multicore if desired later
    CPU_LDFLAGS := # Add CPU-specific linker flags if needed later
    FINAL_FFLAGS := $(CPU_FFLAGS)
    FINAL_LDFLAGS := $(CPU_LDFLAGS)
    BUILD_MODE_MSG := "Building for CPU"
endif

# Directories
SRCDIR = src
COREDIR = $(SRCDIR)/core
NUMERICSDIR = $(SRCDIR)/numerics
MODELSDIR = $(SRCDIR)/models
IODIR = $(SRCDIR)/io
TESTSDIR = $(SRCDIR)/tests
UTILSDIR = $(SRCDIR)/utils
BENCHMARKSDIR = $(SRCDIR)/benchmarks
BUILDDIR = build
BINDIR = bin
MODDIR = $(BUILDDIR)/modules

# Sources
CORE_SRCS = $(wildcard $(COREDIR)/*.f90)
NUMERICS_SRCS = $(wildcard $(NUMERICSDIR)/*.f90)
MODELS_SRCS = $(wildcard $(MODELSDIR)/*.f90)
CONVERTER_SRC = $(SRCDIR)/io/solution_converter.f90
# Filter out solution_converter.f90 from IO_SRCS
IO_SRCS = $(filter-out $(CONVERTER_SRC),$(wildcard $(IODIR)/*.f90))
TESTS_SRCS = $(wildcard $(TESTSDIR)/*.f90)
UTILS_SRCS = $(wildcard $(UTILSDIR)/*.f90)
BENCHMARKS_SRCS = $(wildcard $(BENCHMARKSDIR)/*.f90)
MAIN_SRC = $(SRCDIR)/main.f90

# Objects
CORE_OBJS = $(CORE_SRCS:$(COREDIR)/%.f90=$(BUILDDIR)/%.o)
NUMERICS_OBJS = $(NUMERICS_SRCS:$(NUMERICSDIR)/%.f90=$(BUILDDIR)/%.o)
MODELS_OBJS = $(MODELS_SRCS:$(MODELSDIR)/%.f90=$(BUILDDIR)/%.o)
IO_OBJS = $(IO_SRCS:$(IODIR)/%.f90=$(BUILDDIR)/%.o)
TESTS_OBJS = $(TESTS_SRCS:$(TESTSDIR)/%.f90=$(BUILDDIR)/%.o)
UTILS_OBJS = $(UTILS_SRCS:$(UTILSDIR)/%.f90=$(BUILDDIR)/%.o)
BENCHMARKS_OBJS = $(BENCHMARKS_SRCS:$(BENCHMARKSDIR)/%.f90=$(BUILDDIR)/%.o)
MAIN_OBJ = $(MAIN_SRC:$(SRCDIR)/%.f90=$(BUILDDIR)/%.o)
CONVERTER_OBJ = $(CONVERTER_SRC:$(IODIR)/%.f90=$(BUILDDIR)/%.o)

# All objects
OBJS = $(CORE_OBJS) $(NUMERICS_OBJS) $(MODELS_OBJS) $(IO_OBJS) $(TESTS_OBJS) $(UTILS_OBJS) $(MAIN_OBJ)

# Objects for converter
CONVERTER_ALL_OBJS = $(CORE_OBJS) $(NUMERICS_OBJS) $(MODELS_OBJS) $(IO_OBJS) $(UTILS_OBJS) $(CONVERTER_OBJ)

# Benchmark objects (Assuming these might also be conditional later?)
BENCHMARK_ALL_OBJS = $(CORE_OBJS) $(NUMERICS_OBJS) $(MODELS_OBJS) $(IO_OBJS) $(UTILS_OBJS) $(BENCHMARKS_OBJS)

# Targets
.PHONY: all debug release clean tesla volta ampere turing benchmarks setup build_info converter

# Module directory handling
$(shell mkdir -p $(BINDIR) $(BUILDDIR) $(MODDIR))

all: setup build_info $(BINDIR)/cuCFD $(BINDIR)/solution_converter

benchmarks: setup build_info $(BINDIR)/benchmark

# Add converter target
converter: setup build_info $(BINDIR)/solution_converter

# Modify debug/release to append to FINAL_FFLAGS
debug: FINAL_FFLAGS += $(DEBUG_FLAGS)
debug: all

release: FINAL_FFLAGS += $(RELEASE_FLAGS)
release: all

setup:
	@mkdir -p $(BINDIR) $(BUILDDIR) $(MODDIR)

build_info:
	@echo $(BUILD_MODE_MSG)

# Use FINAL_FFLAGS and FINAL_LDFLAGS for linking
$(BINDIR)/cuCFD: $(OBJS)
	$(FC) $(FINAL_FFLAGS) $(FINAL_LDFLAGS) -module $(MODDIR) -o $@ $^

# Use FINAL_FFLAGS and FINAL_LDFLAGS for solution converter
$(BINDIR)/solution_converter: $(CONVERTER_ALL_OBJS)
	$(FC) $(FINAL_FFLAGS) $(FINAL_LDFLAGS) -module $(MODDIR) -o $@ $^

# Use FINAL_FFLAGS and FINAL_LDFLAGS for benchmark linking
# Link benchmark against all necessary objects
$(BINDIR)/benchmark: $(BENCHMARK_ALL_OBJS)
	$(FC) $(FINAL_FFLAGS) $(FINAL_LDFLAGS) -module $(MODDIR) -o $@ $^

# Module dependencies - needed for proper build order
# Define specific dependencies for modules
$(BUILDDIR)/field.o: $(BUILDDIR)/mesh.o
$(BUILDDIR)/boundary.o: $(BUILDDIR)/mesh.o $(BUILDDIR)/field.o
$(BUILDDIR)/diffusion.o: $(BUILDDIR)/mesh.o $(BUILDDIR)/field.o
$(BUILDDIR)/incompressible_flow.o: $(BUILDDIR)/mesh.o $(BUILDDIR)/field.o $(BUILDDIR)/boundary.o
$(BUILDDIR)/config_parser.o: $(BUILDDIR)/mesh.o
$(BUILDDIR)/solution_io.o: $(BUILDDIR)/mesh.o $(BUILDDIR)/field.o
$(BUILDDIR)/solver_manager.o: $(BUILDDIR)/mesh.o $(BUILDDIR)/field.o $(BUILDDIR)/boundary.o $(BUILDDIR)/config_parser.o $(BUILDDIR)/diffusion.o $(BUILDDIR)/incompressible_flow.o $(BUILDDIR)/solution_io.o
$(MAIN_OBJ): $(BUILDDIR)/mesh.o $(BUILDDIR)/field.o $(BUILDDIR)/config_parser.o $(BUILDDIR)/solver_manager.o $(BUILDDIR)/solution_io.o

# Use FINAL_FFLAGS for compilation rules
# Core module rules
$(BUILDDIR)/%.o: $(COREDIR)/%.f90
	$(FC) $(FINAL_FFLAGS) -module $(MODDIR) -c $< -o $@

# Numerics module rules
$(BUILDDIR)/%.o: $(NUMERICSDIR)/%.f90
	$(FC) $(FINAL_FFLAGS) -module $(MODDIR) -c $< -o $@

# Models module rules
$(BUILDDIR)/%.o: $(MODELSDIR)/%.f90
	$(FC) $(FINAL_FFLAGS) -module $(MODDIR) -c $< -o $@

# IO module rules
$(BUILDDIR)/%.o: $(IODIR)/%.f90
	$(FC) $(FINAL_FFLAGS) -module $(MODDIR) -c $< -o $@

# Tests module rules
$(BUILDDIR)/%.o: $(TESTSDIR)/%.f90
	$(FC) $(FINAL_FFLAGS) -module $(MODDIR) -c $< -o $@

# Utils module rules
$(BUILDDIR)/%.o: $(UTILSDIR)/%.f90
	$(FC) $(FINAL_FFLAGS) -module $(MODDIR) -c $< -o $@

# Benchmarks module rules
$(BUILDDIR)/%.o: $(BENCHMARKSDIR)/%.f90
	$(FC) $(FINAL_FFLAGS) -module $(MODDIR) -c $< -o $@

# Main program rule
$(BUILDDIR)/%.o: $(SRCDIR)/%.f90
	$(FC) $(FINAL_FFLAGS) -module $(MODDIR) -c $< -o $@

# GPU architecture targets
# These still function by setting GPU_ARCH and implicitly calling 'all'
# The conditional logic will pick up the correct GPU_ARCH if USE_GPU=yes
turing: GPU_ARCH = cc75
turing: all

tesla: GPU_ARCH = cc75
tesla: all

volta: GPU_ARCH = cc70
volta: all

ampere: GPU_ARCH = cc80
ampere: all

clean:
	rm -rf $(BUILDDIR) $(BINDIR)
