# Compiler and flags
FC = nvfortran
FCFLAGS = -mp=gpu -cuda -Mpreprocess
OPTFLAGS = -O3 -fast -Mvect
DEBUG_FLAGS = -g -O0 -Minfo=all -Mbounds -Mtrace
RELEASE_FLAGS = $(OPTFLAGS) -Minfo=accel

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

# Default GPU architecture (Tesla - cc75)
GPU_ARCH = cc75

# Sources
CORE_SRCS = $(wildcard $(COREDIR)/*.f90)
NUMERICS_SRCS = $(wildcard $(NUMERICSDIR)/*.f90)
MODELS_SRCS = $(wildcard $(MODELSDIR)/*.f90)
IO_SRCS = $(wildcard $(IODIR)/*.f90)
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

# All objects
OBJS = $(CORE_OBJS) $(NUMERICS_OBJS) $(MODELS_OBJS) $(IO_OBJS) $(TESTS_OBJS) $(UTILS_OBJS) $(MAIN_OBJ)

# Benchmark objects
BENCHMARK_OBJS = $(BENCHMARKS_OBJS)

# Targets
.PHONY: all debug release clean tesla volta ampere turing benchmarks setup

# Module directory handling
$(shell mkdir -p $(BINDIR) $(BUILDDIR) $(MODDIR))

all: setup $(BINDIR)/cuCFD

benchmarks: setup $(BINDIR)/benchmark

debug: FCFLAGS += $(DEBUG_FLAGS)
debug: all

release: FCFLAGS += $(RELEASE_FLAGS)
release: all

setup:
	@mkdir -p $(BINDIR) $(BUILDDIR) $(MODDIR)

$(BINDIR)/cuCFD: $(OBJS)
	$(FC) $(FCFLAGS) -gpu=$(GPU_ARCH) -module $(MODDIR) -o $@ $^

$(BINDIR)/benchmark: $(BENCHMARK_OBJS)
	$(FC) $(FCFLAGS) -gpu=$(GPU_ARCH) -module $(MODDIR) -o $@ $^

# Module dependencies - needed for proper build order
# Define specific dependencies for modules
$(BUILDDIR)/field.o: $(BUILDDIR)/mesh.o
$(BUILDDIR)/boundary.o: $(BUILDDIR)/mesh.o $(BUILDDIR)/field.o
$(BUILDDIR)/diffusion.o: $(BUILDDIR)/mesh.o $(BUILDDIR)/field.o
$(MAIN_OBJ): $(BUILDDIR)/mesh.o $(BUILDDIR)/field.o $(BUILDDIR)/boundary.o $(BUILDDIR)/diffusion.o
$(BUILDDIR)/lid_driven_cavity.o: $(BUILDDIR)/mesh.o $(BUILDDIR)/field.o $(BUILDDIR)/boundary.o

# Core module rules
$(BUILDDIR)/%.o: $(COREDIR)/%.f90
	$(FC) $(FCFLAGS) -gpu=$(GPU_ARCH) -module $(MODDIR) -c $< -o $@

# Numerics module rules
$(BUILDDIR)/%.o: $(NUMERICSDIR)/%.f90
	$(FC) $(FCFLAGS) -gpu=$(GPU_ARCH) -module $(MODDIR) -c $< -o $@

# Models module rules
$(BUILDDIR)/%.o: $(MODELSDIR)/%.f90
	$(FC) $(FCFLAGS) -gpu=$(GPU_ARCH) -module $(MODDIR) -c $< -o $@

# IO module rules
$(BUILDDIR)/%.o: $(IODIR)/%.f90
	$(FC) $(FCFLAGS) -gpu=$(GPU_ARCH) -module $(MODDIR) -c $< -o $@

# Tests module rules
$(BUILDDIR)/%.o: $(TESTSDIR)/%.f90
	$(FC) $(FCFLAGS) -gpu=$(GPU_ARCH) -module $(MODDIR) -c $< -o $@

# Utils module rules
$(BUILDDIR)/%.o: $(UTILSDIR)/%.f90
	$(FC) $(FCFLAGS) -gpu=$(GPU_ARCH) -module $(MODDIR) -c $< -o $@

# Benchmarks module rules
$(BUILDDIR)/%.o: $(BENCHMARKSDIR)/%.f90
	$(FC) $(FCFLAGS) -gpu=$(GPU_ARCH) -module $(MODDIR) -c $< -o $@

# Main program rule
$(BUILDDIR)/%.o: $(SRCDIR)/%.f90
	$(FC) $(FCFLAGS) -gpu=$(GPU_ARCH) -module $(MODDIR) -c $< -o $@

# GPU architecture targets
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
