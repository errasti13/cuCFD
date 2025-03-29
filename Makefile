# Compiler and flags
FC = nvfortran
FCFLAGS = -mp=gpu -cuda -Mpreprocess
OPTFLAGS = -O3 -fast -Mvect
DEBUG_FLAGS = -g -O0 -Minfo=all -Mbounds -Mtrace
RELEASE_FLAGS = $(OPTFLAGS) -Minfo=accel

# Directories
SRCDIR = src
BUILDDIR = build
BINDIR = bin
MODDIR = $(BUILDDIR)/modules

# Default GPU architecture (Turing - GTX 1650)
GPU_ARCH = cc75

# Sources and objects
SRCS = $(wildcard $(SRCDIR)/*.f90)
OBJS = $(SRCS:$(SRCDIR)/%.f90=$(BUILDDIR)/%.o)

# Module directory handling
$(shell mkdir -p $(BINDIR) $(BUILDDIR) $(MODDIR))

# Targets
all: setup $(BINDIR)/program

debug: FCFLAGS += $(DEBUG_FLAGS)
debug: all

release: FCFLAGS += $(RELEASE_FLAGS)
release: all

setup:
	@mkdir -p $(BINDIR) $(BUILDDIR) $(MODDIR)

$(BINDIR)/program: $(OBJS)
	$(FC) $(FCFLAGS) -gpu=$(GPU_ARCH) -module $(MODDIR) -o $@ $^

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

.PHONY: all debug release clean tesla volta ampere turing setup
