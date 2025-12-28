# =============================================================================
# TongueTwister 2.0 Makefile
# Cross-platform build with conditional macOS Accelerate/Metal support
# =============================================================================

# Compiler settings
CXX = g++
OBJCXX = clang++

# Base flags (all platforms)
# Note: -ffast-math enables -ffinite-math-only which breaks std::numeric_limits<>::infinity()
# used by json.hpp. We add -fno-finite-math-only to restore infinity/NaN support while
# keeping other fast-math optimizations (reciprocal-math, no-signed-zeros, etc.)
CXXFLAGS_BASE = -std=c++20 -O3 -march=native -flto -DNDEBUG -ffast-math -fno-finite-math-only -pthread -m64

# Warnings (optional, uncomment for development)
CXXFLAGS_BASE += -Wall -Wextra -Wpedantic

# Target
TARGET = tonguetwister

# =============================================================================
# Platform Detection
# =============================================================================
UNAME_S := $(shell uname -s)

ifeq ($(UNAME_S),Darwin)
    # -----------------------------------------------------------------
    # macOS: Full support with Accelerate framework and Metal GPU
    # -----------------------------------------------------------------
    PLATFORM = macos
    
    # Use clang++ on macOS for better Objective-C++ support
    CXX = clang++
    
    # macOS-specific flags
    # -DACCELERATE_NEW_LAPACK: Use updated ILP64 LAPACK interface (macOS 13.3+)
    # -DHAVE_ACCELERATE: Enable Accelerate framework features
    # -DHAVE_METAL: Enable Metal GPU compute
    CXXFLAGS_PLATFORM = -DHAVE_ACCELERATE -DHAVE_METAL -DACCELERATE_NEW_LAPACK
    
    # Frameworks
    LDFLAGS = -framework Accelerate -framework Metal -framework Foundation
    
    # Objective-C++ flags for .mm files
    OBJCXXFLAGS = $(CXXFLAGS_BASE) $(CXXFLAGS_PLATFORM) -fobjc-arc
    
    # macOS-specific source files
    MACOS_SRCS = MathCacheAccelerated.cpp \
                 CtmcAccelerated.cpp \
                 TransitionProbabilitiesGpu.cpp
    
    # Objective-C++ source files (Metal compute shaders)
    OBJCXX_SRCS = GPUMatrixExponentialBatch.mm
    
else ifeq ($(UNAME_S),Linux)
    # -----------------------------------------------------------------
    # Linux: CPU-only with pthreads
    # -----------------------------------------------------------------
    PLATFORM = linux
    
    # Linux-specific flags
    CXXFLAGS_PLATFORM = -DNO_ACCELERATE -DNO_METAL
    
    # Link against pthread (already in base flags, but explicit for clarity)
    LDFLAGS = -lpthread
    
    # No macOS-specific sources on Linux
    MACOS_SRCS =
    OBJCXX_SRCS =
    
else
    # -----------------------------------------------------------------
    # Other Unix-like systems: Generic fallback
    # -----------------------------------------------------------------
    PLATFORM = generic
    
    CXXFLAGS_PLATFORM = -DNO_ACCELERATE -DNO_METAL
    LDFLAGS = -lpthread
    
    MACOS_SRCS =
    OBJCXX_SRCS =
endif

# Combined C++ flags
CXXFLAGS = $(CXXFLAGS_BASE) $(CXXFLAGS_PLATFORM)

# =============================================================================
# Source Files
# =============================================================================

# Core source files (all platforms)
CORE_SRCS = \
    Alignment.cpp \
    BitMask.cpp \
    BitSet.cpp \
    BranchMapping.cpp \
    Container.cpp \
    Ctmc.cpp \
    JsonData.cpp \
    LikelihoodCalculator.cpp \
    MathCache.cpp \
    MatrixPool.cpp \
    Mcmc.cpp \
    McmcMarginal.cpp \
    McmcOutput.cpp \
    McmcPhase.cpp \
    McmcPosterior.cpp \
    McmcTimer.cpp \
    Model.cpp \
    Msg.cpp \
    Node.cpp \
    NodeComparator.cpp \
    NodeFactory.cpp \
    NodeSet.cpp \
    Parameter.cpp \
    ParameterAlignment.cpp \
    ParameterExchangeabilities.cpp \
    ParameterFrequencies.cpp \
    ParameterIndelRates.cpp \
    ParameterTree.cpp \
    Partition.cpp \
    Probability.cpp \
    RandomVariable.cpp \
    RateMatrix.cpp \
    RateMatrixHelper.cpp \
    States.cpp \
    SteppingStones.cpp \
    Stone.cpp \
    Subset.cpp \
    Threads.cpp \
    TransitionMatrixMap.cpp \
    TransitionProbabilitiesCpu.cpp \
    TransitionProbabilityCalculator.cpp \
    Tree.cpp \
    Update.cpp \
    UpdateAlignment.cpp \
    UpdateBranchLength.cpp \
    UpdateExchangeabilities.cpp \
    UpdateFrequencies.cpp \
    UpdateIndelRates.cpp \
    UpdateManager.cpp \
    UpdateTopology.cpp \
    UserSettings.cpp \
    main.cpp

# All C++ sources for this platform
ALL_CPP_SRCS = $(CORE_SRCS) $(MACOS_SRCS)

# =============================================================================
# Object Files
# =============================================================================

# C++ object files
CPP_OBJS = $(ALL_CPP_SRCS:.cpp=.o)

# Objective-C++ object files (macOS only)
OBJCXX_OBJS = $(OBJCXX_SRCS:.mm=.o)

# All object files
OBJS = $(CPP_OBJS) $(OBJCXX_OBJS)

# =============================================================================
# Header Dependencies (key headers that trigger rebuilds)
# =============================================================================

# Header-only libraries
HEADER_ONLY = json.hpp argparse.hpp

# Core headers
CORE_HEADERS = \
    Container.hpp \
    MathCache.hpp \
    Threads.hpp \
    TransitionProbabilities.hpp

# =============================================================================
# Build Rules
# =============================================================================

.PHONY: all clean info

all: info $(TARGET)

# Link
$(TARGET): $(OBJS)
	@echo "Linking $(TARGET) for $(PLATFORM)..."
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)
	@echo "Build complete: $(TARGET)"

# Compile C++ sources
%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Compile Objective-C++ sources (macOS only)
%.o: %.mm
	$(OBJCXX) $(OBJCXXFLAGS) -c $< -o $@

# Clean
clean:
	rm -f $(OBJS)
	@echo "Cleaned build artifacts"

# Build info
info:
	@echo "=============================================="
	@echo "TongueTwister 2.0 Build Configuration"
	@echo "=============================================="
	@echo "Platform:    $(PLATFORM) ($(UNAME_S))"
	@echo "Compiler:    $(CXX)"
	@echo "C++ Flags:   $(CXXFLAGS)"
ifeq ($(PLATFORM),macos)
	@echo "Frameworks:  Accelerate, Metal, Foundation"
	@echo "Features:    BLAS/LAPACK acceleration, GPU compute"
else
	@echo "Features:    CPU threading only"
endif
	@echo "=============================================="
	@echo ""

# =============================================================================
# Explicit Dependencies (optional - for incremental builds)
# =============================================================================

# Files that depend on Container.hpp (DoubleMatrix, DoubleArray)
Container.o: Container.cpp Container.hpp
MathCache.o: MathCache.cpp MathCache.hpp Container.hpp

# Files that depend on MathCache.hpp
LikelihoodCalculator.o: LikelihoodCalculator.cpp LikelihoodCalculator.hpp MathCache.hpp
TransitionProbabilityCalculator.o: TransitionProbabilityCalculator.cpp TransitionProbabilityCalculator.hpp

# Files that depend on Threads.hpp
Model.o: Model.cpp Model.hpp Threads.hpp
Threads.o: Threads.cpp Threads.hpp MathCache.hpp

# macOS-specific dependencies
ifeq ($(PLATFORM),macos)
MathCacheAccelerated.o: MathCacheAccelerated.cpp MathCacheAccelerated.hpp MathCache.hpp Container.hpp
CtmcAccelerated.o: CtmcAccelerated.cpp CtmcAccelerated.hpp Ctmc.hpp
GPUMatrixExponentialBatch.o: GPUMatrixExponentialBatch.mm GPUMatrixExponentialBatch.hpp Container.hpp
TransitionProbabilitiesGpu.o: TransitionProbabilitiesGpu.cpp TransitionProbabilitiesGpu.hpp TransitionProbabilities.hpp MatrixPool.hpp
endif

# TransitionProbabilities dependencies
TransitionProbabilitiesCpu.o: TransitionProbabilitiesCpu.cpp TransitionProbabilitiesCpu.hpp TransitionProbabilities.hpp MatrixPool.hpp
MatrixPool.o: MatrixPool.cpp MatrixPool.hpp Container.hpp
TransitionMatrixMap.o: TransitionMatrixMap.cpp TransitionMatrixMap.hpp Container.hpp MatrixPool.hpp

# =============================================================================
# Development/Debug Targets
# =============================================================================

.PHONY: debug release

# Debug build (with symbols, no optimization)
debug: CXXFLAGS_BASE = -std=c++20 -g -O0 -pthread -m64 -Wall -Wextra
debug: clean all

# Release build (default, but explicit)
release: all

# =============================================================================
# Installation (optional)
# =============================================================================

PREFIX ?= /usr/local

.PHONY: install uninstall

install: $(TARGET)
	install -d $(PREFIX)/bin
	install -m 755 $(TARGET) $(PREFIX)/bin/

uninstall:
	rm -f $(PREFIX)/bin/$(TARGET)
