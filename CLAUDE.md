# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Build and Test Commands

- **Build**: `./scripts/build.sh` - Builds the project using SimFactory with the Cactus computational toolkit
- **Test**: `./scripts/test.sh` - Runs the test suite with both single-process and multi-process configurations
- **Download dependencies**: `./scripts/download.sh` - Downloads and sets up the Cactus framework and dependencies

The build system uses environment variables:
- `ACCELERATOR`: cpu, cuda, or rocm (determines hardware acceleration)
- `REAL_PRECISION`: real64 (sets floating point precision)  
- `MODE`: debug or optimize (sets compilation mode)

## Architecture Overview

MCTX is a Monte Carlo Transport framework built on top of CarpetX, which is part of the Einstein Toolkit for computational relativity. The codebase follows the Cactus computational toolkit architecture.

### Core Components

**Thorns (Cactus modules):**
- `Particles/`: Main particle simulation thorn containing AMReX-based particle containers for Monte Carlo transport
- `TestParticles/`: Test suite for basic particle functionality
- `TestNuParticles/`: Test suite for neutrino particle simulations

**Key Architecture Patterns:**
- Each thorn has CCL files defining configuration, interfaces, parameters, and scheduling
- C++ source code uses AMReX for particle management and GPU acceleration
- SimFactory handles build configuration and job submission
- Uses CarpetX for adaptive mesh refinement and parallelization

### Thorn Structure
Each thorn follows the standard Cactus format:
- `configuration.ccl`: Dependencies and thorn requirements
- `interface.ccl`: Variable declarations and includes  
- `param.ccl`: Runtime parameters
- `schedule.ccl`: Function scheduling and dependencies
- `src/`: C++ implementation files

### Dependencies
The project requires the full Einstein Toolkit stack including:
- CarpetX for adaptive mesh refinement
- AMReX for particle containers and GPU support  
- Various ExternalLibraries (HDF5, MPI, etc.)

## Development Workflow

When developing new features:
1. Modify thorn CCL files to declare new variables/parameters
2. Implement functionality in C++ source files in `src/`
3. Use the provided test thorns as examples for new functionality
4. Build and test using the provided scripts

The codebase uses GPU-aware code with AMReX, so consider GPU compatibility when making changes.