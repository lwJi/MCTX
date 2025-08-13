# MCTX Documentation

Monte Carlo Transport based on CarpetX - A high-performance particle simulation framework for computational physics.

## Table of Contents

1. [Overview](#overview)
2. [Installation & Setup](#installation--setup)
3. [Architecture](#architecture)
4. [API Reference](#api-reference)
5. [Usage Examples](#usage-examples)
6. [Testing](#testing)
7. [Development](#development)

## Overview

MCTX (Monte Carlo Transport based on CarpetX) is a computational physics framework designed for particle simulations in relativistic environments. Built on top of the Einstein Toolkit's CarpetX infrastructure, it provides:

- **High-performance particle tracking** using AMReX particle containers
- **GPU acceleration** support for CUDA and ROCm
- **Adaptive mesh refinement** through CarpetX integration
- **Relativistic geodesic integration** for particle trajectories
- **Scalable parallel computing** with MPI

The framework is specifically designed for neutrino transport simulations and relativistic particle physics applications.

## Installation & Setup

### Prerequisites

- Einstein Toolkit environment
- CarpetX framework
- AMReX library
- MPI implementation
- Optional: CUDA or ROCm for GPU acceleration

## Architecture

### Core Components

#### Thorns (Cactus Modules)

**Particles** - Main simulation engine
- Location: `Particles/`
- Purpose: Core particle container and physics implementation
- Key files:
  - `NuParticleContainers.hxx/cxx` - Particle container implementation
  - `wolfram/` - Auto-generated physics kernels

**TestNuParticles** - Neutrino-specific tests
- Location: `TestNuParticles/`
- Purpose: Advanced neutrino transport testing

#### Thorn Structure

Each thorn follows the Cactus framework convention:

```
ThornName/
├── configuration.ccl  # Dependencies and requirements
├── interface.ccl      # Variable declarations and includes
├── param.ccl          # Runtime parameters
├── schedule.ccl       # Function scheduling
└── src/              # C++ implementation
    ├── make.code.defn # Build configuration
    └── *.cxx/*.hxx    # Source files
```

### Physics Implementation

#### Particle Container Design

```cpp
struct PIdx {
  enum {
    px = 0, py, pz,        // Physical momentum (SoA)
    x0, y0, z0,           // Saved position at t^n  
    px0, py0, pz0,        // Saved momentum at t^n
    nattribs
  };
};

using Container = amrex::AmrParticleContainer<0, 0, PIdx::nattribs, 0>;
```

#### Integration Scheme

MCTX uses a second-order Runge-Kutta (midpoint) method for geodesic integration:

1. **Substep 1**: Compute k₁ = f(y^n), advance to midpoint y^(n+1/2)
2. **Redistribute**: Move particles across MPI domains/AMR levels
3. **Substep 2**: Compute k₂ = f(y^(n+1/2)), final step y^(n+1) = y^n + dt*k₂

## API Reference

### NuParticleContainer Class

```cpp
class NuParticleContainer : public amrex::AmrParticleContainer<...> {
public:
    // Constructor
    NuParticleContainer(amrex::AmrCore *amr_core);
    
    // Main evolution step
    void PushAndDeposeParticles(const amrex::MultiFab &lapse,
                               const amrex::MultiFab &shift,
                               const amrex::MultiFab &met3d, 
                               CCTK_REAL dt, const int lev);
    
    // Output methods
    void OutputParticlesAscii(CCTK_ARGUMENTS);
    void OutputParticlesPlot(CCTK_ARGUMENTS);
};
```

### Key Methods

#### PushAndDeposeParticles
- **Purpose**: Evolve particles using RK2 geodesic integration
- **Parameters**:
  - `lapse`: ADM lapse function α
  - `shift`: ADM shift vector βⁱ  
  - `met3d`: 3-metric γᵢⱼ
  - `dt`: Time step
  - `lev`: AMR level

#### Field Interpolation
```cpp
void gather_fields_calcrhs_at_pos(VectR &dtmom, VectR &dtpos,
                                  const VectR &mom, const VectR &pos,
                                  /* metric fields */);
```
Trilinear interpolation of metric fields to particle positions with first-order derivative computation.

### Physics Kernels

Generated Wolfram Language code in `wolfram/`:
- `particles_geodesic.wl` - Geodesic equation derivation
- `particles_derivsinline.hxx` - Optimized finite difference stencils

## Usage Examples

### Runtime Configuration

Key parameters from `param.ccl`:
- `out_tsv_every`: ASCII output frequency
- `out_plot_every`: Plot file output frequency  
- Grid setup through CarpetX parameters

### Output Formats

1. **ASCII**: Human-readable particle data
2. **Plot files**: Binary format for visualization
3. **Performance data**: Timing and scaling metrics

## Development

### Adding New Physics

1. **Modify thorn CCL files** to declare new variables/parameters
2. **Implement C++ kernels** in `src/` directory
3. **Update Wolfram derivations** if needed for new physics
4. **Add tests** in appropriate test thorn
5. **Build and validate** using provided scripts

### GPU Development

- Use `CCTK_DEVICE` decorators for GPU kernels
- Leverage AMReX parallel constructs: `amrex::ParallelFor`
- Consider memory bandwidth optimization for particle operations

### Code Generation

Physics kernels are generated from Wolfram Language:
```mathematica
(* particles_geodesic.wl *)
GeodEq = (* Define geodesic equations *);
CForm[GeodEq] (* Generate C++ code *)
```

### Performance Considerations

- **Memory layout**: Structure-of-Arrays (SoA) for optimal vectorization
- **Load balancing**: Particle redistribution across AMR levels
- **Communication**: Minimal ghost cell requirements
- **I/O**: Configurable output frequencies to balance data capture and performance

### Contributing

1. Follow existing code patterns and thorn structure
2. Ensure GPU compatibility for new kernels  
3. Add appropriate tests for new functionality
4. Update documentation for API changes
5. Test across different accelerator backends (CPU/CUDA/ROCm)
