# MCTX

Monte Carlo Transport framework built on CarpetX/AMReX within the Einstein Toolkit (Cactus) for computational relativistic particle transport.

## Tech Stack

- C++ with AMReX for GPU-accelerated particle containers
- Cactus computational toolkit with CarpetX for adaptive mesh refinement
- SimFactory for build configuration and job management
- Mathematica/Wolfram Language for symbolic derivations (`Particles/wolfram/`)

## Project Structure

- `Particles/` — Main thorn: neutrino particle containers, geodesic integration, I/O
- `TestNuParticles/` — Test thorn for neutrino particle simulations
- `TestNuPcsArdBH/` — Test thorn for neutrino particles around black holes
- `scripts/` — Build/test scripts, SimFactory configs, CI support files

Each thorn follows the Cactus format:
- `configuration.ccl`, `interface.ccl`, `param.ccl`, `schedule.ccl` — Thorn metadata
- `src/` — C++ implementation
- `test/` — Parameter files and reference output for the test suite

## Development

```
./scripts/download.sh   # Download Cactus framework and all dependencies
./scripts/build.sh      # Build (requires ACCELERATOR, REAL_PRECISION, MODE env vars)
./scripts/test.sh       # Run test suite (1-proc and 2-proc configurations)
```

Environment variables for build: `ACCELERATOR` (cpu/cuda/rocm), `REAL_PRECISION` (real64), `MODE` (debug/optimize).

CI runs via GitHub Actions (`.github/workflows/ci.yml`) in `einsteintoolkit/carpetx` Docker containers.

## Architecture

- Particles are managed via AMReX `AmrParticleContainer` with SoA attributes for momentum (`NuParticleContainers.hxx`)
- `NuParticleContainer` handles particle push (geodesic integration), deposition, and I/O
- Symbolic expressions in `Particles/wolfram/` generate C++ headers (`.hxx`) — do not edit generated headers directly
- The thornlist `scripts/mctx.th` defines all external dependencies fetched by `download.sh`

## Key Patterns

- GPU compatibility: all particle kernels must work with AMReX GPU abstractions
- New functionality: add parameters in `.ccl` files, implement in `src/`, add tests in a test thorn
- Test thorns use Cactus test suite infrastructure — parameter files in `test/` with reference output
