# MCTX

Monte Carlo Transport framework built on CarpetX/AMReX within the Einstein Toolkit (Cactus) for computational relativistic particle transport.

## Tech Stack

- C++ with AMReX for GPU-accelerated particle containers
- Cactus computational toolkit with CarpetX for adaptive mesh refinement
- Mathematica/Wolfram Language for symbolic derivations (`Particles/wolfram/`)

## Project Structure

- `Particles/` — Main thorn: neutrino particle containers, geodesic integration, I/O
- `TestNuParticles/` — Test thorn for neutrino particle simulations
- `TestNuPcsArdBH/` — Test thorn for neutrino particles around black holes
- `scripts/` — CI support files only

Each thorn follows the Cactus format:
- `configuration.ccl`, `interface.ccl`, `param.ccl`, `schedule.ccl` — Thorn metadata
- `src/` — C++ implementation
- `test/` — Parameter files and reference output for the test suite

## Build & Test

```bash
./agent_scripts/build_and_test.sh
```

## Architecture

- Particles are managed via AMReX `AmrParticleContainer` with SoA attributes for momentum (`NuParticleContainers.hxx`)
- `NuParticleContainer` handles particle push (geodesic integration), deposition, and I/O
- Symbolic expressions in `Particles/wolfram/` generate C++ headers (`.hxx`) — do not edit generated headers directly

## Key Patterns

- GPU compatibility: all particle kernels must work with AMReX GPU abstractions
- New functionality: add parameters in `.ccl` files, implement in `src/`, add tests in a test thorn
- Test thorns use Cactus test suite infrastructure — parameter files in `test/` with reference output
