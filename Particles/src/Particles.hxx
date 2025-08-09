#ifndef PARTICLES_HXX
#define PARTICLES_HXX

#include <AMReX_AmrParticles.H>
#include <AMReX_Particles.H>

namespace Particles {

using ScalR = CCTK_REAL;
using VectR = amrex::GpuArray<CCTK_REAL, 3>;
using SmatR = amrex::GpuArray<CCTK_REAL, 6>;
using dScalR = amrex::GpuArray<CCTK_REAL, 3>;
using dVectR = amrex::GpuArray<amrex::GpuArray<CCTK_REAL, 3>, 3>;
using dSmatR = amrex::GpuArray<amrex::GpuArray<CCTK_REAL, 3>, 6>;

} // namespace Particles

#endif // #ifndef PARTICLES_HXX
