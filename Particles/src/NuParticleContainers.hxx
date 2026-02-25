#ifndef NUPARTICLECONTAINERS_HXX
#define NUPARTICLECONTAINERS_HXX

#include <AMReX_AmrParticles.H>
#include <AMReX_Particles.H>

#include <driver.hxx>

namespace NuParticleContainers {

struct PIdx {
  enum {
    // Positions occupy SoA slots 0..AMREX_SPACEDIM-1 (managed by AMReX)
    // User attributes start at AMREX_SPACEDIM
    px = AMREX_SPACEDIM,
    py,
    pz,
    // saved original state at t^n (persisted across the mid-step Redistribute)
    x0,
    y0,
    z0,
    px0,
    py0,
    pz0,
    nattribs
  };
};

// Pure SoA container with AmrCore tracking
using Container =
    amrex::AmrParticleContainer_impl<amrex::SoAParticle<PIdx::nattribs, 0>,
                                     PIdx::nattribs, 0>;
using ParticleTile = Container::ParticleTileType;
using ParIterType = amrex::ParIterSoA<PIdx::nattribs, 0>;

class NuParticleContainer : public Container {

public:
  NuParticleContainer(amrex::AmrCore *amr_core);

  void PushAndDeposeParticles(const amrex::MultiFab &lapse,
                              const amrex::MultiFab &shift,
                              const amrex::MultiFab &met3d, CCTK_REAL dt,
                              const int lev);

  void OutputParticlesAscii(CCTK_ARGUMENTS);

  void OutputParticlesPlot(CCTK_ARGUMENTS);
};

extern std::vector<std::unique_ptr<NuParticleContainer>> g_nupcs;

} // namespace NuParticleContainers

#endif // #ifndef NUPARTICLECONTAINERS_HXX
