#ifndef NUPARTICLECONTAINERS_HXX
#define NUPARTICLECONTAINERS_HXX

#include <cctk.h>

#include <AMReX_AmrParticles.H>
#include <AMReX_Particles.H>

namespace NuParticleContainers {

struct PIdx {
  enum {
    // Positions occupy SoA slots 0..AMREX_SPACEDIM-1 (managed by AMReX)
    // User attributes start at AMREX_SPACEDIM
    px = AMREX_SPACEDIM,
    py,
    pz,
    nattribs // = 6
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

  void OutputParticlesAscii(const cGH *cctkGH);

  void OutputParticlesPlot(const cGH *cctkGH);

  void OutputParticlesCheckpoint(const cGH *cctkGH);

  void RestartParticles(const std::string &dir);
};

extern std::vector<std::unique_ptr<NuParticleContainer>> g_nupcs;

} // namespace NuParticleContainers

#endif // #ifndef NUPARTICLECONTAINERS_HXX
