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
    time,
    num_neutrinos,
    pt,      // p_t = p_i beta^i - alpha sqrt(gamma^{ij} p_i p_j)
    nattribs // = 9
  };
};

struct PIdxInt {
  enum {
    species = 0,
    cell_id,
    nattribs // = 2
  };
};

// Pure SoA container with AmrCore tracking
using Container = amrex::AmrParticleContainer_impl<
    amrex::SoAParticle<PIdx::nattribs, PIdxInt::nattribs>, PIdx::nattribs,
    PIdxInt::nattribs>;
using ParticleTile = Container::ParticleTileType;
using ParIterType = amrex::ParIterSoA<PIdx::nattribs, PIdxInt::nattribs>;

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
