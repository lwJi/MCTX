#ifndef NUPARTICLECONTAINERS_HXX
#define NUPARTICLECONTAINERS_HXX

#include <AMReX_AmrParticles.H>
#include <AMReX_Particles.H>

#include "../../../CarpetX/CarpetX/src/driver.hxx"

namespace NuParticleContainers {

struct PIdx {
  enum { px = 0, py, pz, nattribs };
};

using Container = amrex::AmrParticleContainer<0, 0, PIdx::nattribs, 0>;
using ParticleTile = Container::ParticleTileType;

class NuParticleContainer : public Container {

public:
  NuParticleContainer(amrex::AmrCore *amr_core);

  // void InitParticles(const amrex::IntVect &a_num_particles_per_cell,
  //                    const amrex::Real a_thermal_momentum_std,
  //                    const amrex::Real a_thermal_momentum_mean,
  //                    const amrex::Real a_density,
  //                    const amrex::RealBox &a_bounds, const int a_problem);

  // void
  // PushAndDeposeParticles(const amrex::MultiFab &Ex, const amrex::MultiFab
  // &Ey,
  //                        const amrex::MultiFab &Ez, const amrex::MultiFab
  //                        &Bx, const amrex::MultiFab &By, const
  //                        amrex::MultiFab &Bz, amrex::MultiFab &jx,
  //                        amrex::MultiFab &jy, amrex::MultiFab &jz,
  //                        amrex::Real dt);

  // void PushParticleMomenta(const amrex::MultiFab &Ex, const amrex::MultiFab
  // &Ey,
  //                          const amrex::MultiFab &Ez, const amrex::MultiFab
  //                          &Bx, const amrex::MultiFab &By, const
  //                          amrex::MultiFab &Bz, amrex::Real dt);

  // void RedistributeLocal() {
  //   const int lev_min = 0;
  //   const int lev_max = 0;
  //   const int nGrow = 0;
  //   const int local = 1;
  //   Redistribute(lev_min, lev_max, nGrow, local);
  // }

  // protected:
  //   int m_species_id;
};

extern std::vector<std::unique_ptr<NuParticleContainer>> g_nupcs;

} // namespace NuParticleContainers

#endif // #ifndef NUPARTICLECONTAINERS_HXX
