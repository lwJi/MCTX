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

class NuParIter : public amrex::ParIter<0, 0, PIdx::nattribs, 0> {
public:
  using amrex::ParIter<0, 0, PIdx::nattribs, 0>::ParIter;

  //    NuParIter (ContainerType& pc, int level);

  const std::array<RealVector, PIdx::nattribs> &GetAttribs() const {
    return GetStructOfArrays().GetRealData();
  }

  std::array<RealVector, PIdx::nattribs> &GetAttribs() {
    return GetStructOfArrays().GetRealData();
  }

  const RealVector &GetAttribs(int comp) const {
    return GetStructOfArrays().GetRealData(comp);
  }

  RealVector &GetAttribs(int comp) {
    return GetStructOfArrays().GetRealData(comp);
  }
};

class NuParticleContainer : public Container {

public:
  NuParticleContainer(amrex::AmrCore *amr_core);

  void PushAndDeposeParticles(CCTK_REAL dt, const int lev);

  void PushParticleMomenta(const amrex::MultiFab &dalp_arr, CCTK_REAL dt,
                           const int lev);

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
