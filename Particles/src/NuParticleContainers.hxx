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

extern std::vector<Container> g_nupcs;

} // namespace NuParticleContainers

#endif // #ifndef NUPARTICLECONTAINERS_HXX
