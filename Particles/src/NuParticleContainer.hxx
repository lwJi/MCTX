#ifndef NUPARTICLECONTAINER_HXX
#define NUPARTICLECONTAINER_HXX

#include <AMReX_AmrParticles.H>
#include <AMReX_Particles.H>

#include "../../../CarpetX/CarpetX/src/driver.hxx"

namespace Particles {

struct PIdx {
  enum { px = 0, py, pz, nattribs };
};

using Container = amrex::AmrParticleContainer<0, 0, PIdx::nattribs, 0>;

extern std::vector<Container> g_nupcs;

} // namespace Particles

#endif
