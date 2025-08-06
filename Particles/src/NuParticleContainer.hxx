#ifndef NUPARTICLECONTAINER_HXX
#define NUPARTICLECONTAINER_HXX

#include <AMReX_Particles.H>

namespace Particles {

struct PIdx {
  enum { px = 0, py, pz, nattribs };
};

using Container = amrex::AmrParticleContainer<0, 0, PIdx::nattribs, 0>;

class NuParticleContainer {

public:
  NuParticleContainer();

protected:
  int m_species_id;
  std::vector<Container> m_containers;
};

} // namespace Particles

#endif
