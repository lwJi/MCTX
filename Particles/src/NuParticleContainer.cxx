
#include "NuParticleContainer.hxx"

namespace Particles {

void NuParticleContainer::NuParticleContainer() {
  for (int patch = 0; patch < ghext->num_patches(); ++patch) {
    const auto &restrict patchdata = ghext->patchdata.at(patch);
    m_containers.push_back(Container(patchdata.amrcore.get()))
  } // for patch
}

} // namespace Particles
