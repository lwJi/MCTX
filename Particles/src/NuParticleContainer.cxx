
#include "NuParticleContainer.hxx"

namespace Particles {

NuParticleContainer::NuParticleContainer() {
  for (int patch = 0; patch < CarpetX::ghext->num_patches(); ++patch) {
    const auto &restrict patchdata = CarpetX::ghext->patchdata.at(patch);
    m_containers.push_back(Container(patchdata.amrcore.get()));
  } // for patch
}

} // namespace Particles
