
#include "NuParticleContainer.hxx"

namespace Particles {

std::vector<Container> g_nupcs;

extern "C" void NuParticleContainer_Setup(CCTK_ARGUMENTS) {
  DECLARE_CCTK_PARAMETERS;

  for (int patch = 0; patch < CarpetX::ghext->num_patches(); ++patch) {
    const auto &restrict patchdata = CarpetX::ghext->patchdata.at(patch);
    g_nupcs.push_back(Container(patchdata.amrcore.get()));
  } // for patch
}

} // namespace Particles
