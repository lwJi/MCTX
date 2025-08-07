#include "NuParticleContainers.hxx"

namespace NuParticleContainers {

std::vector<Container> g_nupcs;

extern "C" void NuParticleContainers_Setup(CCTK_ARGUMENTS) {
  DECLARE_CCTK_PARAMETERS;

  for (int patch = 0; patch < CarpetX::ghext->num_patches(); ++patch) {
    const auto &restrict patchdata = CarpetX::ghext->patchdata.at(patch);
    g_nupcs.push_back(Container(patchdata.amrcore.get()));
  } // for patch
}

} // namespace NuParticleContainers
