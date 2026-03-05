#include "NuParticleContainers.hxx"

#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

#include <driver.hxx>

namespace NuParticleContainers {

std::vector<std::unique_ptr<NuParticleContainer>> g_nupcs;

extern "C" void NuParticleContainers_Setup(CCTK_ARGUMENTS) {
  DECLARE_CCTK_PARAMETERS;

  for (int patch = 0; patch < CarpetX::ghext->num_patches(); ++patch) {
    const auto &patchdata = CarpetX::ghext->patchdata.at(patch);
    g_nupcs.emplace_back(
        std::make_unique<NuParticleContainer>(patchdata.amrcore.get()));
  }
}

extern "C" void NuParticleContainers_Restart(CCTK_ARGUMENTS) {
  DECLARE_CCTK_PARAMETERS;

  const std::string dir(checkpoint_particle_dir);
  if (dir.empty())
    return;

  for (int patch = 0; patch < CarpetX::ghext->num_patches(); ++patch) {
    g_nupcs.at(patch)->RestartParticles(dir);
  }
}

} // namespace NuParticleContainers
