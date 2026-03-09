#include <loop_device.hxx>

#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

#include <AMReX_AmrParticles.H>
#include <AMReX_Particles.H>
#include <AMReX_PlotFileUtil.H>

#include <NuParticleContainers.hxx>

#include <driver.hxx>

namespace TestNuPcsArdBH {
using namespace Loop;
using namespace NuParticleContainers;
using namespace amrex;
using namespace std;

extern "C" void TestNuPcsArdBH_InitFields(CCTK_ARGUMENTS) {
  DECLARE_CCTK_PARAMETERS;
  DECLARE_CCTK_ARGUMENTSX_TestNuPcsArdBH_InitFields;

  const auto rs = 2 * BHmass;

  grid.loop_all_device<0, 0, 0>(
      grid.nghostzones,
      [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
        const auto r = std::sqrt(p.x * p.x + p.y * p.y + p.z * p.z);

        auto f = (r > rs) ? (1 - rs / r) : 1e-10;

        const auto invf = 1.0 / f;

        alp(p.I) = std::sqrt(f);
        betax(p.I) = 0.0;
        betay(p.I) = 0.0;
        betaz(p.I) = 0.0;
        gxx(p.I) = 1.0 + (invf - 1.0) * (p.x * p.x) / (r * r);
        gxy(p.I) = 0.0 + (invf - 1.0) * (p.x * p.y) / (r * r);
        gxz(p.I) = 0.0 + (invf - 1.0) * (p.x * p.z) / (r * r);
        gyy(p.I) = 1.0 + (invf - 1.0) * (p.y * p.y) / (r * r);
        gyz(p.I) = 0.0 + (invf - 1.0) * (p.y * p.z) / (r * r);
        gzz(p.I) = 1.0 + (invf - 1.0) * (p.z * p.z) / (r * r);
      });
}

extern "C" void TestNuPcsArdBH_InitParticles(CCTK_ARGUMENTS) {
  DECLARE_CCTK_PARAMETERS;

  for (int patch = 0; patch < ghext->num_patches(); ++patch) {
    auto &pc = g_nupcs.at(patch);
    const int lev = 0;

    // Phase 1: Use AMReX InitRandom for positions + zero momenta
    // pdata.real_array_data: slots 0-2 are positions (set by InitRandom),
    //                        slots 3-5 are px, py, pz (set to 0, randomized
    //                        below)
    NuParticleContainer::ParticleInitData pdata = {
        {}, {}, {0.0, 0.0, 0.0, 0.0, 0.0, 0.0}, {}};

    // Constrain positions to the target cell near (5.05, 5.05, 0)
    const auto dx = pc->Geom(lev).CellSizeArray();
    RealBox containing_bx({5.0, 5.0, -dx[2]}, {5.0 + dx[0], 5.0 + dx[1], 0.0});

    pc->InitRandom(num_particles, random_seed, pdata, true, containing_bx);

    // Phase 2: Randomize momenta with isotropic sampling
    for (ParIterType pti(*pc, lev); pti.isValid(); ++pti) {
      auto &particle_tile = pti.GetParticleTile();
      auto ptd = particle_tile.getParticleTileData();
      const int np = pti.numParticles();

      amrex::ParallelForRNG(
          np, [=] AMREX_GPU_DEVICE(int i,
                                   amrex::RandomEngine const &engine) noexcept {
            const Real pt = 1.0;
            Real costh = Random(engine) * 2 - 1;
            Real ph = Random(engine) * (2 * M_PI);
            Real sinth = std::sqrt(amrex::max(Real(0), 1 - costh * costh));
            ptd.rdata(PIdx::px)[i] = pt * sinth * std::cos(ph);
            ptd.rdata(PIdx::py)[i] = pt * sinth * std::sin(ph);
            ptd.rdata(PIdx::pz)[i] = pt * costh;
          });
    }
  } // for patch

  // IO
  assert(ghext->num_patches() == 1);
  for (int patch = 0; patch < ghext->num_patches(); ++patch) {
    auto &pc = g_nupcs.at(patch);
    pc->OutputParticlesAscii(cctkGH);
    pc->OutputParticlesPlot(cctkGH);
    pc->OutputParticlesCheckpoint(cctkGH);
  }
}

extern "C" void TestNuPcsArdBH_PushAndDeposeParticles(CCTK_ARGUMENTS) {
  DECLARE_CCTK_PARAMETERS;
  DECLARE_CCTK_ARGUMENTS;

  constexpr int deriv_order = 2;
  for (int d = 0; d < 3; ++d)
    if (cctk_nghostzones[d] < deriv_order / 2)
      CCTK_VERROR("Need at least %d ghost zones for finite difference stencils "
                  "which appear on the rhs of geodesic evolution",
                  deriv_order / 2);

  const CCTK_REAL dt = CCTK_DELTA_TIME;

  const int tl = 0;
  const int gi_lapse = CCTK_GroupIndex("ADMBaseX::lapse");
  const int gi_shift = CCTK_GroupIndex("ADMBaseX::shift");
  const int gi_met3d = CCTK_GroupIndex("ADMBaseX::metric");
  assert(gi_lapse >= 0);
  assert(gi_shift >= 0);
  assert(gi_met3d >= 0);

  for (int patch = 0; patch < ghext->num_patches(); ++patch) {
    auto &pc = g_nupcs.at(patch);

    auto &pd = ghext->patchdata.at(patch);
    for (int lev = 0; lev < pd.leveldata.size(); ++lev) {
      const auto &ld = pd.leveldata.at(lev);
      const auto &gd_lapse = *ld.groupdata.at(gi_lapse);
      const auto &gd_shift = *ld.groupdata.at(gi_shift);
      const auto &gd_met3d = *ld.groupdata.at(gi_met3d);
      const amrex::MultiFab &lapse = *gd_lapse.mfab[tl];
      const amrex::MultiFab &shift = *gd_shift.mfab[tl];
      const amrex::MultiFab &met3d = *gd_met3d.mfab[tl];

      pc->PushAndDeposeParticles(lapse, shift, met3d, dt, lev);
    }
  } // for patch

  // IO
  assert(ghext->num_patches() == 1);

  for (int patch = 0; patch < ghext->num_patches(); ++patch) {
    auto &pc = g_nupcs.at(patch);
    pc->OutputParticlesAscii(cctkGH);
    pc->OutputParticlesPlot(cctkGH);
    pc->OutputParticlesCheckpoint(cctkGH);
  }
}

} // namespace TestNuPcsArdBH
