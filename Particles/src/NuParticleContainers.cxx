#include "NuParticleContainers.hxx"

namespace NuParticleContainers {

std::vector<std::unique_ptr<NuParticleContainer>> g_nupcs;

NuParticleContainer::NuParticleContainer(amrex::AmrCore *amr_core)
    : Container(amr_core) {}

CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline void
push_position(NuParticleContainer::ParticleType &p, CCTK_REAL pxp,
              CCTK_REAL pyp, CCTK_REAL pzp, CCTK_REAL gaminv, CCTK_REAL dt) {
  p.pos(0) += pxp * gaminv * dt;
  p.pos(1) += pyp * gaminv * dt;
  p.pos(2) += pzp * gaminv * dt;
}

void NuParticleContainer::PushAndDeposeParticles(CCTK_REAL dt) {

  const int lev = 0;

  // const auto dxi = Geom(lev).InvCellSizeArray();
  const auto plo = Geom(lev).ProbLoArray();

  for (NuParIter pti(*this, lev); pti.isValid(); ++pti) {
    const int np = pti.numParticles();

    ParticleType *pstruct = &(pti.GetArrayOfStructs()[0]);

    auto &attribs = pti.GetAttribs();
    // CCTK_REAL *wp = attribs[PIdx::w].data();
    CCTK_REAL *pxp = attribs[PIdx::px].data();
    CCTK_REAL *pyp = attribs[PIdx::py].data();
    CCTK_REAL *pzp = attribs[PIdx::pz].data();

    AMREX_FOR_1D(np, i, {
      CCTK_REAL ginv = 1;

      // gather_fields(pstruct[i], Exp, Eyp, Ezp, Bxp, Byp, Bzp, Exarr, Eyarr,
      //               Ezarr, Bxarr, Byarr, Bzarr, plo, dxi);

      // push_momentum_boris(pxp[i], pyp[i], pzp[i], ginv, Exp, Eyp, Ezp, Bxp,
      // Byp,
      //                     Bzp, q, m, dt);

      push_position(pstruct[i], pxp[i], pyp[i], pzp[i], ginv, dt);

      // deposit_current(jxarr, jyarr, jzarr, pstruct[i], pxp[i], pyp[i],
      // pzp[i],
      //                 ginv, wp[i], q, dt, plo, dxi);
    });
  }
}

extern "C" void NuParticleContainers_Setup(CCTK_ARGUMENTS) {
  DECLARE_CCTK_PARAMETERS;

  for (int patch = 0; patch < CarpetX::ghext->num_patches(); ++patch) {
    const auto &restrict patchdata = CarpetX::ghext->patchdata.at(patch);
    g_nupcs.emplace_back(
        std::make_unique<NuParticleContainer>(patchdata.amrcore.get()));
  } // for patch
}

} // namespace NuParticleContainers
