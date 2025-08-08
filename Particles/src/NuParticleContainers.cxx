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

CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline void
push_momentum(CCTK_REAL &pxp, CCTK_REAL &pyp, CCTK_REAL &pzp, CCTK_REAL pxp_rhs,
              CCTK_REAL pyp_rhs, CCTK_REAL pzp_rhs, CCTK_REAL dt) {

  pxp += pxp_rhs * dt;
  pyp += pyp_rhs * dt;
  pzp += pzp_rhs * dt;
}

CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline void
gather_fields(NuParticleContainer::ParticleType const &p, CCTK_REAL &dalpx_p,
              CCTK_REAL &dalpy_p, CCTK_REAL &dalpz_p,
              amrex::Array4<CCTK_REAL const> const &dalpx_arr,
              amrex::Array4<CCTK_REAL const> const &dalpy_arr,
              amrex::Array4<CCTK_REAL const> const &dalpz_arr,
              amrex::GpuArray<CCTK_REAL, AMREX_SPACEDIM> const &plo,
              amrex::GpuArray<CCTK_REAL, AMREX_SPACEDIM> const &dxi) {

  CCTK_REAL x = (p.pos(0) - plo[0]) * dxi[0];
  CCTK_REAL y = (p.pos(1) - plo[1]) * dxi[1];
  CCTK_REAL z = (p.pos(2) - plo[2]) * dxi[2];

  // cell indexes
  int j = std::floor(x);
  int k = std::floor(y);
  int l = std::floor(z);

  // linear interpolation weights
  CCTK_REAL xint = x - j;
  CCTK_REAL yint = y - k;
  CCTK_REAL zint = z - l;
  CCTK_REAL sx[] = {1. - xint, xint};
  CCTK_REAL sy[] = {1. - yint, yint};
  CCTK_REAL sz[] = {1. - zint, zint};

  dalpx_p = 0.;
  dalpy_p = 0.;
  dalpz_p = 0.;
  for (int ll = 0; ll <= 1; ++ll) {
    for (int kk = 0; kk <= 1; ++kk) {
      for (int jj = 0; jj <= 1; ++jj) {
        dalpx_p += sx[jj] * sy[kk] * sz[ll] * dalpx_arr(j + jj, k + kk, l + ll);
        dalpy_p += sx[jj] * sy[kk] * sz[ll] * dalpy_arr(j + jj, k + kk, l + ll);
        dalpz_p += sx[jj] * sy[kk] * sz[ll] * dalpz_arr(j + jj, k + kk, l + ll);
      }
    }
  }
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

      // push_momentum(pxp[i], pyp[i], pzp[i], ginv, Exp, Eyp, Ezp, Bxp,
      // Byp,
      //                     Bzp, q, m, dt);

      push_position(pstruct[i], pxp[i], pyp[i], pzp[i], ginv, dt);

      // deposit_current(jxarr, jyarr, jzarr, pstruct[i], pxp[i], pyp[i],
      // pzp[i],
      //                 ginv, wp[i], q, dt, plo, dxi);
    });
  }
}

void NuParticleContainer::PushParticleMomenta(const amrex::MultiFab &dalpx,
                                              const amrex::MultiFab &dalpy,
                                              const amrex::MultiFab &dalpz,
                                              CCTK_REAL dt) {

  const int lev = 0;

  const auto dxi = Geom(lev).InvCellSizeArray();
  const auto plo = Geom(lev).ProbLoArray();

  for (NuParIter pti(*this, lev); pti.isValid(); ++pti) {
    const int np = pti.numParticles();

    ParticleType const *AMREX_RESTRICT pstruct = &(pti.GetArrayOfStructs()[0]);

    auto &attribs = pti.GetAttribs();
    CCTK_REAL *AMREX_RESTRICT pxp = attribs[PIdx::px].data();
    CCTK_REAL *AMREX_RESTRICT pyp = attribs[PIdx::py].data();
    CCTK_REAL *AMREX_RESTRICT pzp = attribs[PIdx::pz].data();

    auto const dalpx_arr = dalpx.array(pti);
    auto const dalpy_arr = dalpy.array(pti);
    auto const dalpz_arr = dalpz.array(pti);

    AMREX_PARALLEL_FOR_1D(np, i, {
      CCTK_REAL pxp_rhs = 0;
      CCTK_REAL pyp_rhs = 0;
      CCTK_REAL pzp_rhs = 0;

      CCTK_REAL dalpx_p;
      CCTK_REAL dalpy_p;
      CCTK_REAL dalpz_p;

      gather_fields(pstruct[i], dalpx_p, dalpy_p, dalpz_p, dalpx_arr, dalpy_arr,
                    dalpz_arr, plo, dxi);

      push_momentum(pxp[i], pyp[i], pzp[i], pxp_rhs, pyp_rhs, pzp_rhs, dt);
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
