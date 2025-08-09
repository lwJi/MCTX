#include "NuParticleContainers.hxx"
#include "../wolfram/particles_derivsinline.hxx"
#include "../wolfram/particles_geodesic.hxx"

namespace NuParticleContainers {

using namespace Particles;

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

template <typename T>
CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline void
interp_derivs1st(T ws, std::array<T, 3> &dgf_p,
                 amrex::Array4<T const> const &gf_, int j, int k, int l,
                 int comp,
                 amrex::GpuArray<CCTK_REAL, AMREX_SPACEDIM> const &dxi) {
  dgf_p[0] += ws * fd_1_o2<0>(gf_, j, k, l, comp, dxi);
  dgf_p[1] += ws * fd_1_o2<1>(gf_, j, k, l, comp, dxi);
  dgf_p[2] += ws * fd_1_o2<2>(gf_, j, k, l, comp, dxi);
}

CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline void
gather_fields(NuParticleContainer::ParticleType const &p, CCTK_REAL &pxp_rhs,
              CCTK_REAL &pyp_rhs, CCTK_REAL &pzp_rhs,
              amrex::Array4<CCTK_REAL const> const &lapse_arr,
              amrex::Array4<CCTK_REAL const> const &shift_arr,
              amrex::Array4<CCTK_REAL const> const &met3d_arr,
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

  // interp metric and its derivatives
  CCTK_REAL alp_p = 0;
  std::array<CCTK_REAL, 3> beta_p = {0};
  std::array<CCTK_REAL, 6> g_p = {0};

  std::array<CCTK_REAL, 3> dalp_p = {0};
  std::array<std::array<CCTK_REAL, 3>, 3> dbeta_p = {0};
  std::array<std::array<CCTK_REAL, 3>, 6> dg_p = {0};

  for (int ll = 0; ll <= 1; ++ll) {
    for (int kk = 0; kk <= 1; ++kk) {
      for (int jj = 0; jj <= 1; ++jj) {
        const int j0 = j + jj;
        const int k0 = k + kk;
        const int l0 = l + ll;
        const CCTK_REAL ws = sx[jj] * sy[kk] * sz[ll];

        alp_p += ws * lapse_arr(j0, k0, l0, 0);
        interp_derivs1st(ws, dalp_p, lapse_arr, j0, k0, l0, 0, dxi);

        for (int c = 0; c < 3; ++c) {
          beta_p[c] += ws * shift_arr(j0, k0, l0, c);
          interp_derivs1st(ws, dbeta_p[c], shift_arr, j0, k0, l0, c, dxi);
        }

        for (int c = 0; c < 6; ++c) {
          g_p[c] += ws * met3d_arr(j0, k0, l0, c);
          interp_derivs1st(ws, dg_p[c], met3d_arr, j0, k0, l0, c, dxi);
        }
      }
    }
  }

  // calculate inverse of 3-metric

  pxp_rhs = 0.;
  pyp_rhs = 0.;
  pzp_rhs = 0.;
}

void NuParticleContainer::PushAndDeposeParticles(CCTK_REAL dt, const int lev) {

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

void NuParticleContainer::PushParticleMomenta(const amrex::MultiFab &lapse,
                                              const amrex::MultiFab &shift,
                                              const amrex::MultiFab &met3d,
                                              CCTK_REAL dt, const int lev) {

  const auto dxi = Geom(lev).InvCellSizeArray();
  const auto plo = Geom(lev).ProbLoArray();

  for (NuParIter pti(*this, lev); pti.isValid(); ++pti) {
    const int np = pti.numParticles();

    ParticleType const *AMREX_RESTRICT pstruct = &(pti.GetArrayOfStructs()[0]);

    auto &attribs = pti.GetAttribs();
    CCTK_REAL *AMREX_RESTRICT pxp = attribs[PIdx::px].data();
    CCTK_REAL *AMREX_RESTRICT pyp = attribs[PIdx::py].data();
    CCTK_REAL *AMREX_RESTRICT pzp = attribs[PIdx::pz].data();

    auto const lapse_arr = lapse.array(pti);
    auto const shift_arr = shift.array(pti);
    auto const met3d_arr = met3d.array(pti);

    AMREX_PARALLEL_FOR_1D(np, i, {
      CCTK_REAL pxp_rhs = 0;
      CCTK_REAL pyp_rhs = 0;
      CCTK_REAL pzp_rhs = 0;

      gather_fields(pstruct[i], pxp_rhs, pyp_rhs, pzp_rhs, lapse_arr, shift_arr,
                    met3d_arr, plo, dxi);

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
