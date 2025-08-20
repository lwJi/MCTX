#include "NuParticleContainers.hxx"
#include "../wolfram/particles_derivsinline.hxx"
#include "../wolfram/particles_geodesic.hxx"
#include "Particles.hxx"

namespace NuParticleContainers {

using namespace Loop;
using namespace Particles;

std::vector<std::unique_ptr<NuParticleContainer>> g_nupcs;

NuParticleContainer::NuParticleContainer(amrex::AmrCore *amr_core)
    : Container(amr_core) {}

template <typename T>
CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline void
interp_derivs1st(T ws, dScalR &dgf_p, amrex::Array4<T const> const &gf_, int j,
                 int k, int l, int comp, VectR const &dxi) {
  dgf_p[0] += ws * fd_1_o2<0>(gf_, j, k, l, comp, dxi);
  dgf_p[1] += ws * fd_1_o2<1>(gf_, j, k, l, comp, dxi);
  dgf_p[2] += ws * fd_1_o2<2>(gf_, j, k, l, comp, dxi);
}

CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline void
gather_fields_calcrhs_at_pos(VectR &dtmom, VectR &dtpos, const VectR &mom,
                             const VectR &pos,
                             amrex::Array4<CCTK_REAL const> const &lapse_arr,
                             amrex::Array4<CCTK_REAL const> const &shift_arr,
                             amrex::Array4<CCTK_REAL const> const &met3d_arr,
                             VectR const &plo, VectR const &dxi) {

  CCTK_REAL x = (pos[0] - plo[0]) * dxi[0];
  CCTK_REAL y = (pos[1] - plo[1]) * dxi[1];
  CCTK_REAL z = (pos[2] - plo[2]) * dxi[2];

  // cell indexes
  int j = amrex::Math::floor(x);
  int k = amrex::Math::floor(y);
  int l = amrex::Math::floor(z);

  // linear interpolation weights
  CCTK_REAL xint = x - j;
  CCTK_REAL yint = y - k;
  CCTK_REAL zint = z - l;
  CCTK_REAL sx[] = {CCTK_REAL(1) - xint, xint};
  CCTK_REAL sy[] = {CCTK_REAL(1) - yint, yint};
  CCTK_REAL sz[] = {CCTK_REAL(1) - zint, zint};

  // interp metric and its derivatives
  ScalR alp_p = {0};
  VectR beta_p = {0};
  SmatR g_p = {0};

  dScalR dalp_p = {0};
  dVectR dbeta_p = {0};
  dSmatR dg_p = {0};

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

  // calculate rhs of position and momentum
  calc_rhs_geodesic(dtmom, dtpos, mom, alp_p, beta_p, g_p, dalp_p, dbeta_p,
                    dg_p);
}

void NuParticleContainer::PushAndDeposeParticles(const amrex::MultiFab &lapse,
                                                 const amrex::MultiFab &shift,
                                                 const amrex::MultiFab &met3d,
                                                 CCTK_REAL dt, const int lev) {
  const auto plo0 = Geom(0).ProbLoArray();
  const auto phi0 = Geom(0).ProbHiArray();

  const auto dxi = Geom(lev).InvCellSizeArray();
  const auto plo = Geom(lev).ProbLoArray();
  const CCTK_REAL half_dt = CCTK_REAL(0.5) * dt;

  // RK2 (midpoint): substep 1
  // Save y^n in SoA, compute k1 as temps, advance to y^{n+1/2} in-place.
  for (NuParIter pti(*this, lev); pti.isValid(); ++pti) {
    const int np = pti.numParticles();
    ParticleType *AMREX_RESTRICT pstruct = &(pti.GetArrayOfStructs()[0]);

    auto &attribs = pti.GetAttribs();
    // current (midpoint will be written in-place)
    CCTK_REAL *AMREX_RESTRICT pxp = attribs[PIdx::px].data();
    CCTK_REAL *AMREX_RESTRICT pyp = attribs[PIdx::py].data();
    CCTK_REAL *AMREX_RESTRICT pzp = attribs[PIdx::pz].data();
    // original state y^n (persist across Redistribute)
    CCTK_REAL *AMREX_RESTRICT x0a = attribs[PIdx::x0].data();
    CCTK_REAL *AMREX_RESTRICT y0a = attribs[PIdx::y0].data();
    CCTK_REAL *AMREX_RESTRICT z0a = attribs[PIdx::z0].data();
    CCTK_REAL *AMREX_RESTRICT px0a = attribs[PIdx::px0].data();
    CCTK_REAL *AMREX_RESTRICT py0a = attribs[PIdx::py0].data();
    CCTK_REAL *AMREX_RESTRICT pz0a = attribs[PIdx::pz0].data();

    auto const lapse_arr = lapse.array(pti);
    auto const shift_arr = shift.array(pti);
    auto const met3d_arr = met3d.array(pti);

    amrex::ParallelFor(np, [=] CCTK_DEVICE(int i) noexcept {
      // save y^n
      x0a[i] = pstruct[i].pos(0);
      y0a[i] = pstruct[i].pos(1);
      z0a[i] = pstruct[i].pos(2);
      px0a[i] = pxp[i];
      py0a[i] = pyp[i];
      pz0a[i] = pzp[i];

      const VectR pos0{x0a[i], y0a[i], z0a[i]};
      const VectR mom0{px0a[i], py0a[i], pz0a[i]};

      // k1 = f(y^n)  (temps only)
      VectR k1_m{}, k1_x{};
      gather_fields_calcrhs_at_pos(k1_m, k1_x, mom0, pos0, lapse_arr, shift_arr,
                                   met3d_arr, plo, dxi);

      // y^{n+1/2} = y^n + (dt/2)*k1   (write midpoint in-place)
      pstruct[i].pos(0) = x0a[i] + half_dt * k1_x[0];
      pstruct[i].pos(1) = y0a[i] + half_dt * k1_x[1];
      pstruct[i].pos(2) = z0a[i] + half_dt * k1_x[2];
      pxp[i] = px0a[i] + half_dt * k1_m[0];
      pyp[i] = py0a[i] + half_dt * k1_m[1];
      pzp[i] = pz0a[i] + half_dt * k1_m[2];
    });
  }

  // move particles that crossed tiles/ranks during substep 1
  this->Redistribute();

  // RK2 (midpoint): substep 2
  // Compute k2 at midpoint; finish with y^{n+1} = y^n + dt*k2.
  for (NuParIter pti(*this, lev); pti.isValid(); ++pti) {
    const int np = pti.numParticles();
    ParticleType *AMREX_RESTRICT pstruct = &(pti.GetArrayOfStructs()[0]);

    auto &attribs = pti.GetAttribs();
    // current (midpoint state on entry)
    CCTK_REAL *AMREX_RESTRICT pxp = attribs[PIdx::px].data();
    CCTK_REAL *AMREX_RESTRICT pyp = attribs[PIdx::py].data();
    CCTK_REAL *AMREX_RESTRICT pzp = attribs[PIdx::pz].data();
    // saved y^n
    CCTK_REAL *AMREX_RESTRICT x0a = attribs[PIdx::x0].data();
    CCTK_REAL *AMREX_RESTRICT y0a = attribs[PIdx::y0].data();
    CCTK_REAL *AMREX_RESTRICT z0a = attribs[PIdx::z0].data();
    CCTK_REAL *AMREX_RESTRICT px0a = attribs[PIdx::px0].data();
    CCTK_REAL *AMREX_RESTRICT py0a = attribs[PIdx::py0].data();
    CCTK_REAL *AMREX_RESTRICT pz0a = attribs[PIdx::pz0].data();

    auto const lapse_arr = lapse.array(pti);
    auto const shift_arr = shift.array(pti);
    auto const met3d_arr = met3d.array(pti);

    amrex::ParallelFor(np, [=] CCTK_DEVICE(int i) noexcept {
      ParticleType &p = pstruct[i];

      // midpoint state
      const VectR posh{p.pos(0), p.pos(1), p.pos(2)};
      const VectR momh{pxp[i], pyp[i], pzp[i]};

      // k2 = f(y^{n+1/2})  (temps only)
      VectR k2_m{}, k2_x{};
      gather_fields_calcrhs_at_pos(k2_m, k2_x, momh, posh, lapse_arr, shift_arr,
                                   met3d_arr, plo, dxi);

      // y^{n+1} = y^n + dt*k2  (use saved y^n)
      p.pos(0) = x0a[i] + dt * k2_x[0];
      p.pos(1) = y0a[i] + dt * k2_x[1];
      p.pos(2) = z0a[i] + dt * k2_x[2];
      pxp[i] = px0a[i] + dt * k2_m[0];
      pyp[i] = py0a[i] + dt * k2_m[1];
      pzp[i] = pz0a[i] + dt * k2_m[2];

      // Depose
      if (p.pos(0) > phi0[0] || p.pos(0) < plo0[0] || p.pos(1) > phi0[1] ||
          p.pos(1) < plo0[1] || p.pos(2) > phi0[2] || p.pos(2) < plo0[2])
        p.id() = -1;
    });
  }

  // tidy after the full step
  this->Redistribute();
}

void NuParticleContainer::OutputParticlesAscii(CCTK_ARGUMENTS) {
  DECLARE_CCTK_PARAMETERS;

  const int it = cctkGH->cctk_iteration;
  if (out_tsv_every > 0 && it % out_tsv_every == 0) {
    const std::string name =
        std::string(out_dir) + "/" + amrex::Concatenate("ptcl_asc_", it);
    amrex::Print() << "  Writing ascii file " << name << "\n";

    this->WriteAsciiFile(name);
  }
}

void NuParticleContainer::OutputParticlesPlot(CCTK_ARGUMENTS) {
  DECLARE_CCTK_PARAMETERS;

  const int it = cctkGH->cctk_iteration;
  if (out_plot_every > 0 && it % out_plot_every == 0) {
    const std::string name =
        std::string(out_dir) + "/" + amrex::Concatenate("ptcl_plt_", it);
    amrex::Print() << "  Writing plot file " << name << "\n";

    this->WritePlotFile(name, "particles");
  }
}

extern "C" void NuParticleContainers_Setup(CCTK_ARGUMENTS) {
  DECLARE_CCTK_PARAMETERS;

  for (int patch = 0; patch < CarpetX::ghext->num_patches(); ++patch) {
    const auto &patchdata = CarpetX::ghext->patchdata.at(patch);
    g_nupcs.emplace_back(
        std::make_unique<NuParticleContainer>(patchdata.amrcore.get()));
  } // for patch
}

} // namespace NuParticleContainers
