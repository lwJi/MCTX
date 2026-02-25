#include "NuParticleContainers.hxx"
#include "../wolfram/particles_derivsinline.hxx"
#include "../wolfram/particles_geodesic.hxx"
#include "Particles.hxx"

#include <fstream>

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
  for (ParIterType pti(*this, lev); pti.isValid(); ++pti) {
    const int np = pti.numParticles();

    auto &soa = pti.GetStructOfArrays();
    // positions
    auto *AMREX_RESTRICT xp = soa.GetRealData(0).data();
    auto *AMREX_RESTRICT yp = soa.GetRealData(1).data();
    auto *AMREX_RESTRICT zp = soa.GetRealData(2).data();
    // momenta
    auto *AMREX_RESTRICT pxp = soa.GetRealData(PIdx::px).data();
    auto *AMREX_RESTRICT pyp = soa.GetRealData(PIdx::py).data();
    auto *AMREX_RESTRICT pzp = soa.GetRealData(PIdx::pz).data();
    // saved state y^n
    auto *AMREX_RESTRICT x0a = soa.GetRealData(PIdx::x0).data();
    auto *AMREX_RESTRICT y0a = soa.GetRealData(PIdx::y0).data();
    auto *AMREX_RESTRICT z0a = soa.GetRealData(PIdx::z0).data();
    auto *AMREX_RESTRICT px0a = soa.GetRealData(PIdx::px0).data();
    auto *AMREX_RESTRICT py0a = soa.GetRealData(PIdx::py0).data();
    auto *AMREX_RESTRICT pz0a = soa.GetRealData(PIdx::pz0).data();

    auto const lapse_arr = lapse.array(pti);
    auto const shift_arr = shift.array(pti);
    auto const met3d_arr = met3d.array(pti);

    amrex::ParallelFor(np, [=] CCTK_DEVICE(int i) noexcept {
      // save y^n
      x0a[i] = xp[i];
      y0a[i] = yp[i];
      z0a[i] = zp[i];
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
      xp[i] = x0a[i] + half_dt * k1_x[0];
      yp[i] = y0a[i] + half_dt * k1_x[1];
      zp[i] = z0a[i] + half_dt * k1_x[2];
      pxp[i] = px0a[i] + half_dt * k1_m[0];
      pyp[i] = py0a[i] + half_dt * k1_m[1];
      pzp[i] = pz0a[i] + half_dt * k1_m[2];
    });
  }

  // move particles that crossed tiles/ranks during substep 1
  this->Redistribute();

  // RK2 (midpoint): substep 2
  // Compute k2 at midpoint; finish with y^{n+1} = y^n + dt*k2.
  for (ParIterType pti(*this, lev); pti.isValid(); ++pti) {
    const int np = pti.numParticles();

    auto &soa = pti.GetStructOfArrays();
    auto *AMREX_RESTRICT xp = soa.GetRealData(0).data();
    auto *AMREX_RESTRICT yp = soa.GetRealData(1).data();
    auto *AMREX_RESTRICT zp = soa.GetRealData(2).data();
    auto *AMREX_RESTRICT pxp = soa.GetRealData(PIdx::px).data();
    auto *AMREX_RESTRICT pyp = soa.GetRealData(PIdx::py).data();
    auto *AMREX_RESTRICT pzp = soa.GetRealData(PIdx::pz).data();
    auto *AMREX_RESTRICT x0a = soa.GetRealData(PIdx::x0).data();
    auto *AMREX_RESTRICT y0a = soa.GetRealData(PIdx::y0).data();
    auto *AMREX_RESTRICT z0a = soa.GetRealData(PIdx::z0).data();
    auto *AMREX_RESTRICT px0a = soa.GetRealData(PIdx::px0).data();
    auto *AMREX_RESTRICT py0a = soa.GetRealData(PIdx::py0).data();
    auto *AMREX_RESTRICT pz0a = soa.GetRealData(PIdx::pz0).data();
    auto *AMREX_RESTRICT idcpu_arr = soa.GetIdCPUData().data();

    auto const lapse_arr = lapse.array(pti);
    auto const shift_arr = shift.array(pti);
    auto const met3d_arr = met3d.array(pti);

    amrex::ParallelFor(np, [=] CCTK_DEVICE(int i) noexcept {
      // midpoint state
      const VectR posh{xp[i], yp[i], zp[i]};
      const VectR momh{pxp[i], pyp[i], pzp[i]};

      // k2 = f(y^{n+1/2})  (temps only)
      VectR k2_m{}, k2_x{};
      gather_fields_calcrhs_at_pos(k2_m, k2_x, momh, posh, lapse_arr, shift_arr,
                                   met3d_arr, plo, dxi);

      // y^{n+1} = y^n + dt*k2  (use saved y^n)
      xp[i] = x0a[i] + dt * k2_x[0];
      yp[i] = y0a[i] + dt * k2_x[1];
      zp[i] = z0a[i] + dt * k2_x[2];
      pxp[i] = px0a[i] + dt * k2_m[0];
      pyp[i] = py0a[i] + dt * k2_m[1];
      pzp[i] = pz0a[i] + dt * k2_m[2];

      // Depose
      if (xp[i] > phi0[0] || xp[i] < plo0[0] || yp[i] > phi0[1] ||
          yp[i] < plo0[1] || zp[i] > phi0[2] || zp[i] < plo0[2])
        amrex::ParticleIDWrapper<uint64_t>{idcpu_arr[i]} = -1;
    });
  }

  // tidy after the full step
  this->Redistribute();
}

void NuParticleContainer::OutputParticlesAscii(CCTK_ARGUMENTS) {
  DECLARE_CCTK_PARAMETERS;

  const int it = cctkGH->cctk_iteration;
  if (out_tsv_every > 0 && it % out_tsv_every == 0) {
    const std::string name = std::string(out_dir) + "/" +
                             amrex::Concatenate("ptcl_asc_", it) + ".tsv";
    amrex::Print() << "  Writing ascii file " << name << "\n";

    // Pure SoA replacement for WriteAsciiFile (which requires AoS).
    // Produces identical output format: positions, id, cpu, then user SoA
    // attributes (momenta + saved state), excluding the position SoA slots.
    const auto nparticles = this->TotalNumberOfParticles(true, false);
    constexpr int n_user_reals = PIdx::nattribs - AMREX_SPACEDIM;

    // Header (I/O processor only)
    if (amrex::ParallelDescriptor::IOProcessor()) {
      std::ofstream File(name, std::ios::out | std::ios::trunc);
      if (!File.good())
        amrex::FileOpenFailed(name);
      File << nparticles << '\n'
           << 0 << '\n'            // NStructReal
           << 0 << '\n'            // NStructInt
           << n_user_reals << '\n' // user SoA reals (positions excluded)
           << 0 << '\n';           // NArrayInt
      File.close();
    }

    amrex::ParallelDescriptor::Barrier();

    // Each processor appends its particles sequentially
    const int MyProc = amrex::ParallelDescriptor::MyProc();
    for (int proc = 0; proc < amrex::ParallelDescriptor::NProcs(); ++proc) {
      if (MyProc == proc) {
        std::ofstream File(name, std::ios::out | std::ios::app);
        File.precision(15);
        if (!File.good())
          amrex::FileOpenFailed(name);

        for (int lev = 0; lev <= finestLevel(); ++lev) {
          using PinnedTile =
              amrex::ParticleTile<amrex::SoAParticle<PIdx::nattribs, 0>,
                                  PIdx::nattribs, 0,
                                  amrex::PolymorphicArenaAllocator>;

          for (const auto &kv : GetParticles(lev)) {
            PinnedTile pinned;
            pinned.define(NumRuntimeRealComps(), NumRuntimeIntComps(), nullptr,
                          nullptr, amrex::The_Pinned_Arena());
            pinned.resize(kv.second.numParticles());
            amrex::copyParticles(pinned, kv.second);

            auto &soa = pinned.GetStructOfArrays();
            const int np = pinned.numParticles();
            auto *xp = soa.GetRealData(0).data();
            auto *yp = soa.GetRealData(1).data();
            auto *zp = soa.GetRealData(2).data();
            auto *idcpu = soa.GetIdCPUData().data();

            for (int i = 0; i < np; ++i) {
              if (amrex::ParticleIDWrapper<>(idcpu[i]).is_valid()) {
                File << xp[i] << ' ' << yp[i] << ' ' << zp[i] << ' ';
                File << int(amrex::ParticleIDWrapper<>(idcpu[i])) << ' ';
                File << int(amrex::ParticleCPUWrapper(idcpu[i])) << ' ';
                for (int c = AMREX_SPACEDIM; c < PIdx::nattribs; ++c) {
                  File << soa.GetRealData(c)[i] << ' ';
                }
                File << '\n';
              }
            }
          }
        }

        File.flush();
        File.close();
        if (!File.good())
          amrex::Abort("OutputParticlesAscii: problem writing file");
      }
      amrex::ParallelDescriptor::Barrier();
    }
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
