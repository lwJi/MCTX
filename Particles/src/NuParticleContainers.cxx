#include "NuParticleContainers.hxx"
#include "../wolfram/particles_derivsinline.hxx"
#include "../wolfram/particles_geodesic.hxx"
#include "Particles.hxx"

#include <driver.hxx>

#include <fstream>

namespace NuParticleContainers {

using namespace Loop;
using namespace Particles;

std::vector<std::unique_ptr<NuParticleContainer>> g_nupcs;

NuParticleContainer::NuParticleContainer(amrex::AmrCore *amr_core)
    : Container(amr_core) {
  SetSoACompileTimeNames({"x", "y", "z", "px", "py", "pz"}, {});
}

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

  // RK2 (midpoint): both substeps in a single tile loop
  for (ParIterType pti(*this, lev); pti.isValid(); ++pti) {
    const int np = pti.numParticles();
    auto ptd = pti.GetParticleTile().getParticleTileData();

    auto const lapse_arr = lapse.array(pti);
    auto const shift_arr = shift.array(pti);
    auto const met3d_arr = met3d.array(pti);

    // Tile-local scratch: packed [x0|y0|z0|px0|py0|pz0], each segment np long
    amrex::Gpu::DeviceVector<CCTK_REAL> scratch(6 * np);
    auto *const x0 = scratch.data();
    auto *const y0 = scratch.data() + np;
    auto *const z0 = scratch.data() + 2 * np;
    auto *const px0 = scratch.data() + 3 * np;
    auto *const py0 = scratch.data() + 4 * np;
    auto *const pz0 = scratch.data() + 5 * np;

    // Substep 1: save y^n, compute k1, advance to y^{n+1/2}
    amrex::ParallelFor(np, [=] CCTK_DEVICE(int i) noexcept {
      // save y^n into tile-local scratch
      x0[i] = ptd.pos(0, i);
      y0[i] = ptd.pos(1, i);
      z0[i] = ptd.pos(2, i);
      px0[i] = ptd.rdata(PIdx::px)[i];
      py0[i] = ptd.rdata(PIdx::py)[i];
      pz0[i] = ptd.rdata(PIdx::pz)[i];

      const VectR pos0{ptd.pos(0, i), ptd.pos(1, i), ptd.pos(2, i)};
      const VectR mom0{ptd.rdata(PIdx::px)[i], ptd.rdata(PIdx::py)[i],
                       ptd.rdata(PIdx::pz)[i]};

      // k1 = f(y^n)
      VectR k1_m{}, k1_x{};
      gather_fields_calcrhs_at_pos(k1_m, k1_x, mom0, pos0, lapse_arr, shift_arr,
                                   met3d_arr, plo, dxi);

      // y^{n+1/2} = y^n + (dt/2)*k1
      ptd.pos(0, i) = x0[i] + half_dt * k1_x[0];
      ptd.pos(1, i) = y0[i] + half_dt * k1_x[1];
      ptd.pos(2, i) = z0[i] + half_dt * k1_x[2];
      ptd.rdata(PIdx::px)[i] = px0[i] + half_dt * k1_m[0];
      ptd.rdata(PIdx::py)[i] = py0[i] + half_dt * k1_m[1];
      ptd.rdata(PIdx::pz)[i] = pz0[i] + half_dt * k1_m[2];
    });

    // Substep 2: compute k2 at midpoint, write y^{n+1} = y^n + dt*k2
    amrex::ParallelFor(np, [=] CCTK_DEVICE(int i) noexcept {
      const VectR posh{ptd.pos(0, i), ptd.pos(1, i), ptd.pos(2, i)};
      const VectR momh{ptd.rdata(PIdx::px)[i], ptd.rdata(PIdx::py)[i],
                       ptd.rdata(PIdx::pz)[i]};

      // k2 = f(y^{n+1/2})
      VectR k2_m{}, k2_x{};
      gather_fields_calcrhs_at_pos(k2_m, k2_x, momh, posh, lapse_arr, shift_arr,
                                   met3d_arr, plo, dxi);

      // y^{n+1} = y^n + dt*k2
      ptd.pos(0, i) = x0[i] + dt * k2_x[0];
      ptd.pos(1, i) = y0[i] + dt * k2_x[1];
      ptd.pos(2, i) = z0[i] + dt * k2_x[2];
      ptd.rdata(PIdx::px)[i] = px0[i] + dt * k2_m[0];
      ptd.rdata(PIdx::py)[i] = py0[i] + dt * k2_m[1];
      ptd.rdata(PIdx::pz)[i] = pz0[i] + dt * k2_m[2];

      // Depose
      if (ptd.pos(0, i) > phi0[0] || ptd.pos(0, i) < plo0[0] ||
          ptd.pos(1, i) > phi0[1] || ptd.pos(1, i) < plo0[1] ||
          ptd.pos(2, i) > phi0[2] || ptd.pos(2, i) < plo0[2])
        ptd.id(i) = -1;
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

            const int np = pinned.numParticles();
            auto ptd = pinned.getConstParticleTileData();

            for (int i = 0; i < np; ++i) {
              if (ptd.id(i).is_valid()) {
                File << ptd.pos(0, i) << ' ' << ptd.pos(1, i) << ' '
                     << ptd.pos(2, i) << ' ';
                File << int(ptd.id(i)) << ' ';
                File << int(ptd.cpu(i)) << ' ';
                for (int c = AMREX_SPACEDIM; c < PIdx::nattribs; ++c) {
                  File << ptd.rdata(c)[i] << ' ';
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

void NuParticleContainer::OutputParticlesCheckpoint(CCTK_ARGUMENTS) {
  DECLARE_CCTK_PARAMETERS;

  const int it = cctkGH->cctk_iteration;
  if (out_checkpoint_every > 0 && it % out_checkpoint_every == 0) {
    const std::string dir =
        std::string(out_dir) + "/" + amrex::Concatenate("ptcl_chk_", it);
    amrex::Print() << "  Writing particle checkpoint " << dir << "\n";

    // Names for user SoA attributes (positions excluded for pure SoA)
    const amrex::Vector<std::string> real_comp_names = {"px", "py", "pz"};
    const amrex::Vector<std::string> int_comp_names = {};

    this->Checkpoint(dir, "particles", real_comp_names, int_comp_names);
  }
}

void NuParticleContainer::RestartParticles(const std::string &dir) {
  amrex::Print() << "  Restarting particles from " << dir << "\n";
  this->Restart(dir, "particles");
}

extern "C" void NuParticleContainers_Setup(CCTK_ARGUMENTS) {
  DECLARE_CCTK_PARAMETERS;

  for (int patch = 0; patch < CarpetX::ghext->num_patches(); ++patch) {
    const auto &patchdata = CarpetX::ghext->patchdata.at(patch);
    g_nupcs.emplace_back(
        std::make_unique<NuParticleContainer>(patchdata.amrcore.get()));
  } // for patch
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
