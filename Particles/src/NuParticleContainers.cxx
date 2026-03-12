#include "NuParticleContainers.hxx"
#include "../wolfram/particles_geodesic.hxx"
#include "Particles.hxx"
#include "VertexLinearInterpolator.hxx"

#include <fstream>

namespace NuParticleContainers {

using namespace Particles;

NuParticleContainer::NuParticleContainer(amrex::AmrCore *amr_core)
    : Container(amr_core) {
  SetSoACompileTimeNames(
      {"x", "y", "z", "px", "py", "pz", "time", "num_neutrinos"},
      {"species", "cell_id"});
}

CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline void
gather_fields_calcrhs_at_pos(VectR &dtmom, VectR &dtpos, const VectR &mom,
                             const VectR &pos,
                             amrex::Array4<CCTK_REAL const> const &lapse_arr,
                             amrex::Array4<CCTK_REAL const> const &shift_arr,
                             amrex::Array4<CCTK_REAL const> const &met3d_arr,
                             VectR const &plo, VectR const &dxi) {

  const VertexLinear interp(pos, plo, dxi);

  // Interpolate metric fields and their gradients
  const ScalR alp_p = interp.gather(lapse_arr, 0);
  const dScalR dalp_p = interp.gather_deriv(lapse_arr, 0);

  VectR beta_p;
  dVectR dbeta_p;
  for (int c = 0; c < 3; ++c) {
    beta_p[c] = interp.gather(shift_arr, c);
    dbeta_p[c] = interp.gather_deriv(shift_arr, c);
  }

  SmatR g_p;
  dSmatR dg_p;
  for (int c = 0; c < 6; ++c) {
    g_p[c] = interp.gather(met3d_arr, c);
    dg_p[c] = interp.gather_deriv(met3d_arr, c);
  }

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

    // Record starting cell index before the RK2 step
    amrex::ParallelFor(np, [=] CCTK_DEVICE(int i) noexcept {
      // Compute cell indices from particle position
      const int ix = static_cast<int>(
          amrex::Math::floor((ptd.pos(0, i) - plo[0]) * dxi[0]));
      const int iy = static_cast<int>(
          amrex::Math::floor((ptd.pos(1, i) - plo[1]) * dxi[1]));
      const int iz = static_cast<int>(
          amrex::Math::floor((ptd.pos(2, i) - plo[2]) * dxi[2]));
      // Linearize: cellid = (ix * ny + iy) * nz + iz
      const int nx =
          static_cast<int>((phi0[0] - plo0[0]) * dxi[0] + CCTK_REAL(0.5));
      const int ny =
          static_cast<int>((phi0[1] - plo0[1]) * dxi[1] + CCTK_REAL(0.5));
      const int nz =
          static_cast<int>((phi0[2] - plo0[2]) * dxi[2] + CCTK_REAL(0.5));
      ptd.idata(PIdxInt::cell_id)[i] = (ix * ny + iy) * nz + iz;
    });

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

    // Advance particle time
    amrex::ParallelFor(np, [=] CCTK_DEVICE(int i) noexcept {
      ptd.rdata(PIdx::time)[i] += dt;
    });
  }
}

void NuParticleContainer::OutputParticlesAscii(const cGH *cctkGH) {
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
           << 0 << '\n'                  // NStructReal
           << 0 << '\n'                  // NStructInt
           << n_user_reals << '\n'       // user SoA reals (positions excluded)
           << PIdxInt::nattribs << '\n'; // NArrayInt
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
          using PinnedTile = amrex::ParticleTile<
              amrex::SoAParticle<PIdx::nattribs, PIdxInt::nattribs>,
              PIdx::nattribs, PIdxInt::nattribs,
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
                for (int c = 0; c < PIdxInt::nattribs; ++c) {
                  File << ptd.idata(c)[i] << ' ';
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

void NuParticleContainer::OutputParticlesPlot(const cGH *cctkGH) {
  DECLARE_CCTK_PARAMETERS;

  const int it = cctkGH->cctk_iteration;
  if (out_plot_every > 0 && it % out_plot_every == 0) {
    const std::string name =
        std::string(out_dir) + "/" + amrex::Concatenate("ptcl_plt_", it);
    amrex::Print() << "  Writing plot file " << name << "\n";

    this->WritePlotFile(name, "particles");
  }
}

void NuParticleContainer::OutputParticlesCheckpoint(const cGH *cctkGH) {
  DECLARE_CCTK_PARAMETERS;

  const int it = cctkGH->cctk_iteration;
  if (out_checkpoint_every > 0 && it % out_checkpoint_every == 0) {
    const std::string dir =
        std::string(out_dir) + "/" + amrex::Concatenate("ptcl_chk_", it);
    amrex::Print() << "  Writing particle checkpoint " << dir << "\n";

    // Names for user SoA attributes (positions excluded for pure SoA)
    const amrex::Vector<std::string> real_comp_names = {
        "px", "py", "pz", "time", "num_neutrinos"};
    const amrex::Vector<std::string> int_comp_names = {"species", "cell_id"};

    this->Checkpoint(dir, "particles", real_comp_names, int_comp_names);
  }
}

void NuParticleContainer::RestartParticles(const std::string &dir) {
  amrex::Print() << "  Restarting particles from " << dir << "\n";
  this->Restart(dir, "particles");
}

} // namespace NuParticleContainers
