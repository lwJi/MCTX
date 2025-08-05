#include <loop_device.hxx>

#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

#include <AMReX_AmrParticles.H>
#include <AMReX_Particles.H>

#include "../../../CarpetX/CarpetX/src/driver.hxx"

namespace TestPaticles {
using namespace amrex;
using namespace std;

struct PIdx {
  enum { ux = 0, uy, uz, nattribs };
};

CCTK_HOST CCTK_DEVICE void
get_position_unit_cell(Real *r, const array<int, 3> &nppc, int i_part) {
  int nx = nppc[0];
  int ny = nppc[1];
  int nz = nppc[2];

  int ix_part = i_part / (ny * nz);
  int iy_part = (i_part % (ny * nz)) % ny;
  int iz_part = (i_part % (ny * nz)) / ny;

  r[0] = (0.5 + ix_part) / nx;
  r[1] = (0.5 + iy_part) / ny;
  r[2] = (0.5 + iz_part) / nz;
}

extern "C" void TestParticles_Init(CCTK_ARGUMENTS) {
  DECLARE_CCTK_PARAMETERS;

  const array<int, 3> nppc{4, 4, 4};

  using Container = amrex::AmrParticleContainer<0, 0, PIdx::nattribs, 0>;
  using ParticleTile = Container::ParticleTileType;
  std::vector<Container> containers(ghext->num_patches());

  for (int patch = 0; patch < ghext->num_patches(); ++patch) {
    const auto &restrict patchdata = ghext->patchdata.at(patch);
    containers.at(patch) = Container(patchdata.amrcore.get());
    const auto &pc = containers.at(patch);

    // Init Particles
    const int lev = 0;

    const amrex::Geometry &geom = ghext->patchdata.at(patch).amrcore->Geom(lev);
    const auto dx = geom.CellSizeArray();
    const auto plo = geom.ProbLoArray();

    const int num_ppc = nppc[0] * nppc[1] * nppc[2];
    const Real scale_fac = dx[0] * dx[1] * dx[2] / num_ppc;

    for (MFIter mfi = pc.MakeMFIter(lev); mfi.isValid(); ++mfi) {
      const Box &tile_box = mfi.tilebox();

      const auto lo = amrex::lbound(tile_box);
      const auto hi = amrex::ubound(tile_box);

      Gpu::ManagedVector<unsigned int> counts(tile_box.numPts(), 0);
      unsigned int *pcount = counts.dataPtr();

      Gpu::ManagedVector<unsigned int> offsets(tile_box.numPts());
      unsigned int *poffset = offsets.dataPtr();

      amrex::ParallelFor(
          tile_box, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
            for (int i_part = 0; i_part < num_ppc; i_part++) {
              Real r[3];

              get_position_unit_cell(r, nppc, i_part);

              Real x = plo[0] + (i + r[0]) * dx[0];
              Real y = plo[1] + (j + r[1]) * dx[1];
              Real z = plo[2] + (k + r[2]) * dx[2];

              if (x >= a_bounds.hi(0) || x < a_bounds.lo(0) ||
                  y >= a_bounds.hi(1) || y < a_bounds.lo(1) ||
                  z >= a_bounds.hi(2) || z < a_bounds.lo(2))
                continue;

              int ix = i - lo.x;
              int iy = j - lo.y;
              int iz = k - lo.z;
              int nx = hi.x - lo.x + 1;
              int ny = hi.y - lo.y + 1;
              int nz = hi.z - lo.z + 1;
              unsigned int uix = amrex::min(nx - 1, amrex::max(0, ix));
              unsigned int uiy = amrex::min(ny - 1, amrex::max(0, iy));
              unsigned int uiz = amrex::min(nz - 1, amrex::max(0, iz));
              unsigned int cellid = (uix * ny + uiy) * nz + uiz;
              pcount[cellid] += 1;
            }
          });

      Gpu::exclusive_scan(counts.begin(), counts.end(), offsets.begin());

      int num_to_add =
          offsets[tile_box.numPts() - 1] + counts[tile_box.numPts() - 1];

      auto &particles = GetParticles(lev);
      auto &particle_tile =
          particles[std::make_pair(mfi.index(), mfi.LocalTileIndex())];

      auto old_size = particle_tile.GetArrayOfStructs().size();
      auto new_size = old_size + num_to_add;
      particle_tile.resize(new_size);

      if (num_to_add == 0)
        continue;

      ParticleType *pstruct = particle_tile.GetArrayOfStructs()().data();

      auto arrdata = particle_tile.GetStructOfArrays().realarray();

      int procID = ParallelDescriptor::MyProc();

      amrex::ParallelForRNG(
          tile_box,
          [=] AMREX_GPU_DEVICE(int i, int j, int k,
                               amrex::RandomEngine const &engine) noexcept {
            int ix = i - lo.x;
            int iy = j - lo.y;
            int iz = k - lo.z;
            int nx = hi.x - lo.x + 1;
            int ny = hi.y - lo.y + 1;
            int nz = hi.z - lo.z + 1;
            unsigned int uix = amrex::min(nx - 1, amrex::max(0, ix));
            unsigned int uiy = amrex::min(ny - 1, amrex::max(0, iy));
            unsigned int uiz = amrex::min(nz - 1, amrex::max(0, iz));
            unsigned int cellid = (uix * ny + uiy) * nz + uiz;

            int pidx = poffset[cellid] - poffset[0];

            for (int i_part = 0; i_part < num_ppc; i_part++) {
              Real r[3];
              Real u[3];

              get_position_unit_cell(r, nppc, i_part);

              Real x = plo[0] + (i + r[0]) * dx[0];
              Real y = plo[1] + (j + r[1]) * dx[1];
              Real z = plo[2] + (k + r[2]) * dx[2];

              u[0] = 0.01;
              u[1] = 0.0;
              u[2] = 0.0;

              if (x >= a_bounds.hi(0) || x < a_bounds.lo(0) ||
                  y >= a_bounds.hi(1) || y < a_bounds.lo(1) ||
                  z >= a_bounds.hi(2) || z < a_bounds.lo(2))
                continue;

              ParticleType &p = pstruct[pidx];
              p.id() = pidx + 1;
              p.cpu() = procID;
              p.pos(0) = x;
              p.pos(1) = y;
              p.pos(2) = z;

              arrdata[PIdx::ux][pidx] = u[0];
              arrdata[PIdx::uy][pidx] = u[1];
              arrdata[PIdx::uz][pidx] = u[2];

              ++pidx;
            }
          });
    } // for mfi
  } // for patch
}

// extern "C" void TestParticles_Update(CCTK_ARGUMENTS) {
//   DECLARE_CCTK_PARAMETERS;
//   DECLARE_CCTK_ARGUMENTSX_TestParticles_Update;
// }

} // namespace TestPaticles
