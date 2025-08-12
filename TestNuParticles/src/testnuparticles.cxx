#include <loop_device.hxx>

#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

#include <AMReX_AmrParticles.H>
#include <AMReX_Particles.H>
#include <AMReX_PlotFileUtil.H>

#include <NuParticleContainers.hxx>

#include "../../../CarpetX/CarpetX/src/driver.hxx"

namespace TestNuParticles {
using namespace Loop;
using namespace NuParticleContainers;
using namespace amrex;
using namespace std;

using ParticleType = NuParticleContainer::ParticleType;

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

CCTK_HOST void OutputParticles(const int it) {
  const std::string &plotfilename = amrex::Concatenate("plt", it);
  amrex::Print() << "  Writing plotfile " << plotfilename << "\n";

  for (int patch = 0; patch < ghext->num_patches(); ++patch) {
    auto &pc = g_nupcs.at(patch);

    pc->WriteAsciiFile(plotfilename);

    // pc->WritePlotFile(plotfilename, "particles");
  } // for patch
}

extern "C" void TestNuParticles_InitFields(CCTK_ARGUMENTS) {
  DECLARE_CCTK_PARAMETERS;
  DECLARE_CCTK_ARGUMENTSX_TestNuParticles_InitFields;

  grid.loop_all_device<0, 0, 0>(grid.nghostzones,
                                [=] CCTK_DEVICE(const PointDesc &p)
                                    CCTK_ATTRIBUTE_ALWAYS_INLINE {
                                      alp(p.I) = 1.0;
                                      betax(p.I) = 0.0;
                                      betay(p.I) = 0.0;
                                      betaz(p.I) = 0.0;
                                      gxx(p.I) = 1.0;
                                      gxy(p.I) = 0.0;
                                      gxz(p.I) = 0.0;
                                      gyy(p.I) = 1.0;
                                      gyz(p.I) = 0.0;
                                      gzz(p.I) = 1.0;
                                    });
}

extern "C" void TestNuParticles_InitParticles(CCTK_ARGUMENTS) {
  DECLARE_CCTK_PARAMETERS;

  const array<int, 3> nppc{4, 4, 4};

  for (int patch = 0; patch < ghext->num_patches(); ++patch) {
    auto &pc = g_nupcs.at(patch);

    // Init Particles
    const int lev = 0;

    const auto dx = pc->Geom(lev).CellSizeArray();
    const auto p_lo = pc->Geom(lev).ProbLoArray();
    const auto p_hi = pc->Geom(lev).ProbHiArray();

    const int num_ppc = nppc[0] * nppc[1] * nppc[2];

    for (MFIter mfi = pc->MakeMFIter(lev); mfi.isValid(); ++mfi) {
      const Box &tile_box = mfi.tilebox();

      const auto lo = amrex::lbound(tile_box);
      const auto hi = amrex::ubound(tile_box);

      Gpu::ManagedVector<unsigned int> counts(tile_box.numPts(), 0);
      unsigned int *pcount = counts.dataPtr();

      Gpu::ManagedVector<unsigned int> offsets(tile_box.numPts());
      unsigned int *poffset = offsets.dataPtr();

      // Counting pass: figure out exactly how many particles need to be created
      // in each grid cell
      amrex::ParallelFor(
          tile_box, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
            // Select cells
            Real xc = p_lo[0] + (i + 0.5) * dx[0];
            Real yc = p_lo[1] + (j + 0.5) * dx[1];
            Real zc = p_lo[2] + (k + 0.5) * dx[2];
            Real rc = std::sqrt(xc * xc + yc * yc + zc * zc);
            if (rc > 1.0)
              return;

            for (int i_part = 0; i_part < num_ppc; i_part++) {
              // Calculates a unique 1D index (cellid) from the 3D cell index
              // (i, j, k). This maps the 3D grid cell to a 1D memory location
              // in the counts array.
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

              // One valid particle should created here.
              pcount[cellid] += 1;
            }
          });

      // Prefix sum: this is the key that lets the second pass write data in
      // parallel without conflicts
      Gpu::exclusive_scan(counts.begin(), counts.end(), offsets.begin());

      int num_to_add =
          offsets[tile_box.numPts() - 1] + counts[tile_box.numPts() - 1];

      auto &particles = pc->GetParticles(lev);
      auto &particle_tile =
          particles[std::make_pair(mfi.index(), mfi.LocalTileIndex())];

      // Determines the current size and the required new size
      auto old_size = particle_tile.GetArrayOfStructs().size();
      auto new_size = old_size + num_to_add;

      // Crucially, this resizes the container only once, which is much more
      // efficient than adding particles one by one.
      particle_tile.resize(new_size);

      if (num_to_add == 0)
        continue;

      // Gets raw pointers to the two different ways particle data is stored for
      // performance reasons: Array of Struct (AoS) and Struct of Arrays (SoA)
      ParticleType *pstruct = particle_tile.GetArrayOfStructs()().data();
      auto arrdata = particle_tile.GetStructOfArrays().realarray();

      int procID = ParallelDescriptor::MyProc();

      // Placement pass:
      amrex::ParallelForRNG(
          tile_box,
          [=] AMREX_GPU_DEVICE(int i, int j, int k,
                               amrex::RandomEngine const &engine) noexcept {
            // Select cells
            Real xc = p_lo[0] + (i + 0.5) * dx[0];
            Real yc = p_lo[1] + (j + 0.5) * dx[1];
            Real zc = p_lo[2] + (k + 0.5) * dx[2];
            Real rc = std::sqrt(xc * xc + yc * yc + zc * zc);
            if (rc > 1.0)
              return;

            // Calculate cellid
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

            // Retrievers the starting write index (pidx) for the current cell
            // (i, j, k) from the offsets array that was calculated by the
            // exclusive_scan
            int pidx = old_size + poffset[cellid];

            for (int i_part = 0; i_part < num_ppc; i_part++) {
              Real ratio[3];
              ratio[0] = Random(engine);
              ratio[1] = Random(engine);
              ratio[2] = Random(engine);

              Real x = p_lo[0] + (i + ratio[0]) * dx[0];
              Real y = p_lo[1] + (j + ratio[1]) * dx[1];
              Real z = p_lo[2] + (k + ratio[2]) * dx[2];

              // Sampling initial momentum
              const Real pt = 1.0;
              Real costh = Random(engine) * 2 - 1;
              Real ph = Random(engine) * (2 * M_PI);
              Real sinth = std::sqrt(amrex::max(Real(0), 1 - costh * costh));
              Real cosph = std::cos(ph);
              Real sinph = std::sin(ph);

              // The core particle properties are written to the Array of
              // Structs (AoS) memory layout
              ParticleType &p = pstruct[pidx];
              p.id() = pidx + 1;
              p.cpu() = procID;
              p.pos(0) = x;
              p.pos(1) = y;
              p.pos(2) = z;

              // Write the remaining physical properties to the Struct of Arrays
              // (SoA) memory layout
              arrdata[PIdx::px][pidx] = pt * sinth * cosph;
              arrdata[PIdx::py][pidx] = pt * sinth * sinph;
              arrdata[PIdx::pz][pidx] = pt * costh;

              ++pidx;
            }
          });
    } // for mfi
  } // for patch

  // IO
  OutputParticles(cctkGH->cctk_iteration);
}

extern "C" void TestNuParticles_PushAndDeposeParticles(CCTK_ARGUMENTS) {
  DECLARE_CCTK_PARAMETERS;
  DECLARE_CCTK_ARGUMENTS;

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

  // Redistribute particles
  for (auto &pc : g_nupcs) {
    pc->Redistribute();
  }

  // IO
  OutputParticles(cctkGH->cctk_iteration);
}

} // namespace TestNuParticles
