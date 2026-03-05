#ifndef VERTEXLINEARINTERPOLATOR_HXX
#define VERTEXLINEARINTERPOLATOR_HXX

#include <cctk.h>

#include <AMReX_Array4.H>
#include <AMReX_GpuQualifiers.H>
#include <AMReX_Math.H>

#include "Particles.hxx"

namespace Particles {

// Trilinear interpolator for vertex-centered (node-centered) data.
// Modeled on amrex::ParticleInterpolator::Linear but without the +0.5
// cell-centered offset, and extended with analytic derivative support.
struct VertexLinear {
  int idx[3];                       // lower-corner cell indices
  CCTK_REAL sx[2], sy[2], sz[2];    // value weights
  CCTK_REAL dsx[2], dsy[2], dsz[2]; // derivative weights

  CCTK_DEVICE CCTK_HOST CCTK_ATTRIBUTE_ALWAYS_INLINE
  VertexLinear(const VectR &pos, const VectR &plo, const VectR &dxi) {
    CCTK_REAL x = (pos[0] - plo[0]) * dxi[0];
    CCTK_REAL y = (pos[1] - plo[1]) * dxi[1];
    CCTK_REAL z = (pos[2] - plo[2]) * dxi[2];

    idx[0] = amrex::Math::floor(x);
    idx[1] = amrex::Math::floor(y);
    idx[2] = amrex::Math::floor(z);

    CCTK_REAL xint = x - idx[0];
    CCTK_REAL yint = y - idx[1];
    CCTK_REAL zint = z - idx[2];

    // Value weights
    sx[0] = CCTK_REAL(1) - xint;
    sx[1] = xint;
    sy[0] = CCTK_REAL(1) - yint;
    sy[1] = yint;
    sz[0] = CCTK_REAL(1) - zint;
    sz[1] = zint;

    // Analytic derivative weights: d(sx)/d(pos) = d(sx)/d(xint) * dxi
    dsx[0] = -dxi[0];
    dsx[1] = dxi[0];
    dsy[0] = -dxi[1];
    dsy[1] = dxi[1];
    dsz[0] = -dxi[2];
    dsz[1] = dxi[2];
  }

  // Interpolate a single component of a grid function to the particle position.
  CCTK_DEVICE CCTK_HOST CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_REAL
  gather(amrex::Array4<CCTK_REAL const> const &arr, int comp) const {
    CCTK_REAL val = 0;
    for (int ll = 0; ll <= 1; ++ll)
      for (int kk = 0; kk <= 1; ++kk)
        for (int jj = 0; jj <= 1; ++jj)
          val += sx[jj] * sy[kk] * sz[ll] *
                 arr(idx[0] + jj, idx[1] + kk, idx[2] + ll, comp);
    return val;
  }

  // Compute the spatial gradient of a single component via analytic
  // differentiation of the trilinear interpolant.
  CCTK_DEVICE CCTK_HOST CCTK_ATTRIBUTE_ALWAYS_INLINE dScalR
  gather_deriv(amrex::Array4<CCTK_REAL const> const &arr, int comp) const {
    dScalR grad = {0, 0, 0};
    for (int ll = 0; ll <= 1; ++ll)
      for (int kk = 0; kk <= 1; ++kk)
        for (int jj = 0; jj <= 1; ++jj) {
          const CCTK_REAL f = arr(idx[0] + jj, idx[1] + kk, idx[2] + ll, comp);
          grad[0] += dsx[jj] * sy[kk] * sz[ll] * f;
          grad[1] += sx[jj] * dsy[kk] * sz[ll] * f;
          grad[2] += sx[jj] * sy[kk] * dsz[ll] * f;
        }
    return grad;
  }
};

} // namespace Particles

#endif // VERTEXLINEARINTERPOLATOR_HXX
