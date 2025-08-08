/* particles_derivsinline.hxx */
/* Produced with Generato */

#ifndef PARTICLES_DERIVSINLINE_HXX
#define PARTICLES_DERIVSINLINE_HXX

#include <loop_device.hxx>

#include <array>
#include <cmath>

#include "particles_powerinline.hxx"

namespace Particles {
using namespace Loop;

template <int D, typename T>
CCTK_DEVICE CCTK_HOST CCTK_ATTRIBUTE_ALWAYS_INLINE inline T
fd_1_o2(amrex::Array4<T const> const &gf, int i, int j, int k, int comp, amrex::GpuArray<T, 3> const &dxi) {
  const T m1 = gf(i + (D == 0 ? -1 : 0), j + (D == 1 ? -1 : 0), k + (D == 2 ? -1 : 0));
  const T p1 = gf(i + (D == 0 ? 1 : 0), j + (D == 1 ? 1 : 0), k + (D == 2 ? 1 : 0));
  return
    (dxi[D]*(-m1 + p1))/2.;
}

template <int D, typename T>
CCTK_DEVICE CCTK_HOST CCTK_ATTRIBUTE_ALWAYS_INLINE inline T
fd_1_o4(amrex::Array4<T const> const &gf, int i, int j, int k, int comp, amrex::GpuArray<T, 3> const &dxi) {
  const T m2 = gf(i + (D == 0 ? -2 : 0), j + (D == 1 ? -2 : 0), k + (D == 2 ? -2 : 0));
  const T m1 = gf(i + (D == 0 ? -1 : 0), j + (D == 1 ? -1 : 0), k + (D == 2 ? -1 : 0));
  const T p1 = gf(i + (D == 0 ? 1 : 0), j + (D == 1 ? 1 : 0), k + (D == 2 ? 1 : 0));
  const T p2 = gf(i + (D == 0 ? 2 : 0), j + (D == 1 ? 2 : 0), k + (D == 2 ? 2 : 0));
  return
    (dxi[D]*(-8*m1 + m2 + 8*p1 - p2))/12.;
}

template <int D, typename T>
CCTK_DEVICE CCTK_HOST CCTK_ATTRIBUTE_ALWAYS_INLINE inline T
fd_1_o6(amrex::Array4<T const> const &gf, int i, int j, int k, int comp, amrex::GpuArray<T, 3> const &dxi) {
  const T m3 = gf(i + (D == 0 ? -3 : 0), j + (D == 1 ? -3 : 0), k + (D == 2 ? -3 : 0));
  const T m2 = gf(i + (D == 0 ? -2 : 0), j + (D == 1 ? -2 : 0), k + (D == 2 ? -2 : 0));
  const T m1 = gf(i + (D == 0 ? -1 : 0), j + (D == 1 ? -1 : 0), k + (D == 2 ? -1 : 0));
  const T p1 = gf(i + (D == 0 ? 1 : 0), j + (D == 1 ? 1 : 0), k + (D == 2 ? 1 : 0));
  const T p2 = gf(i + (D == 0 ? 2 : 0), j + (D == 1 ? 2 : 0), k + (D == 2 ? 2 : 0));
  const T p3 = gf(i + (D == 0 ? 3 : 0), j + (D == 1 ? 3 : 0), k + (D == 2 ? 3 : 0));
  return
    (dxi[D]*(-45*m1 + 9*m2 - m3 + 45*p1 - 9*p2 + p3))/60.;
}

template <int D, typename T>
CCTK_DEVICE CCTK_HOST CCTK_ATTRIBUTE_ALWAYS_INLINE inline T
fd_1_o8(amrex::Array4<T const> const &gf, int i, int j, int k, int comp, amrex::GpuArray<T, 3> const &dxi) {
  const T m4 = gf(i + (D == 0 ? -4 : 0), j + (D == 1 ? -4 : 0), k + (D == 2 ? -4 : 0));
  const T m3 = gf(i + (D == 0 ? -3 : 0), j + (D == 1 ? -3 : 0), k + (D == 2 ? -3 : 0));
  const T m2 = gf(i + (D == 0 ? -2 : 0), j + (D == 1 ? -2 : 0), k + (D == 2 ? -2 : 0));
  const T m1 = gf(i + (D == 0 ? -1 : 0), j + (D == 1 ? -1 : 0), k + (D == 2 ? -1 : 0));
  const T p1 = gf(i + (D == 0 ? 1 : 0), j + (D == 1 ? 1 : 0), k + (D == 2 ? 1 : 0));
  const T p2 = gf(i + (D == 0 ? 2 : 0), j + (D == 1 ? 2 : 0), k + (D == 2 ? 2 : 0));
  const T p3 = gf(i + (D == 0 ? 3 : 0), j + (D == 1 ? 3 : 0), k + (D == 2 ? 3 : 0));
  const T p4 = gf(i + (D == 0 ? 4 : 0), j + (D == 1 ? 4 : 0), k + (D == 2 ? 4 : 0));
  return
    (dxi[D]*(-672*m1 + 168*m2 - 32*m3 + 3*m4 + 672*p1 - 168*p2 + 32*p3 - 3*p4))/840.;
}

} // namespace Particles

#endif // #ifndef PARTICLES_DERIVSINLINE_HXX

/* particles_derivsinline.hxx */
