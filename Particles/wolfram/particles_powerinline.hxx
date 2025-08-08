#ifndef PARTICLES_POWERINLINE_HXX
#define PARTICLES_POWERINLINE_HXX

#include <loop_device.hxx>

#include <array>
#include <cmath>
#include <type_traits>

namespace Particles {

template <typename T, typename U>
CCTK_DEVICE CCTK_HOST constexpr T Power(T x, U y) {
  // If y is an integer, use exponentiation by squaring
  if constexpr (std::is_integral<U>::value) {
    if (y < 0) {
      return static_cast<T>(1) / Power(x, -y); // Handle negative exponents
    }
    T result = 1;
    while (y) {
      if (y % 2) {
        result *= x;
      }
      x *= x;
      y /= 2;
    }
    return result;
  } else {
    // Use std::pow for non-integer exponents
    return std::pow(x, y);
  }
}

} // namespace Particles

#endif // #ifndef PARTICLES_POWERINLINE_HXX
