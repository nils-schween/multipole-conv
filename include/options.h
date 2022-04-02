#ifndef OPTIONS_H
#define OPTIONS_H

namespace multipole_conv {
enum class MPOptions {
  none = 0,
  complex = 1 << 0,        // use complex solid harmonics
  complex_conjugate = 1 << 1,  // use the complex conjugate of the solid harmonics 
  normalisation = 1 << 2,  // normalise (real/complex) solid harmonics
  remove_condon_shortley_phase = 1 << 3,  // remove the condon shortley phase
  // include the addition theorem factor in the definition of the spherical
  // multipole moments: (2l + 1)/4pi (with normalisation) and (2/(1 + \delta_m0)
  // * (l - m)!/(l + m)! (without normalisation)
  include_addition_theorem = 1 << 4,
  split_addition_theorem = 1 << 5,  // Include the square root of the addition
                                    // theorem factor in the definiton spherical
                                    // multipole moments
  cartesian = 1 << 6,               // cartesian mps in terms of spherical mps
  dependent_components = 1 << 7,    // compute the dependent components of the
                                    // cartesian mps as well
  include_l_factorial = 1 << 8,     //  Include l factorial in the definition of
                                    //  the cartesian mps
  default_case = 1 << 9		    // no options set: real solid harmonics
				    // (without normalisation)
};

constexpr MPOptions operator|(MPOptions f1, MPOptions f2) {
  return static_cast<MPOptions>(static_cast<int>(f1) | static_cast<int>(f2));
}

constexpr MPOptions operator&(MPOptions f1, MPOptions f2) {
  return static_cast<MPOptions>(static_cast<int>(f1) & static_cast<int>(f2));
}

template <typename StreamType>
inline StreamType &operator<<(StreamType &s, MPOptions f) {
  s << "Used options: \n";
  if ((f & MPOptions::complex) != MPOptions::none)
    s << " - Complex solid harmonics are used\n";
  if ((f & MPOptions::complex_conjugate) != MPOptions::none)
    s << " - Complex conjugate solid harmonics are used\n";
  if ((f & MPOptions::normalisation) != MPOptions::none)
    s << " - The solid harmonics are normalised\n";
  if ((f & MPOptions::remove_condon_shortley_phase) != MPOptions::none) {
    s << " - The Condon-Shortley phase is removed from the spherical "
         "harmonics\n";
  }
  if ((f & MPOptions::include_addition_theorem) != MPOptions::none) {
    s << " - An additional factor coming from the addition theorem is "
         "included in the definition of the spherical multipole moments\n";
  }
  if ((f & MPOptions::split_addition_theorem) != MPOptions::none) {
    s << " - The square root of the additional factor coming from the "
         "addition theorem is included in the definition of the spherical "
         "multipole moments\n";
  }
  if ((f & MPOptions::cartesian) != MPOptions::none) {
    s << " - The Cartesian multipole moments are expressed as a linear "
         "combination of the spherical multipole moments.\n";
  }
  if ((f & MPOptions::dependent_components) != MPOptions::none) {
    s << " - The dependent components of the Cartesian multipole moments "
         "are computed as well.\n";
  }
  if ((f & MPOptions::include_l_factorial) != MPOptions::none) {
    s << " - The l! from the Taylor expansion is included in the "
         "defintion of the Cartesian multipole moments\n";
  }
  return s;
}
}  // namespace multipole_conv
#endif
