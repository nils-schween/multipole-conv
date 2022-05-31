#ifndef CONVENTIONS_H
#define CONVENTIONS_H
#include "options.h"

namespace multipole_conv {
namespace conventions {
constexpr MPOptions solid_harmonics =
    MPOptions::complex | MPOptions::normalisation;

constexpr MPOptions jackson = MPOptions::complex | MPOptions::normalisation |
                              MPOptions::complex_conjugate |
                              MPOptions::dependent_components;
constexpr MPOptions johnston =
    MPOptions::cartesian | MPOptions::include_addition_theorem |
    MPOptions::remove_condon_shortley_phase | MPOptions::include_l_factorial |
    MPOptions::dependent_components;

constexpr MPOptions real_solid_harmonics = MPOptions::normalisation;
}  // namespace conventions
}  // namespace multipole_conv
#endif
