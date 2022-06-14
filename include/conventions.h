/*
 * Copyright (C) 2022 Schween, Nils W. and Reville, Brian
 *
 * This file is part of multipole-conv.
 *
 * multipole-conv is free software: you can redistribute it and/or modify it
 * under the terms of the GNU Affero General Public License as published by the
 * Free Software Foundation, either version 3 of the License, or (at your
 * option) any later version.
 *
 * multipole-conv is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the GNU Affero General Public License for more
 * details.
 *
 * You should have received a copy of the GNU Affero General Public License
 * along with multipole-conv. If not, see <https://www.gnu.org/licenses/>.
 */

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
