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

#include <boost/program_options.hpp>
#include <cmath>
#include <cstddef>
#include <iostream>
#include <string>
#include <utility>

#include "conventions.h"
#include "options.h"
#include "version.h"

namespace multipole_conv {
MPOptions set_convention(const std::string &convention) {
  MPOptions options = MPOptions::none;
  if (convention == "solid_harmonics")
    options = conventions::solid_harmonics;
  else if (convention == "jackson")
    options = conventions::jackson;
  else if (convention == "johnston")
    options = conventions::johnston;
  else if (convention == "real_solid_harmonics")
    options = conventions::real_solid_harmonics;
  else
    std::cout << "No matching convention was found.\n";
  return options;
}

std::pair<std::size_t, MPOptions> cmd_parser(int argc, char *argv[]) {
  using namespace boost::program_options;
  std::pair<std::size_t, MPOptions> cmd_options{0, MPOptions::none};
  try {
    options_description desc{"Options"};
    desc.add_options()("help,h", "Help screen")("version,v",
                                                "Displays the version number.")(
        "degree,d", value<int>(),
        "Degree of spherical or Cartesian multipole moments.")(
        "convention,c", value<std::string>(),
        "Conventions are predefined set of options. Possible values are "
        "solid_harmonics, jackson, johnston, real_solid_harmonics.")(
        "complex", value<bool>(), "Use complex solid harmonics.")(
        "complex-conjugate", value<bool>(),
        "Use the complex conjugate of the solid harmonics.")(
        "normalisation", value<bool>(),
        "Normalise the (real or complex) solid harmonics.")(
        "remove-csp", value<bool>(), "Remove the Condon-Shortley phase.")(
        "include-addition-thm", value<bool>(),
        "Include the addition theorem factor in the definition of the "
        "spherical multipole moments.")(
        "split-addition-thm", value<bool>(),
        "Include the square root of the addition theorem factor in the "
        "definition of the spherical multipole moments. (Requires that the "
        "option \"inlcude-addition-thm\" is set.)")(
        "cartesian", value<bool>(), "Compute the Cartesian multipole moments.")(
        "dependent-components", value<bool>(),
        "Compute the dependent components of the Cartesian multipole moments.")(
        "include-l-factorial", value<bool>(),
        "Include l! from the Taylor expansion in the definition of the "
        "Cartesian multipole moment");

    variables_map vm;
    store(parse_command_line(argc, argv, desc), vm);

    if (vm.count("help")) {
      std::cout << desc << '\n';
      return cmd_options;
    }
    if (vm.count("version")) {
      std::cout << "multipole-conv version " << version_major << "."
                << version_minor << "\n\n"
                << "Copyright (C) 2022 Schween, Nils W."
                << "\n"
                << "License AGPLv3+: GNU AGPL version 3 or later "
                   "<https://gnu.org/licenses/gpl.html>"
                << "\n"
                << "This is free software: you are free to change and "
                   "redistribute it."
                << "\n"
                << "This program comes with ABSOLUTELY NO WARRANTY"
                << "\n\n";
      return cmd_options;
    }
    if (vm.count("degree")) {
      std::cout
          << "The degree of the (spherical/Cartesian) multipole moment is: "
          << vm["degree"].as<int>() << '\n';
      cmd_options.first = static_cast<std::size_t>(vm["degree"].as<int>());
      cmd_options.second = cmd_options.second | MPOptions::default_case;
    }
    if (vm.count("convention"))
      cmd_options.second = set_convention(vm["convention"].as<std::string>());
    if (vm.count("complex"))
      if (vm["complex"].as<bool>())
        cmd_options.second = cmd_options.second | MPOptions::complex;
    if (vm.count("complex-conjugate"))
      if (vm["complex-conjugate"].as<bool>())
        cmd_options.second = cmd_options.second | MPOptions::complex_conjugate;
    if (vm.count("normalisation"))
      if (vm["normalisation"].as<bool>())
        cmd_options.second = cmd_options.second | MPOptions::normalisation;
    if (vm.count("remove-csp"))
      if (vm["remove-csp"].as<bool>())
        cmd_options.second =
            cmd_options.second | MPOptions::remove_condon_shortley_phase;
    if (vm.count("include-addition-thm"))
      if (vm["include-addition-thm"].as<bool>())
        cmd_options.second =
            cmd_options.second | MPOptions::include_addition_theorem;
    if (vm.count("split-addition-thm"))
      if (vm["split-addition-thm"].as<bool>())
        cmd_options.second =
            cmd_options.second | MPOptions::split_addition_theorem;
    if (vm.count("cartesian"))
      if (vm["cartesian"].as<bool>())
        cmd_options.second = cmd_options.second | MPOptions::cartesian;
    if (vm.count("dependent-components"))
      if (vm["dependent-components"].as<bool>())
        cmd_options.second =
            cmd_options.second | MPOptions::dependent_components;
    if (vm.count("include-l-factorial"))
      if (vm["include-l-factorial"].as<bool>())
        cmd_options.second =
            cmd_options.second | MPOptions::include_l_factorial;

  } catch (const error &ex) {
    std::cerr << ex.what() << '\n';
  }
  return cmd_options;
}
}  // namespace multipole_conv
