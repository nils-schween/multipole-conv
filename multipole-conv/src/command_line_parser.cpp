#include <boost/program_options.hpp>
#include <cmath>
#include <cstddef>
#include <iostream>
#include <utility>

#include "options.h"

namespace multipole_conv {
std::pair<std::size_t, MPOptions> cmd_parser(int argc, char *argv[]) {
  using namespace boost::program_options;
  std::pair<std::size_t, MPOptions> options{0, MPOptions::none};
  try {
    options_description desc{"Multipole Converter\nOptions"};
    desc.add_options()("help,h", "Help screen")(
        "degree,d", value<int>()->default_value(2),
        "Degree of spherical or Cartesian multipole moments.")(
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
    // notify(vm);
    if (vm.count("help")) std::cout << desc << '\n';
    if (vm.count("degree")) {
      std::cout
          << "The degree of the (spherical/Cartesian) multipole moment is: "
          << vm["degree"].as<int>() << '\n';
      options.first = static_cast<std::size_t>(vm["degree"].as<int>());
    }
    if (vm.count("complex"))
      if (vm["complex"].as<bool>())
        options.second = options.second | MPOptions::complex;
    if (vm.count("complex-conjugate"))
      if (vm["complex-conjugate"].as<bool>())
        options.second = options.second | MPOptions::complex_conjugate;
    if (vm.count("normalisation"))
      if (vm["normalisation"].as<bool>())
        options.second = options.second | MPOptions::normalisation;
    if (vm.count("remove-csp"))
      if (vm["remove-csp"].as<bool>())
        options.second =
            options.second | MPOptions::remove_condon_shortley_phase;
    if (vm.count("include-addition-thm"))
      if (vm["include-addition-thm"].as<bool>())
        options.second = options.second | MPOptions::include_addition_theorem;
    if (vm.count("split-addition-thm"))
      if (vm["split-addition-thm"].as<bool>())
        options.second = options.second | MPOptions::split_addition_theorem;
    if (vm.count("cartesian"))
      if (vm["cartesian"].as<bool>())
        options.second = options.second | MPOptions::cartesian;
    if (vm.count("dependent-components"))
      if (vm["dependent-components"].as<bool>())
        options.second = options.second | MPOptions::dependent_components;
    if (vm.count("include-l-factorial"))
      if (vm["include-l-factorial"].as<bool>())
        options.second = options.second | MPOptions::include_l_factorial;

  } catch (const error &ex) {
    std::cerr << ex.what() << '\n';
  }
  return options;
}
}  // namespace multipole_conv
