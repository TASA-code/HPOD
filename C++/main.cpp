#include "Orbit.h"
#include "Coordinate.h"
#include "Time.h"
#include <iostream>

// add -lboost_program_options to compile
#include <boost/program_options.hpp>
namespace po = boost::program_options;

int main(int argc, char *argv[]) {
    po::options_description opts("Allowed options");
    opts.add_options()("h", "produce help message")

    ("SMA", po::value<double>()->required(), "Semi-Major Axis (double)")

    ("e", po::value<double>()->required(), "Eccentricity (double)")

    ("i", po::value<double>()->required(), "inclination (deg.)")
                  
    ("M", po::value<double>()->required(),"Mean Anomaly (deg.)")

    ("w", po::value<double>()->required(),"Argument of Periapsis (deg.)")

    ("RAAN", po::value<double>()->required(),"Right Ascension Ascending Node (deg.)");

  po::variables_map vm;

  // Handles errors with the passed arguments and provides helpful error
  // messages.
  try {
    po::store(po::parse_command_line(argc, argv, opts), vm);
    // If the user asks for help, or no arguments are supplied, provide
    // a help message and list the possible arguments.
    if (vm.count("h") || argc == 1) {
      std::cout << "Please provide a value for all the required arguments "
                << std::endl
                << opts << std::endl;
      return 0; // Exiting with 0 as error was properly handled.
    }

    po::notify(vm);
  } catch (std::exception &e) {
    std::cout << "Error: " << e.what() << std::endl;
    std::cout << opts << std::endl;
    return 0;
  }

  // Extracts the parameter values given to variables using the appropriate
  // datatype.
  const double SMA = vm["SMA"].as<double>();
  const double e = vm["e"].as<double>();
  const double i = vm["i"].as<double>();
  const double M = vm["M"].as<double>();
  const double w = vm["w"].as<double>();
  const double RAAN = vm["RAAN"].as<double>();

  Orbit orbit_prop;

  orbit_prop.SetParameter(SMA, e, i, M, w, RAAN);

  orbit_prop.integrate();

  // orbit_prop.~Orbit();
}
