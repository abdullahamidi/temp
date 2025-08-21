#include <boost/log/trivial.hpp>
#include <boost/program_options.hpp>
#include <cgnslib.h>
#include <boost/program_options/options_description.hpp>
#include <filesystem>
#include <iostream>
#include <string>

#include "cgns/source.hpp"
#include "openfoam/writer.hpp"


namespace po = boost::program_options;
namespace fs = std::filesystem;

struct program_options {
    std::string input_file;
    std::string output;
    std::string base_name;
    std::string zone_name;
    bool verbose = false;
};

bool validate_arguments(program_options const &opts) {

    // Check if input file exists
    if (! fs::exists(opts.input_file)) {
        BOOST_LOG_TRIVIAL(fatal)
          << "Error: Input file '" << opts.input_file << "' does not exist.";

        return false;
    }

    // Check if input file has .cgns extension
    if (fs::path(opts.input_file).extension() != ".cgns") {
        BOOST_LOG_TRIVIAL(fatal)
          << "Warning: Input file does not have .cgns extension.";
    }

    // Create output directory if it doesn't exist
    try {
        if (! fs::exists(opts.output)) {
            fs::create_directories(opts.output);
            BOOST_LOG_TRIVIAL(trace)
              << "Created output directory: " << opts.output;
        }
    } catch (fs::filesystem_error const &e) {

        BOOST_LOG_TRIVIAL(fatal)
          << "Error: Cannot create output directory '" << opts.output
          << "': " << e.what();
        return false;
    }

    return true;
}

int do_convert(int argc, char *argv[]) {
    program_options opts;

    try {
        // Define command line options
        po::options_description
          desc = po::options_description("Convert CGNS mesh to OpenFOAM "
                                         "polyMesh");

        desc.add_options()("help,h", "Show this help message")(
          "input,i",
          po::value<std::string>(&opts.input_file)->required(),
          "Path to input CGNS file")(
          "output,o",
          po::value<std::string>(&opts.output)->required(),
          "Output path")("base,b",
                         po::value<std::string>(&opts.base_name),
                         "Base name to read")(
          "zone,z",
          po::value<std::string>(&opts.zone_name),
          "Zone name to read")("verbose,v",
                               po::bool_switch(&opts.verbose),
                               "Enable verbose output");

        // Parse command line arguments
        po::variables_map vm;
        po::store(po::parse_command_line(argc, argv, desc), vm);

        // Handle help option
        if (vm.contains("help")) {
            BOOST_LOG_TRIVIAL(info) << desc << std::endl;
            return 0;
        }

        // Notify to trigger required argument validation
        po::notify(vm);

        // Validate arguments
        if (! validate_arguments(opts)) { return 1; }


        BOOST_LOG_TRIVIAL(trace)
          << "Command line parsing completed successfully.";


        CGNS::OpenfoamSource
          source(opts.input_file, opts.base_name, opts.zone_name);
        source.load();
        openfoam::Writer writer(source, opts.output);
        writer.write();


    } catch (po::error const &e) {
        BOOST_LOG_TRIVIAL(fatal) << "Error: " << e.what();
        BOOST_LOG_TRIVIAL(fatal) << "Use --help for usage information.";
        return 1;
    } catch (std::exception const &e) {
        BOOST_LOG_TRIVIAL(fatal) << "Error: " << e.what();
        return 1;
    }

    return 0;
}
