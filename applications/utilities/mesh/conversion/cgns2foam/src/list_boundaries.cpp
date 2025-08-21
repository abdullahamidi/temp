
#include <boost/log/trivial.hpp>
#include <boost/program_options.hpp>
#include <cgnslib.h>
#include <boost/program_options/options_description.hpp>
#include <filesystem>
#include <iostream>
#include <string>
#include <fstream>


#include "cgns/boundary_list.hpp"

namespace {
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
    if (! opts.output.empty()) {

        auto output_dir = std::filesystem::path(opts.output).parent_path();
        try {
            if (! fs::exists(output_dir)) {
                fs::create_directories(output_dir);
                BOOST_LOG_TRIVIAL(trace)
                  << "Created output directory: " << output_dir;
            }
        } catch (fs::filesystem_error const &e) {

            BOOST_LOG_TRIVIAL(fatal)
              << "Error: Cannot create output directory '" << output_dir
              << "': " << e.what();
            return false;
        }
    }

    return true;
}
}

int do_list_boundaries(int argc, char *argv[]) {
    program_options opts;

    try {
        // Define command line options
        po::options_description desc = po::
          options_description("List boundaries provided in CGNS mesh file");

        desc.add_options()("help,h", "Show this help message")(
          "input,i",
          po::value<std::string>(&opts.input_file)->required(),
          "Path to input CGNS file")(
          "output,o",
          po::value<std::string>(&opts.output),
          "Output file path")("base,b",
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


        if (opts.verbose) {
            BOOST_LOG_TRIVIAL(trace)
              << "Command line parsing completed successfully.";
            BOOST_LOG_TRIVIAL(trace)
              << "Ready for CGNS mesh processing implementation.";
        }

        CGNS::BoundaryList bl(opts.input_file, opts.base_name, opts.zone_name);

        std::string json = boost::json::serialize(to_json(bl));

        if (opts.output.empty()) {
            std::cout << json << std::endl;
        } else {
            std::ofstream out(opts.output);
            out << json;
            out.close();
        }

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
