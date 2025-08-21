
#include <cgnslib.h>
#include <iostream>
#include <string>
#include <boost/log/trivial.hpp>
#include <boost/program_options.hpp>
#include <boost/program_options/options_description.hpp>
#include "macro.hpp"


namespace po = boost::program_options;
namespace fs = std::filesystem;


static constexpr char const *subcommands[] = { "convert", "list-boundaries" };

enum SUBCOMMAND {
    convert = 0,
    list_boundaries = 1
};

int do_convert(int argc, char *argv[]);

int do_list_boundaries(int argc, char *argv[]);

void print_commands(int _, char *argv[]) {
    std::string joined;
    for (auto const &word: subcommands) {
        if (! joined.empty()) { joined += ", "; }
        joined += word;
    }

    std::cout << "Usage: " << argv[0] << " <subcommand> [options]" << std::endl
              << "\tsubcommand options: " << joined << std::endl;
}

int main(int argc, char *argv[]) {
    if (argc < 2) {
        print_commands(argc, argv);
        return 1;
    }

    std::string subcommand = argv[1];

    // Print CGNS library version info
    BOOST_LOG_TRIVIAL(trace)
      << "CGNS library version: " + std::string(QUOTE(CGNS_DOTVERS));

    if (subcommand == subcommands[SUBCOMMAND::convert]) {
        return do_convert(argc, argv);
    }
    if (subcommand == subcommands[SUBCOMMAND::list_boundaries]) {
        return do_list_boundaries(argc, argv);
    }

    print_commands(argc, argv);
    return 1;
}
