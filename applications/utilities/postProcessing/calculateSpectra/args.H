/*---------------------------------------------------------------------------*\

                        Open Source CFD Group

\*---------------------------------------------------------------------------*/

#ifndef PROCESS_SIGNAL_ARGS_H
#define PROCESS_SIGNAL_ARGS_H

#define PRINT_MSG_TERMINATE(msg) { Foam::Info << msg << Foam::nl; std::terminate(); }
#define RUNTIME_CHECK(x, msg) if(!(x)) { PRINT_MSG_TERMINATE(msg) }


static constexpr const char* usage =
"====================================================================\n"
"Open Source CFD Group\n"
"====================================================================\n"
"Usage :\n"
"calculateSpectra (input file type) (input file path) [option value] \n"
"\t --begin-time \t\t: Processing the file from the specified time value \n"
"\t --window-count \t: How many windows going to be applied to the source data before FFT \n"
"\t --window-overlap \t: FFT calculations going to be made by windowing the source data \n"
"\t --max-frequency \t: Print FFT results until the specified frequency \n"
"\t --output-prefix \t: Prefix for output files \n";

enum E_TYPE_ARG : long long int {
    PRESSURE,
    VELOCITY
} type;

std::string working_directory;
std::string path;
std::string output_prefix;
Foam::scalar begin_time = 0;
Foam::scalar max_frequency = std::numeric_limits<Foam::scalar>::max();
Foam::label window_count = 1;
Foam::scalar window_overlap = 0;
Foam::scalar pref_sqr = 4e-10;

void set_type(std::string selected_type) {
    if (selected_type == "pressure")
        type = E_TYPE_ARG::PRESSURE;
    else if (selected_type == "velocity")
        type = E_TYPE_ARG::VELOCITY;
    else PRINT_MSG_TERMINATE(usage);
}

void set_option(std::string option, std::string value) {
    if (option == "--begin-time") {
        begin_time = std::stod(value);
        RUNTIME_CHECK(begin_time >= 0, "Invalid begin-time arg");
    } else if (option == "--max-frequency") {
        max_frequency = std::stod(value);
        RUNTIME_CHECK(max_frequency >= 0, "Invalid max-frequency arg");
    } else if (option == "--window-count") {
        window_count = std::stoull(value);
        RUNTIME_CHECK(window_count > 0, "Invalid window count arg");
    } else if (option == "--window-overlap") {
        window_overlap = std::stod(value);
        RUNTIME_CHECK(0 <= window_overlap && window_overlap <= 1, "Invalid windowing overlap arg");
    } else if (option == "--pref-square") {
        pref_sqr = std::stod(value);
    } else if (option == "--output-prefix") {
        output_prefix = value;
    } else {
        Foam::Info << "Unknown option " << option << " passed" << Foam::nl;
    }
}


#endif
