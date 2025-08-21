#include "messageStream.H"
#include "scalarField.H"
#include "fileName.H"
#include "OFstream.H"
#include <iostream>
#include <fstream>
#include <vector>
#include <limits>
#include <cstdlib>
#include "args.H"
#include "writeCSV.H"
#include "turbulentKineticEnergySpectra.H"
#include "probeData.H"
#include "signalProcess.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void processPressure();

int main(int argc, char* argv[]) {
    static constexpr int ARG_CWD_INDEX = 0;
    static constexpr int ARG_TYPE_INDEX = 1;
    static constexpr int ARG_PATH_INDEX = 2;
    static constexpr int ARG_C_MANDOTARY = 3;

    RUNTIME_CHECK(argc >= ARG_C_MANDOTARY, usage);

    working_directory = argv[ARG_CWD_INDEX];
    set_type(argv[ARG_TYPE_INDEX]);
    path = argv[ARG_PATH_INDEX];

    for (int i = ARG_C_MANDOTARY; i < argc; i = i + 2) {
        set_option(argv[i], argv[i + 1]);
    }

    if (type == E_TYPE_ARG::PRESSURE) {
        processPressure();
    } else if (type == E_TYPE_ARG::VELOCITY) {
        Foam::scalar delta_t;
        Foam::label count;
        Foam::Field<Foam::vectorField> vdata = Foam::readVectorProbe(path, begin_time, delta_t, count);
        Foam::functionObjects::turbulentKineticEnergySpectra::do_calculate(vdata, count, delta_t, window_count, window_overlap, path + output_prefix + "-psd.csv");
    }

    return 0;
}

void processPressure() {
    Foam::scalar delta_t;
    Foam::label count;
    Foam::Field<Foam::scalarField> data = Foam::readScalarProbe(path, begin_time, delta_t, count);
    Foam::Field<Foam::scalarField> spl(data.size());
    Foam::scalarField oaspl1(data.size());
    Foam::scalarField oaspl2(data.size());
    Foam::scalarField frequencies = Foam::fft_frequencies(count, window_count, window_overlap, delta_t);
    Foam::scalarField time(count);
    Foam::labelList indices(data.size());

    std::iota(indices.begin(), indices.end(), 0);
    std::iota(time.begin(), time.end(), 0);

    time *= delta_t;

    Foam::Info << "Square of reference pressure is " << pref_sqr << Foam::nl;

    Foam::Info << count << " data going to be processed with " << window_count << " windows and " << window_overlap * 100 << "\% overlapping" << Foam::nl;


    forAll(data, i) {
        Foam::scalar avg = Foam::average(data[i]);

        oaspl1[i] = 10.0 * std::log10(Foam::average(Foam::sqr(data[i] - avg)) / pref_sqr);

        Foam::label data_count = 0;

        auto spectrum = single_sided_spectrum(data[i], window_count, window_overlap, data_count);

        spl[i] = 10.0 * Foam::log10(sqr(spectrum) / (2 * pref_sqr));
        spl[i][0] = 0;

        oaspl2[i] = 10.0 * std::log10(sum(Foam::pow(10, spl[i] / 10.0)));
    }

    Foam::Info << "Exporting FFT data until " << max_frequency << " Hz" << Foam::nl;

    Foam::label nfreq;
    for (nfreq = 0; nfreq < frequencies.size(); nfreq++)
        if (frequencies[nfreq] > max_frequency)
            break;

    {
        std::string firstHeader = "Frequency [Hz]";
        Foam::scalarField& firstColumn = frequencies;
        std::vector<std::string> headers(spl.size());
        for (Foam::label i = 0; i < spl.size(); ++i) headers[i] = "SPL @ Probe " + std::to_string(i) + " [dB]";

        write_csv(path + output_prefix + "-spl.csv", nfreq, firstHeader, firstColumn, headers, spl);
    }

    {
        std::string firstHeader = "Probe Index";
        std::vector<std::string> headers(1);
        headers[0] = "OASPL [dB]";
        Foam::Field<Foam::scalarField> data(1);
        data[0] = oaspl1;

        write_csv(path + output_prefix + "-oaspl1.csv", oaspl1.size(), firstHeader, indices, headers, data);
    }

    {
        std::string firstHeader = "Probe Index";
        std::vector<std::string> headers(1);
        headers[0] = "OASPL [dB]";
        Foam::Field<Foam::scalarField> data(1);
        data[0] = oaspl2;

        write_csv(path + output_prefix + "-oaspl2.csv", oaspl2.size(), firstHeader, indices, headers, data);
    }
}