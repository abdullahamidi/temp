/*---------------------------------------------------------------------------*\

    OSCFD Group

\*---------------------------------------------------------------------------*/

#include "soundPressureLevel.H"
#include "fvMesh.H"
#include "OFstream.H"
#include "IFstream.H"
#include "IOmanip.H"
#include "mathematicalConstants.H"
#include "volFields.H"
#include "addToRunTimeSelectionTable.H"
#include "probeData.H"
#include "signalProcess.H"
#include "writeCSV.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam {
namespace functionObjects {

defineTypeNameAndDebug(soundPressureLevel, 0);
addToRunTimeSelectionTable(functionObject, soundPressureLevel, dictionary);


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void soundPressureLevel::writeFileHeader(Ostream& os) {
    writeHeader(os, "Turbulent Kinetic Energy Spectra");

    writeCommented(os, "kappa E(kappa)");

    os << endl;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

soundPressureLevel::soundPressureLevel
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
    :
    fvMeshFunctionObject(name, runTime, dict),
    writeFile(mesh_, name),
    UProbeName_(""),
    begin_time_(0),
    window_count_(1),
    window_overlap_(0),
    delta_t_(0),
data_count_(0) {
    read(dict);
    Info << "soundPressureLevel" << nl;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool soundPressureLevel::read(const dictionary& dict) {
    fvMeshFunctionObject::read(dict);
    writeFile::read(dict);

    if (!dict.found("probeName"))
        FatalError << "Couldn't find entry probeName for " << name() << "function object" << nl;
    dict.readEntry("probeName", UProbeName_);
    dict.readIfPresent("beginTime", begin_time_);
    dict.readIfPresent("windowCount", window_count_);
    dict.readIfPresent("windowOverlap", window_overlap_);

    if (mesh_.time().startTime().value() > begin_time_)
        Warning << "Specified begin time older than simulation restart time. Only the current data will be accounted." << nl;

    std::ostringstream buf;
    buf.setf(ios_base::fmtflags(), ios_base::floatfield);
    buf.precision(6);
    buf << mesh_.time().startTime().value();

    fileName probeSubDir = UProbeName_;

    if (mesh_.name() != polyMesh::defaultRegion) {
        probeSubDir = probeSubDir / mesh_.name();
    }

    probeDir_ = fileName
    (
        mesh_.time().globalPath()
        / functionObject::outputPrefix
        / probeSubDir
        / buf.str()
    );
    probeDir_.clean();

    return true;
}


bool soundPressureLevel::execute() {
    return true;
}

void soundPressureLevel::read_probe() {
    Info << "Reading " << probeDir_ / "U" << nl;
    vdata_ = readVectorProbe(probeDir_ / "U", begin_time_, delta_t_, data_count_);
}

void soundPressureLevel::do_calculate(const Field<vectorField>& vdata, label data_count, scalar delta_t, label window_count, scalar window_overlap, fileName output) {
    if (data_count == 0) return;
    
    Field<scalarField> u_prime(vdata.size());
    Field<scalarField> spectra(vdata.size());
    scalarField frequencies = fft_frequencies(data_count, window_count, window_overlap, delta_t);
    label fft_result_size = frequencies.size();

    forAll(vdata, i) {
        scalar Umean = average(mag(vdata[i]));
        
        u_prime[i] = 0.5 * magSqr(vdata[i] - average(vdata[i]));
        u_prime[i] = (2.0 / 3.0) * u_prime[i];
        u_prime[i] = sqrt(u_prime[i]);

        Info << "Probe " << Foam::setw(2) << i << ": " << "Umean = " << Foam::setw(7) << Umean << "    u' = " << Foam::setw(7) << average(u_prime[i]) << nl;
    }
    
    forAll(u_prime, i) {
        label data_count;
        spectra[i] = spectral_density(u_prime[i], delta_t, window_count, window_overlap, data_count);
    }

    {
        std::string firstHeader = "Frequency";
        std::vector<std::string> headers(spectra.size());
        for (label i = 0; i < spectra.size(); ++i) {
            headers[i] = "Probe " + std::to_string(i);
        }

        Foam::write_csv(output, fft_result_size, firstHeader, frequencies, headers, spectra);
    }

    Info << nl << nl << nl;
}

bool soundPressureLevel::write() {
    read_probe();
    do_calculate(vdata_, data_count_, delta_t_, window_count_, window_overlap_, probeDir_ / name());

    return true;
}

}
}


// ************************************************************************* //
