/*---------------------------------------------------------------------------*\

    V7320 - Özgün Mühendislik Yazılımları

\*---------------------------------------------------------------------------*/

#if OPENFOAM > 2206
#include "ISstream.H"
#include "Pstream.H"
#include "PstreamReduceOps.H"
#include "Time.H"
#include "UPstream.H"
#include "addToRunTimeSelectionTable.H"
#include "csv.H"
#include "instantList.H"
#include "int32.H"
#include "noiseModel.H"
#include "objectRegistry.H"
#include "scalar.H"
#include "scalarField.H"
#include "stdFoam.H"
#include "vectorField.H"
#include "volumeNoise.H"
#include <algorithm>
#include <utility>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam {
namespace noiseModels {

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(volumeNoise, 0);
addToRunTimeSelectionTable(noiseModel, volumeNoise, dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

volumeNoise::volumeNoise(dictionary const &dict, objectRegistry const &obr,
                         word const &name)
    : noiseModel(dict, obr, name, false), pName_("p"), startTimeIndex_(0),
      fftWriteInterval_(1), volumeAverage_(false), currentTime_(obr.time()),
      times_(), deltaT_(0), parallel_(Pstream::is_parallel()) {
  read(dict);
}

fileName volumeNoise::getOutputDir(Time const &time,
                                   word const &functionObjectName) {
  return time.timeName() + "/" + functionObjectName;
}

fileName volumeNoise::getFftOutputDir(Time const &time,
                                      word const &functionObjectName) {
  return getOutputDir(time, functionObjectName) / "fft";
}

fileName volumeNoise::getFftOutputFilename(word freq) {
  return freq + "hz.csv";
}

word volumeNoise::getFreqFromFftOutputFilename(fileName const &file) {
  return file.name().substr(0, file.name().find("hz.csv"));
}

fileName volumeNoise::get13OctaveOutputDir(Time const &time,
                                           word const &functionObjectName) {
  return getOutputDir(time, functionObjectName) / "13octave";
}

fileName volumeNoise::get13OctaveOutpuFilename(word freq) {
  return freq + "hz.csv";
}

word volumeNoise::getFreqFrom13OctaveOutpuFilename(fileName const &file) {
  return file.name().substr(0, file.name().find("hz.csv"));
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool volumeNoise::read(dictionary const &dict) {
  if (noiseModel::read(dict)) {

    dict.readEntry("filteredFieldFunctionObjectNames",
                   filteredFieldFunctionObjectNames_);

    dict.readIfPresent("p", pName_);
    dict.readIfPresent("fftWriteInterval", fftWriteInterval_);

    Info << this->type() << nl << "    Pressure field name: " << pName_ << nl
         << "    FFT write interval: " << fftWriteInterval_ << nl;

    dict.readIfPresent("volumeAverage", volumeAverage_);

    if (volumeAverage_) {
      Info << "    Averaging: volume weighted" << endl;
    } else {
      Info << "    Averaging: ensemble" << endl;
    }

    Info << endl;

    return true;
  }

  return false;
}

bool volumeNoise::readData(label functionObjectIndex) {
  Info << "Reading data for "
       << filteredFieldFunctionObjectNames_[functionObjectIndex] << endl;

  instantList allTimes = currentTime_.times();

#if OPENFOAM > 2206
  startTimeIndex_ = instant::findStart(allTimes, startTime_);
#else
  startTimeIndex_ = findStartTimeIndex(allTimes, startTime_);
#endif

  if (allTimes[startTimeIndex_].value() < startTime_) {
    WarningInFunction << "Could not find start time: " << startTime_ << endl;
    return false;
  }

  label nAvailableTimes = allTimes.size() - startTimeIndex_;

  label nRequiredTimes = windowModelPtr_->validate(nAvailableTimes);

  times_.setSize(nRequiredTimes);

  forAll(times_, timeI) {
    times_[timeI] = allTimes[timeI + startTimeIndex_].value();
  }

  deltaT_ = checkUniformTimeStep(times_);

  localIndices.clear();
  globalIndices.clear();
  centroids.clear();
  volumes.clear();
  ptList.clear();

  for (Foam::label t = startTimeIndex_; t < allTimes.size(); ++t) {
    IOobject read_from(
        filteredFieldFunctionObjectNames_[functionObjectIndex] + ".csv",
        allTimes[t].name(), fileObr_, IOobject::NO_READ, IOobject::NO_WRITE);

    labelList _localIndices;
    labelList _globalIndices;
    vectorField _centroids;
    scalarField _volumes;
    scalarField _field;

    if (!isFile(read_from.objectPath())) {
      WarningInFunction << "Cannot find file: " << read_from.objectPath()
                        << endl;
      continue;
    }

    readCsv(
        read_from, std::make_pair("index", &_localIndices),
        std::make_pair("global_index", parallel_ ? &_globalIndices : nullptr),
        std::make_pair("centroid", &_centroids),
        std::make_pair("volume", &_volumes),
        std::make_pair(pName_.c_str(), &_field));

    if (parallel_) {
      _globalIndices = _localIndices.clone();
    }

    if (globalIndices.empty()) {
      localIndices.resize(_localIndices.size());
      std::copy(_localIndices.begin(), _localIndices.end(),
                localIndices.begin());
      globalIndices.resize(_globalIndices.size());
      std::copy(_globalIndices.begin(), _globalIndices.end(),
                globalIndices.begin());
      centroids.resize(_centroids.size());
      std::copy(_centroids.begin(), _centroids.end(), centroids.begin());
      volumes.resize(_volumes.size());
      std::copy(_volumes.begin(), _volumes.end(), volumes.begin());
      ptList.resize(_field.size());
      forAll(ptList, i) {
        ptList.resize(_field.size());
        ptList[i].resize(nAvailableTimes);
        ptList[i][t - startTimeIndex_] = _field[i];
      }
    } else {
      // Remove array indices from indices, centroids, volumes, and pt_list if
      // indices is not in _global_indices
      auto new_size = globalIndices.size();

      for (label i = 0; i < new_size; i++) {
        auto pos = std::find(_globalIndices.begin(), _globalIndices.end(),
                             globalIndices[i]);
        auto index = pos - _globalIndices.begin();
        if (pos == _globalIndices.end()) {
          std::swap(localIndices[i], localIndices[--new_size]);
          std::swap(globalIndices[i], globalIndices[new_size]);
          std::swap(centroids[i], centroids[new_size]);
          std::swap(volumes[i], volumes[new_size]);
          std::swap(ptList[i], ptList[new_size]);
        }
        ptList[i][t - startTimeIndex_] = _field[index];
      }

      localIndices.resize(new_size);
      globalIndices.resize(new_size);
      centroids.resize(new_size);
      volumes.resize(new_size);
      ptList.resize(new_size);
    }
  }

  Info << "Read data for " << centroids.size() << " points." << endl;

  return true;
}

void volumeNoise::doCalculate(label functionObjectNameIndex) {
  scalarField const freq1(uniformFrequencies(deltaT_, true));

  scalar const maxFreq1 = max(freq1);

  // Reset desired frequency range if outside actual frequency range
  fLower_ = min(fLower_, maxFreq1);
  fUpper_ = min(fUpper_, maxFreq1);

  // Storage for FFT data
  label const nPoints = ptList.size();
  label const nFFT = ceil(freq1.size() / scalar(fftWriteInterval_));

  windowModel const &win = windowModelPtr_();

  List<scalarField> reducedPrmsf(nFFT);
  List<scalarField> reducedPSDf(nFFT);

  forAll(reducedPrmsf, freqI) {
    reducedPrmsf[freqI].setSize(nPoints);
    reducedPSDf[freqI].setSize(nPoints);
  }

  // Storage for 1/3 octave data
  labelList octave13BandIDs;
  scalarField octave13FreqCentre;
  setOctaveBands(freq1, fLower_, fUpper_, 3, octave13BandIDs,
                 octave13FreqCentre);

  label bandSize = 0;
  if (octave13BandIDs.empty()) {
    WarningInFunction << "Octave band calculation failed (zero sized). "
                      << "Please check your input data" << endl;
  } else {
    bandSize = octave13BandIDs.size() - 1;
  }

  Info << "Creating noise FFTs" << endl;

  List<scalarField> reducedPrms13f(bandSize);
  forAll(reducedPrms13f, freqI) { reducedPrms13f[freqI].setSize(nPoints); }

  {
    forAll(ptList, point_index) {
      scalarField const &p = ptList[point_index];

      // Generate the FFT-based data
      scalarField const prmsf(RMSmeanPf(p));
      scalarField const psdf(PSDf(p, deltaT_));

      // Store the frequency results in slot for volume
      forAll(reducedPrmsf, i) {
        label freqI = i * fftWriteInterval_;
        reducedPrmsf[i][point_index] = prmsf[freqI];
        reducedPSDf[i][point_index] = psdf[freqI];
      }

      if (!octave13BandIDs.empty()) {
        // Integrated PSD = P(rms)^2 [Pa^2]
        scalarField const Prms13f(octaves(psdf, freq1, octave13BandIDs));

        // Store the 1/3 octave results in slot for  volume
        forAll(reducedPrms13f, freqI) {
          reducedPrms13f[freqI][point_index] = Prms13f[freqI];
        }
      }
    }
  }

  scalar const deltaf = 1.0 / (deltaT_ * win.nSamples());

  Info << "Writing FFT volume data at every " << fftWriteInterval_
       << " frequency points" << endl;
  {
    // Determine frequency range of interest
    // Note: frequencies have fixed interval, and are in the range
    //       0 to fftWriteInterval_*(n-1)*deltaf
    label f0 = ceil(fLower_ / deltaf / scalar(fftWriteInterval_));
    label f1 = floor(fUpper_ / deltaf / scalar(fftWriteInterval_));
    label nFreq = f1 - f0;

    scalarField PrmsfAve(nFreq, Zero);
    scalarField PSDfAve(nFreq, Zero);
    scalarField fOut(nFreq, Zero);

    if (nFreq == 0) {
      WarningInFunction
          << "No volume data available using a fftWriteInterval of "
          << fftWriteInterval_ << endl;
    } else {
      forAll(fOut, i) {
        label freqI = (i + f0) * fftWriteInterval_;
        fOut[i] = freq1[freqI];

        IOobject write_to(
            getFftOutputFilename(std::to_string(fOut[i])),
            getFftOutputDir(
                currentTime_,
                filteredFieldFunctionObjectNames_[functionObjectNameIndex]),
            fileObr_, IOobject::NO_READ, IOobject::NO_WRITE);

        auto psd = PSD(reducedPSDf[i + f0]);
        auto spl = SPL(reducedPSDf[i + f0] * deltaf, freq1[freqI]);

        writeCsv(
            write_to, localIndices.size(),
            std::make_pair("index", &localIndices),
            std::make_pair("global_index",
                           parallel_ ? &globalIndices : nullptr),
            std::make_pair("centroid", &centroids),
            std::make_pair("volume", &volumes),
            std::make_pair("Prmsf",
                           writePrmsf_ ? &reducedPrmsf[i + f0] : nullptr),
            std::make_pair("PSDf", writePSDf_ ? &reducedPSDf[i + f0] : nullptr),
            std::make_pair("PSD", writePSD_ ? psd.ptr() : nullptr),
            std::make_pair("SPL", writeSPL_ ? spl.ptr() : nullptr));

        PrmsfAve[i] = volumeAverage(reducedPrmsf[i + f0]);
        PSDfAve[i] = volumeAverage(reducedPSDf[i + f0]);
      }
    }

    if (UPstream::master()) {
      Info << "Writing volume averages for FFT" << endl;

      IOobject write_to(
          "fft_avg.csv",
          getOutputDir(
              currentTime_,
              filteredFieldFunctionObjectNames_[functionObjectNameIndex]),
          fileObr_, IOobject::NO_READ, IOobject::NO_WRITE);

      writeCsv(
          write_to, fOut.size(), std::make_pair("f [Hz]", &fOut),
          std::make_pair("P(f) [Pa]", &PrmsfAve),
          std::make_pair("PSD(f) [PaPa_Hz]", &PSDfAve),
          std::make_pair("PSD(f) [dB_Hz]", PSD(PSDfAve).ptr()),
          std::make_pair("SPL(f) [dB]", SPL(PSDfAve * deltaf, fOut).ptr()));
    }
  }

  if (!reducedPrms13f.empty()) {
    Info << "Writing one-third octave volume data" << endl;

    scalarField Prms13fAve(reducedPrms13f.size(), Zero);

    forAll(reducedPrms13f, i) {

      IOobject write_to(
          get13OctaveOutpuFilename(std::to_string(octave13FreqCentre[i])),
          get13OctaveOutputDir(
              currentTime_,
              filteredFieldFunctionObjectNames_[functionObjectNameIndex]),
          fileObr_, IOobject::NO_READ, IOobject::NO_WRITE);

      writeCsv(
          write_to, localIndices.size(), std::make_pair("index", &localIndices),
          std::make_pair("global_index", parallel_ ? &globalIndices : nullptr),
          std::make_pair("centroid", &centroids),
          std::make_pair("volume", &volumes),
          std::make_pair(
              "SPL13", writeOctaves_
                           ? SPL(reducedPrms13f[i], octave13FreqCentre[i]).ptr()
                           : nullptr));

      Prms13fAve[i] = volumeAverage(reducedPrms13f[i]);
    }
    if (UPstream::master()) {
      Info << "Writing one-third octave volume averages" << endl;

      IOobject write_to(
          "13octave_avg.csv",
          getOutputDir(
              currentTime_,
              filteredFieldFunctionObjectNames_[functionObjectNameIndex]),
          fileObr_, IOobject::NO_READ, IOobject::NO_WRITE);

      writeCsv(write_to, octave13FreqCentre.size(),
               std::make_pair("f [Hz]", &octave13FreqCentre),
               std::make_pair("P(fm) [Pa]", &Prms13fAve),
               std::make_pair("SPL(fm) [dB]",
                              SPL(Prms13fAve, octave13FreqCentre).ptr()));
    }
  }
}

scalar volumeNoise::volumeAverage(scalarField const &field) const {
  scalar localWeightedSum;

  scalar localWeightSum;

  if (volumeAverage_) {
    localWeightedSum = sum(field * volumes);
    localWeightSum = sum(volumes);
  } else {
    localWeightedSum = sum(field);
    localWeightSum = field.size();
  }

  label reduceRequest[2];

  reduce(localWeightedSum, sumOp<scalar>(), Pstream::msgType(),
         UPstream::worldComm, reduceRequest[0]);

  reduce(localWeightSum, sumOp<scalar>(), Pstream::msgType(),
         UPstream::worldComm, reduceRequest[1]);

  Pstream::waitRequest(reduceRequest[0]);
  Pstream::waitRequest(reduceRequest[1]);

  return localWeightedSum / localWeightSum;
}

void volumeNoise::calculate() {
  forAll(filteredFieldFunctionObjectNames_, i) {
    if (readData(i)) {
      doCalculate(i);
    }
  }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace noiseModels
} // End namespace Foam

// ************************************************************************* //

#else

#warning("OpenFOAM version not supported with noiseModel classes")

#endif
