/*---------------------------------------------------------------------------*\

    V7320 - Özgün Mühendislik Yazılımları

-------------------------------------------------------------------------------
Class
    volumeNoise

Description

Authors
    Barış Cumhur
        Turkish Aerospace Industries, Turkey

\*---------------------------------------------------------------------------*/

#ifndef FOAM_VOLUME_NOISE_H
#define FOAM_VOLUME_NOISE_H

#include "ISstream.H"
#include "List.H"
#include "Time.H"
#include "TimeState.H"
#include "UPstream.H"
#include "autoPtr.H"
#include "fileName.H"
#include "noiseModel.H"
#include "scalarField.H"
#include "vectorField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam {

namespace noiseModels {

/*---------------------------------------------------------------------------*\
                        Class volumeNoise Declaration
\*---------------------------------------------------------------------------*/

class volumeNoise : public noiseModel {

protected:
  //- Input file names
  List<word> filteredFieldFunctionObjectNames_;

  //- Name of exported field
  word pName_;

  //- Start time index
  label startTimeIndex_;

  //- Frequency data output interval, default = 1
  //  nSamples/2 data points are returned from the FFT, which can
  //  result in a very large number of output files (1 per frequency)
  label fftWriteInterval_;

  //- Apply volume average; default = no (ensemble average) for backwards
  //- compatibility
  bool volumeAverage_;

  bool parallel_;

  // Protected Member Functions

  bool readData(label functionObjectNameIndex);

  void doCalculate(label functionObjectNameIndex);

  scalar volumeAverage(scalarField const &field) const;

public:
  //- Runtime type information
  TypeName("volumeNoise");

  //- Constructor
  volumeNoise(dictionary const &dict, objectRegistry const &obr,
              word const &name = typeName);

  //- Destructor
  virtual ~volumeNoise() = default;

  static fileName getOutputDir(Time const &time,
                               word const &functionObjectName);

  static fileName getFftOutputDir(Time const &time,
                                  word const &functionObjectName);

  static fileName getFftOutputFilename(word freq);

  static word getFreqFromFftOutputFilename(fileName const &file);

  static fileName get13OctaveOutputDir(Time const &time,
                                       word const &functionObjectName);

  static fileName get13OctaveOutpuFilename(word freq);

  static word getFreqFrom13OctaveOutpuFilename(fileName const &file);

  // Public Member Functions

  //- Read from dictionary
  virtual bool read(dictionary const &dict);

  //- Calculate
  virtual void calculate();

private:
  Time const &currentTime_;

  labelList localIndices;

  labelList globalIndices;

  vectorField centroids;

  scalarField volumes;

  List<scalarField> ptList;

  scalarList times_;

  scalar deltaT_;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace noiseModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
