

/*---------------------------------------------------------------------------*\

                        Open Source CFD Group

\*---------------------------------------------------------------------------*/

#ifndef probeData_H
#define probeData_H

#include "fileName.H"
#include "label.H"
#include "scalarField.H"
#include "vectorField.H"

namespace Foam {

Field<scalarField> readScalarProbe(fileName file, scalar begin_time, scalar& delta_t, label& count);

Field<vectorField> readVectorProbe(fileName file, scalar begin_time, scalar& delta_t, label& count);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
