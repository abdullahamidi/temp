

/*---------------------------------------------------------------------------*\

                        Open Source CFD Group

\*---------------------------------------------------------------------------*/

#ifndef processSignal_H
#define processSignal_H

#include "label.H"
#include "scalarField.H"
#include "vectorField.H"

namespace Foam {

scalarField fftmagsqr(const scalar* _in, label n);

scalarField fftmag(const scalar* in, label n);

scalarField fftmagsqr(const scalarField& in);

scalarField fftmag(const scalarField& in);

scalarField fftmagsqr(const scalarField& in, label offset, label count);

scalarField fftmag(const scalarField& in, label offset, label count);

scalarField fftmagsqr(const scalarField& in, label nWindow, scalar overlap, label& data_count_per_window);

scalarField fftmag(const scalarField& in, label nWindow, scalar overlap, label& data_count_per_window);

scalarField fft_frequencies(label count, label nWindow, scalar overlap, scalar delta_T);

scalarField single_sided_spectrum(const scalarField& in, label nWindow, scalar overlap, label& data_count_per_window);

scalarField single_sided_psd(const scalarField& in, scalar dt, label nWindow, scalar overlap, label& data_count_per_window);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
