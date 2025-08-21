/*---------------------------------------------------------------------------*\

                        Open Source CFD Group

\*---------------------------------------------------------------------------*/

#include "signalProcess.H"
#include "Field.H"
#include "SubField.H"
#include <fftw3.h>
#include <fstream>
#include <math.h>
#include <limits>

namespace Foam {

scalarField fftmagsqr(const scalar* _in, label n) {
    const label nBy2 = n / 2;
    std::vector<double> in(n);
    std::vector<double> out(n);

    for (size_t i = 0; i < in.size(); ++i) {
        in[i] = _in[i];
    }

    // Using real to half-complex fftw 'kind'
    fftw_plan plan =
        fftw_plan_r2r_1d
        (
            n,
            in.data(),
            out.data(),
            FFTW_R2HC,
            FFTW_ESTIMATE
        );

    ::fftw_execute(plan);

    scalarField result(nBy2 + 1);

    // 0 th value = DC component
    // nBy2 th value = real only if n is even
    result[0] = out[0] * out[0];
    result[nBy2] = out[nBy2] * out[nBy2];
    for (label i = 1; i < nBy2; ++i) {
        const auto re = out[i];
        const auto im = out[n - i];
        result[i] = (re * re + im * im);
    }

    fftw_destroy_plan(plan);

    return result;
}

scalarField fftmag(const scalar* in, label n) {
    return Foam::sqrt(fftmagsqr(in, n)).ref();
}

scalarField fftmagsqr(const scalarField& in) {
    return fftmagsqr(in.cdata(), in.size());
}

scalarField fftmag(const scalarField& in) {
    return Foam::sqrt(fftmagsqr(in)).ref();
}

scalarField fftmagsqr(const scalarField& in, label offset, label count) {
    if (offset + count > in.size()) {
        FatalErrorInFunction
            << "Expected " << offset + count
            << " but got" << in.size()
            << abort(FatalError);
    }
    
    return fftmagsqr(in.cdata() + offset, count);
}

scalarField fftmag(const scalarField& in, label offset, label count) {
    return Foam::sqrt(fftmagsqr(in, offset, count)).ref();
}

scalarField fftmagsqr(const scalarField& in, label nWindow, scalar overlap, label& data_count_per_window) {
    if (in.size() == 0) return {};
    
    label data_count = in.size();
    data_count_per_window = std::floor(scalar(data_count) / (scalar(nWindow) * (1.0 - overlap) + overlap));
    label overlap_count = std::ceil(data_count_per_window * overlap);
    label fft_result_size = data_count_per_window / 2 + 1;
    //scalar actual_overlap = scalar(overlap_count) / scalar(data_count_per_window);

    // Info << overlap * 100 << "\% overlap specified but must apply " << actual_overlap * 100 << "\% to take " << data_count_per_window << " data per window" << nl;

    scalarField result(fft_result_size);

    for (label i = 0; i < nWindow; ++i) {
        label data_offset = i * (data_count_per_window - overlap_count);
        label data_count = data_count_per_window;
        
        // Info << "Processing " << i + 1 << "th window in range [" << data_offset << ", " << data_offset + data_count_per_window << ")" << nl;

        if (data_offset + data_count > in.size()) {
            FatalErrorInFunction
                << "Expected at least " << data_offset + data_count
                << " but got " << in.size() << " data"
                << abort(FatalError);
        }
        
        result += fftmagsqr(in.cdata() + data_offset, data_count_per_window);
    }

    result /= scalar(nWindow);

    return result;
}

scalarField fftmag(const scalarField& in, label nWindow, scalar overlap, label& data_count_per_window) {
    return Foam::sqrt(fftmagsqr(in, nWindow, overlap, data_count_per_window)).ref();
}

scalarField fft_frequencies(label data_count, label nWindow, scalar overlap, scalar delta_T)
{
    label data_count_per_window = std::floor(scalar(data_count) / (scalar(nWindow) * (1.0 - overlap) + overlap));
    label fft_result_size = data_count_per_window / 2 + 1;

    scalarField f(fft_result_size);

    const scalar delta_f = 1 / (data_count_per_window * delta_T);

    for (label i = 0; i < f.size(); ++i) f[i] = i * delta_f;

    return f;
}

scalarField single_sided_spectrum(const scalarField& in, label nWindow, scalar overlap, label& data_count_per_window) {
    auto result = fftmag(in, nWindow, overlap, data_count_per_window);

    result *= 2.0;
    result[0] *= 0.5;
    result[result.size() - 1] *= 0.5;
    
    return (result / data_count_per_window).ref();
}

scalarField single_sided_psd(const scalarField& in, scalar dt, label nWindow, scalar overlap, label& data_count_per_window) {
    return (sqr(single_sided_spectrum(in, nWindow, overlap, data_count_per_window)) * dt).ref();
}

}


// ************************************************************************* //
