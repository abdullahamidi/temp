/*---------------------------------------------------------------------------*\

    HiSA: High Speed Aerodynamic solver

    Copyright (C) 2014-2018 Oliver Oxtoby - CSIR, South Africa
    Copyright (C) 2014-2018 Johan Heyns - CSIR, South Africa
    Copyright (C) 2004-2010 OpenCFD Ltd.

-------------------------------------------------------------------------------
License
    This file is part of HiSA.

    HiSA is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    HiSA is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with HiSA.  If not, see <http://www.gnu.org/licenses/>.

Class
    Foam::wPureFirstOrderLimiter

Description
    Class with limiter function zero. Pure first order, use the 
    cell-center values directly on the face-center.
        
SourceFiles
    wPureFirstOrder.C

Authors
    Saleh Abuhanieh

\*---------------------------------------------------------------------------*/

#ifndef wPureFirstOrder_H
#define wPureFirstOrder_H

#include "vector.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class SuperBeeLimiter Declaration
\*---------------------------------------------------------------------------*/

template<class LimiterFunc>
class wPureFirstOrderLimiter
:
    public LimiterFunc
{

public:

    wPureFirstOrderLimiter(Istream&)
    {}

    scalar limiter
    (
        const scalar cdWeight,
        const scalar faceFlux,
        const typename LimiterFunc::phiType& phiP,
        const typename LimiterFunc::phiType& phiN,
        const typename LimiterFunc::gradPhiType& gradcP,
        const typename LimiterFunc::gradPhiType& gradcN,
        const vector& d
    ) const
    {
        //scalar r = LimiterFunc::r
        //(
            //faceFlux, phiP, phiN, gradcP, gradcN, d
        //);

        //scalar cd = (faceFlux >= 0 ? 1.0-cdWeight : cdWeight);
        return 0.0; //return max(max(min(1.0/cd*r, 1), min(r, 1.0/cd)), 0);
    }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
