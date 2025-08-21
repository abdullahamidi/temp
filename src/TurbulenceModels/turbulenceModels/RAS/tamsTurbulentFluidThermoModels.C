/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2013-2016 OpenFOAM Foundation
    Copyright (C) 2022 OpenCFD Ltd.
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "turbulentFluidThermoModels.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


// -------------------------------------------------------------------------- //
// Laminar models
// -------------------------------------------------------------------------- //

//#include "Stokes.H"
//makeLaminarModel(Stokes);

//#include "generalizedNewtonian.H"
//makeLaminarModel(generalizedNewtonian);


// -------------------------------------------------------------------------- //
// RAS models
// -------------------------------------------------------------------------- //

//#include "SpalartAllmaras.H"
//makeRASModel(SpalartAllmaras);

//#include "kEpsilon.H"
//makeRASModel(kEpsilon);

#include "kOmegaSSTCC.H"
makeRASModel(kOmegaSSTCC);

#include "kOmegaSST_tamsOverset.H"
makeRASModel(kOmegaSST_tamsOverset);

#include "kOmegaSSTLM_tamsOverset.H"
makeRASModel(kOmegaSSTLM_tamsOverset);

#include "NASA_SST2003.H"
makeRASModel(NASA_SST2003);

#include "NASA_SST2003m.H"
makeRASModel(NASA_SST2003m);

#include "NASA_SST2003Vm.H"
makeRASModel(NASA_SST2003Vm);

#include "SpalartAllmarasRCC.H"
makeRASModel(SpalartAllmarasRCC);

// -------------------------------------------------------------------------- //
// LES models
// -------------------------------------------------------------------------- //

//#include "Smagorinsky.H"
//makeLESModel(Smagorinsky);

//#include "WALE.H"
//makeLESModel(WALE);


// ************************************************************************* //
