/*---------------------------------------------------------------------------*\

    TAMS-AERO: Is part of the Turkish Aerospace Multiphysics Solver
               which is dedicated for the aerodynamic applications.
               
-------------------------------------------------------------------------------
License
    This file is part of TAMS-AERO.

\*---------------------------------------------------------------------------*/

#include "lusgsTams.H"
#include "emptyFvPatch.H"
#include "transformFvPatchFields.H"
#include "fvjOperators.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePreconditionerType(lusgsTams, 2, 1);
}


// ************************************************************************* //
