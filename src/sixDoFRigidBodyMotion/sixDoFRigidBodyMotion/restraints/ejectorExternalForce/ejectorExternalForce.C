/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
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

#include "ejectorExternalForce.H"
#include "addToRunTimeSelectionTable.H"
#include "sixDoFRigidBodyMotion.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace sixDoFRigidBodyMotionRestraints
{
    defineTypeNameAndDebug(ejectorExternalForce, 0);

    addToRunTimeSelectionTable
    (
        sixDoFRigidBodyMotionRestraint,
        ejectorExternalForce,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::sixDoFRigidBodyMotionRestraints::ejectorExternalForce::ejectorExternalForce
(
    const word& name,
    const dictionary& sDoFRBMRDict
)
:
    sixDoFRigidBodyMotionRestraint(name, sDoFRBMRDict)
{
    read(sDoFRBMRDict);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::sixDoFRigidBodyMotionRestraints::ejectorExternalForce::~ejectorExternalForce()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::sixDoFRigidBodyMotionRestraints::ejectorExternalForce::restrain
(
    const sixDoFRigidBodyMotion& motion,
    vector& restraintPosition,
    vector& restraintForce,
    vector& restraintMoment
) const
{
    restraintPosition = motion.transform(location_);

    restraintForce = force_;

    restraintMoment = Zero;    //restraintMoment = location_ ^ force_;

    vector displacement = motion.centreOfMass() - motion.initialCentreOfMass();

    scalar dispMagInForceDir = (displacement & direction_); 

    if (dispMagInForceDir > strokeLength_)
    {
        restraintForce = Zero;
        //restraintMoment = Zero; 
    }

    if (motion.report())
    {
        Info<< " force location " << restraintPosition
            << " force " << restraintForce
            //<< " moment " << restraintMoment
            << " dispMagInForceDir " << dispMagInForceDir
            << " is the force still applied ? " << (dispMagInForceDir < strokeLength_)
            << endl;
    }
}


bool Foam::sixDoFRigidBodyMotionRestraints::ejectorExternalForce::read
(
    const dictionary& sDoFRBMRDict
)
{
    sixDoFRigidBodyMotionRestraint::read(sDoFRBMRDict);

    sDoFRBMRCoeffs_.readEntry("location", location_);
    sDoFRBMRCoeffs_.readEntry("force", force_);
    sDoFRBMRCoeffs_.readEntry("strokeLength", strokeLength_);
    sDoFRBMRCoeffs_.readEntry("direction", direction_);

    return true;
}


void Foam::sixDoFRigidBodyMotionRestraints::ejectorExternalForce::write
(
    Ostream& os
) const
{
    os.writeEntry("location", location_);
    os.writeEntry("force", force_);
    os.writeEntry("strokeLength", strokeLength_);
    os.writeEntry("direction", direction_);
}

// ************************************************************************* //
