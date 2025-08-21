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

#include "propulsionExternalForce.H"
#include "addToRunTimeSelectionTable.H"
#include "sixDoFRigidBodyMotion.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace sixDoFRigidBodyMotionRestraints
{
    defineTypeNameAndDebug(propulsionExternalForce, 0);

    addToRunTimeSelectionTable
    (
        sixDoFRigidBodyMotionRestraint,
        propulsionExternalForce,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::sixDoFRigidBodyMotionRestraints::propulsionExternalForce::propulsionExternalForce
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

Foam::sixDoFRigidBodyMotionRestraints::propulsionExternalForce::~propulsionExternalForce()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::sixDoFRigidBodyMotionRestraints::propulsionExternalForce::restrain
(
    const sixDoFRigidBodyMotion& motion,
    vector& restraintPosition,
    vector& restraintForce,
    vector& restraintMoment
) const
{

    const scalar t = motion.time().value();
    
    restraintPosition = motion.transform(location_); //move the point (where we apply the force) with the body

    bool forceOff = t < forceStartTime_ || t >= forceEndTime_;
    if (forceOff)
    {
        restraintForce = Zero;
        restraintMoment = Zero; 
    }
    else
    {
        restraintForce = motion.transform(force_->value(t)); //move the force with the body
        restraintMoment = (restraintPosition - motion.centreOfMass()) ^ restraintForce;
    }    

    if (motion.report())
    {
        Info<< " propulsion force location " << restraintPosition
            << " propulsion force " << restraintForce
            << " propulsion moment " << restraintMoment
            //<< " is the force active? " <<std::noboolalpha<< forceOff            
            << endl;
    }
}


bool Foam::sixDoFRigidBodyMotionRestraints::propulsionExternalForce::read
(
    const dictionary& sDoFRBMRDict
)
{
    sixDoFRigidBodyMotionRestraint::read(sDoFRBMRDict);

    sDoFRBMRCoeffs_.readEntry("location", location_);
    force_ = Function1<vector>::New("force", sDoFRBMRDict);
    sDoFRBMRCoeffs_.readEntry("forceStartTime", forceStartTime_);
    sDoFRBMRCoeffs_.readEntry("forceEndTime", forceEndTime_);

    if (forceStartTime_ < 0 || forceEndTime_ < 0)
    {
        FatalIOErrorInFunction(sDoFRBMRDict)
            << "Force start time and end time shall be positive scalars." << exit(FatalIOError);
    }
    if (forceStartTime_ > forceEndTime_)
    {
        FatalIOErrorInFunction(sDoFRBMRDict)
            << "Force end time can't be less than the force start time." << exit(FatalIOError);
    }    
    

    return true;
}


void Foam::sixDoFRigidBodyMotionRestraints::propulsionExternalForce::write
(
    Ostream& os
) const
{
    os.writeEntry("location", location_);
    force_().writeData(os);
    os.writeEntry("forceStartTime", forceStartTime_);
    os.writeEntry("forceEndTime", forceEndTime_);
}

// ************************************************************************* //
