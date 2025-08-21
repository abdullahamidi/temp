/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
     \\/     M anipulation  |
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



    Saleh Abuhanieh 23.10.21
    - I started from the code written by Ajit Kumar (Shiv Nadar University)
    - Updated to OpenFOAM 2006
    - Corrected the implementaion to match the required usual/standard 6DOF function
      The original implementaion was not providing the correct septernion representation
      Thus, the motion was not correct.
    - Adding the restart capability, it was restarting from the original state only {10.11.21}

\*---------------------------------------------------------------------------*/

#include "sixDoFDynamicMotion.H"
#include "addToRunTimeSelectionTable.H"

#include "polyMesh.H"
#include "pointPatchDist.H"
#include "pointConstraints.H"
#include "uniformDimensionedFields.H"
#include "forces.H"
#include "mathematicalConstants.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace solidBodyMotionFunctions
{
    defineTypeNameAndDebug(sixDoFDynamicMotion, 0);
    addToRunTimeSelectionTable
    (
        solidBodyMotionFunction,
        sixDoFDynamicMotion,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::solidBodyMotionFunctions::sixDoFDynamicMotion::sixDoFDynamicMotion
(
    const dictionary& SBMFCoeffs,
    const Time& runTime
)
:
    solidBodyMotionFunction(SBMFCoeffs, runTime),
    dict_(SBMFCoeffs),
    //motion_(SBMFCoeffs,SBMFCoeffs, runTime),
    motion_
    (
        SBMFCoeffs,
        IOobject
        (
            "sixDoFRigidBodyMotionState",
            time_.timeName(),
            "uniform",
            time_
        ).typeHeaderOk<IOdictionary>(true)
      ? IOdictionary
        (
            IOobject
            (
                "sixDoFRigidBodyMotionState",
                time_.timeName(),
                "uniform",
                time_,
                IOobject::READ_IF_PRESENT,
                IOobject::NO_WRITE,
                false
            )
        )
      : SBMFCoeffs,
        time_
    ),
    curTimeIndex_(-1),
    patches_(wordRes(SBMFCoeffs.lookup("patches"))),
    rhoInf_(1.0),
    pRef_(SBMFCoeffs.lookupOrDefault<scalar>("pRef", 0.0)),
    rhoName_(SBMFCoeffs.lookupOrDefault<word>("rho", "rho")),
    test_(SBMFCoeffs.getOrDefault("test", false)),
    CofGvelocity_(SBMFCoeffs.getOrDefault<word>("CofGvelocity", "none")),
    CoM0_(dict_.lookup("centreOfMass"))
{
    read(SBMFCoeffs);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::solidBodyMotionFunctions::sixDoFDynamicMotion::~sixDoFDynamicMotion()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::septernion
Foam::solidBodyMotionFunctions::sixDoFDynamicMotion::transformation() const
{
    

    // Store the motion state at the beginning of the time-step
    bool firstIter = false;
    if (curTimeIndex_ != time_.time().timeIndex())
    {
        motion_.newTime();
        curTimeIndex_ = time_.time().timeIndex();
        firstIter = true;
    }

    dimensionedVector g("g", dimAcceleration, Zero);
    if (dict_.found("g"))
    {
        dict_.lookup("g") >> g;
    }
    //Info << "g = " << g << endl;

    // const scalar ramp = min(max((this->db().time().value() - 5)/10, 0), 1);
    const scalar ramp = 1.0;

    if (test_)
    {
        motion_.update
        (
            firstIter,
            ramp*(motion_.mass()*g.value()),
            ramp*(motion_.mass()*(motion_.momentArm() ^ g.value())),
            time_.deltaTValue(),
            time_.deltaT0Value()
        );
    }
    else
    {
        dictionary forcesDict;

        forcesDict.add("type", functionObjects::forces::typeName);
        forcesDict.add("patches", patches_);
        forcesDict.add("rhoInf", rhoInf_);
        forcesDict.add("pRef", pRef_);
        forcesDict.add("rho", rhoName_);
        forcesDict.add("CofR", motion_.centreOfRotation());

        functionObjects::forces f("forces", time_, forcesDict);

        f.calcForcesMoments();

        motion_.update
        (
            firstIter,
            ramp*(f.forceEff() + motion_.mass()*g.value()),
            ramp
           *(
               f.momentEff()
             + motion_.mass()*(motion_.momentArm() ^ g.value())
            ),
            time_.deltaTValue(),
            time_.deltaT0Value()
        );

        if (CofGvelocity_ != "none")
        {
            if
            (
                time_.time().foundObject<uniformDimensionedVectorField>
                (
                    CofGvelocity_
                )
            )
            {
                uniformDimensionedVectorField& vel =
                    time_.time().lookupObjectRef<uniformDimensionedVectorField>
                    (
                        CofGvelocity_
                    );
                vel = motion_.v();
            }
        }
    }

    vector CoM = motion_.state().centreOfRotation();
    quaternion R(motion_.state().Q());
    septernion TR(septernion(-CoM0_- (CoM - CoM0_))*R*septernion(CoM0_));

    DebugInFunction << "Time = " << time_.value() << " transformation: " << TR << endl;

    //printing the motion state for restart
    if(time_.writeTime())
    {
        IOdictionary dict
        (
            IOobject
            (
                "sixDoFRigidBodyMotionState",
                time_.timeName(),
                "uniform",
                time_,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            )
        );

        motion_.state().write(dict);
        dict.regIOobject::write();
    }


    return TR;
}


bool Foam::solidBodyMotionFunctions::sixDoFDynamicMotion::read
(
    const dictionary& SBMFCoeffs
)
{
    solidBodyMotionFunction::read(SBMFCoeffs);

    return true;
}


// ************************************************************************* //
