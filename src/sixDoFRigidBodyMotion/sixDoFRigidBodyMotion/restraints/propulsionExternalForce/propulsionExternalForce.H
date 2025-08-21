/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2013 OpenFOAM Foundation
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

Class
    Foam::sixDoFRigidBodyMotionRestraints::propulsionExternalForce

Description

    For applying a propulsion force to a moving body in the CFD/global frame
    while moving according to the 6dof equations.
    The value of forces can be function of time or constant, the forces
    is applied at the force start time and stopes and the propulsion end time.

    S.A 26.09.24

Usage
    Example applying a propulsion force to a certain region:
    \verbatim
    restraints
    {
        force_propulsion
        {
            sixDoFRigidBodyMotionRestraint propulsionExternalForce;
            location          (4.2614 3.302 -1.03478);
            force             (4500 0 0);
            forceStartTime    0.06;
            forceEndTime      1.0;;
        }
    }
    \endverbatim    

SourceFiles
    propulsionExternalForce.C

\*---------------------------------------------------------------------------*/

#ifndef propulsionExternalForce_H
#define propulsionExternalForce_H

#include "sixDoFRigidBodyMotionRestraint.H"
#include "point.H"
#include "Function1.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace sixDoFRigidBodyMotionRestraints
{

/*---------------------------------------------------------------------------*\
                          Class linearSpring Declaration
\*---------------------------------------------------------------------------*/

class propulsionExternalForce
:
    public sixDoFRigidBodyMotionRestraint
{
    // Private data
        //- Point of application of the force (m)
        vector location_;

        //- Force vector (N)
        //vector force_;
        autoPtr<Function1<vector>> force_;

        //- Force start time (s)
        scalar forceStartTime_;

        //- Force end time (s)
        scalar forceEndTime_;

public:

    //- Runtime type information
    TypeName("propulsionExternalForce");


    // Constructors

        //- Construct from components
        propulsionExternalForce
        (
            const word& name,
            const dictionary& sDoFRBMRDict
        );

        //- Construct and return a clone
        virtual autoPtr<sixDoFRigidBodyMotionRestraint> clone() const
        {
            return autoPtr<sixDoFRigidBodyMotionRestraint>
            (
                new propulsionExternalForce(*this)
            );
        }


    //- Destructor
    virtual ~propulsionExternalForce();


    // Member Functions

        //- Calculate the restraint position, force and moment.
        //  Global reference frame vectors.
        virtual void restrain
        (
            const sixDoFRigidBodyMotion& motion,
            vector& restraintPosition,
            vector& restraintForce,
            vector& restraintMoment
        ) const;

        //- Update properties from given dictionary
        virtual bool read(const dictionary& sDoFRBMRCoeff);

        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace solidBodyMotionFunctions
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
