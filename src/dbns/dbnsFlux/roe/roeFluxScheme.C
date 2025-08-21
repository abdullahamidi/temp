/*---------------------------------------------------------------------------*\

    HiSA: High Speed Aerodynamic solver

    Copyright (C) 2014-2018 Oliver Oxtoby - CSIR, South Africa
    Copyright (C) 2014-2018 Johan Heyns - CSIR, South Africa
    Copyright (C) 1991-2008 OpenCFD Ltd.

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

\*---------------------------------------------------------------------------*/

#include "roeFluxScheme.H"
#include "addToRunTimeSelectionTable.H"
#include "Switch.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "surfaceInterpolate.H"
#include "linear.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

defineTypeNameAndDebug(roeFluxScheme, 0);
addToRunTimeSelectionTable(fluxScheme, roeFluxScheme, dictionary);


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

roeFluxScheme::roeFluxScheme
(
    const dictionary& dict,
    const fluidThermo& thermo,
    const volScalarField& rho,
    const volVectorField& U,
    const volVectorField& rhoU,
    const volScalarField& rhoE
)
:
    fluxScheme(typeName, dict),
    mesh_(U.mesh()),
    thermo_(thermo),
    rho_(rho),
    U_(U),
    rhoU_(rhoU),
    rhoE_(rhoE),

    pos_(surfaceScalarField
    (
        IOobject
        (
            "pos",
            mesh_.time().timeName(),
            mesh_
        ),
        mesh_,
        dimensionedScalar("pos", dimless, 1.0)
    )),

    neg_(surfaceScalarField
    (
        IOobject
        (
            "neg",
            mesh_.time().timeName(),
            mesh_
        ),
        mesh_,
        dimensionedScalar("neg", dimless, -1.0)
    ))

{}


// * * * * * * * * * * * * * * * * Destructors * * * * * * * * * * * * * * * //

roeFluxScheme::~roeFluxScheme()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::roeFluxScheme::calcFlux(surfaceScalarField& phi, surfaceVectorField& phiUp, surfaceScalarField& phiEp, surfaceVectorField& Up)
{

    surfaceVectorField n(mesh_.Sf()/mesh_.magSf());
    surfaceScalarField magSf(mesh_.magSf());
#if OPENFOAM >= 1712
    // Prevent oriented/unoriented incompatibility below
    n.setOriented(false);
#endif

    // Check flux relative to mesh movement
    if (mesh_.moving())
    {
        FatalErrorInFunction
            << "Roe does not support moving meshes for now! "
            << endl << endl
            << exit(FatalError);
    }

    // Left and right states
    tmp< surfaceScalarField > rho_l (fvc::interpolate(rho_, pos_, "reconstruct(rho)"));
    tmp< surfaceScalarField > rho_r (fvc::interpolate(rho_, neg_, "reconstruct(rho)"));

    tmp< volScalarField > p = thermo_.p();
    tmp< surfaceScalarField > p_l = fvc::interpolate(p(), pos_, "reconstruct(rho)");
    tmp< surfaceScalarField > p_r = fvc::interpolate(p(), neg_, "reconstruct(rho)");
    p.clear();

    surfaceVectorField U_l (IOobject("U_l", mesh_.time().timeName(), mesh_), mesh_, dimensionedVector("zero", U_.dimensions(), vector::zero));
    surfaceVectorField U_r (IOobject("U_r", mesh_.time().timeName(), mesh_), mesh_, dimensionedVector("zero", U_.dimensions(), vector::zero));
    U_l.replace(0, fvc::interpolate(U_.component(0), pos_, "reconstruct(U)"));
    U_l.replace(1, fvc::interpolate(U_.component(1), pos_, "reconstruct(U)"));
    U_l.replace(2, fvc::interpolate(U_.component(2), pos_, "reconstruct(U)"));
    U_r.replace(0, fvc::interpolate(U_.component(0), neg_, "reconstruct(U)"));
    U_r.replace(1, fvc::interpolate(U_.component(1), neg_, "reconstruct(U)"));
    U_r.replace(2, fvc::interpolate(U_.component(2), neg_, "reconstruct(U)"));

    // Acoustic velocity - c = sqrt(\gamma R T)
    tmp< volScalarField > gamma = thermo_.gamma();
    tmp< volScalarField > psi = thermo_.psi();
    dimensionedScalar c0("c0",dimVelocity,VSMALL);
    tmp< volScalarField > c = max(sqrt(gamma()/psi()),c0);
    tmp< surfaceScalarField > c_l = fvc::interpolate(c(), pos_, "reconstruct(T)");
    tmp< surfaceScalarField > c_r = fvc::interpolate(c(), neg_, "reconstruct(T)");
    c.clear();
    psi.clear();

    tmp< volScalarField > E(rhoE_/rho_);
    tmp< surfaceScalarField > E_l (fvc::interpolate(E(), pos_, "reconstruct(T)"));
    tmp< surfaceScalarField > E_r (fvc::interpolate(E(), neg_, "reconstruct(T)"));

    // NOTE: Literature suggest enthalpy should be interpolated seperately and
    // not be assembled using left and right states of energy and pressure
//    tmp< surfaceScalarField > H_l (E_l() + p_l()/rho_l());
//    tmp< surfaceScalarField > H_r (E_r() + p_r()/rho_r());
    tmp<volScalarField > H (E()+ p()/rho_);
    tmp<surfaceScalarField > H_l (fvc::interpolate(H(), pos_, "reconstruct(T)"));
    tmp<surfaceScalarField > H_r (fvc::interpolate(H(), neg_, "reconstruct(T)"));
    E.clear();
    H.clear();

    //// Step I: Compute the Roe average values

    dimensionedScalar rho0("rho0",dimDensity,VSMALL);
    tmp< surfaceScalarField > rhoAvg = sqrt(rho_l()*rho_r());
    tmp< surfaceScalarField > coefR = sqrt(max(rho0,rho_r())/max(rho0,rho_l()));
    tmp< surfaceVectorField > uAvg = (coefR()*U_r + U_l)/(coefR() + 1.0);
    tmp< surfaceScalarField > HAvg = (coefR()*H_r() + H_l())/(coefR() + 1.0);
    tmp< surfaceScalarField > cAvg = sqrt((fvc::interpolate(gamma())-1.0)*(HAvg() - 0.5*magSqr(uAvg())));
    gamma.clear();
    coefR.clear();

        // Contravariant velocity
    tmp< surfaceScalarField > uMag_l = U_l&n;
    tmp< surfaceScalarField > uMag_r = U_r&n;
    tmp< surfaceScalarField > uMagAvg = uAvg()&n;

        // Compute primitive differences
    const surfaceScalarField deltaP (p_r() - p_l());
    const surfaceScalarField deltaRho(rho_r() - rho_l());
    const surfaceVectorField deltaU(U_r - U_l);
    const surfaceScalarField deltaContrV(deltaU & n);

    //// Step II: Compute the eigenvalues

    tmp< surfaceScalarField > lambda1_tmp = mag(uMagAvg() - cAvg());
    surfaceScalarField lambda1 = lambda1_tmp.ref();
    tmp< surfaceScalarField > lambda2_tmp = mag(uMagAvg());
    surfaceScalarField lambda2 = lambda2_tmp.ref();
    tmp< surfaceScalarField > lambda3_tmp = mag(uMagAvg() + cAvg());
    surfaceScalarField lambda3 = lambda3_tmp.ref();

        // Harten entropy correction (according to S. Chun et al. 2013)
    tmp< surfaceScalarField > eps_tmp = 0.1*cAvg(); //adjustable parameter
    surfaceScalarField eps = eps_tmp.ref();
    forAll(eps, facei)
    {   scalar eps_ = 0.05*max(mag(lambda1[facei]), mag(lambda3[facei]));
        if (lambda1[facei] < eps_)
        {
            lambda1[facei] = (sqr(lambda1[facei]) + sqr(eps_))/(2.0*eps_);
        }
        if (lambda2[facei] < eps_)
        {
            lambda2[facei] = (sqr(lambda2[facei]) + sqr(eps_))/(2.0*eps_);
        }
        if (lambda3[facei] < eps_)
        {
            lambda3[facei] = (sqr(lambda3[facei]) + sqr(eps_))/(2.0*eps_);
        }
    }

    forAll(eps.boundaryField(), patchi)
    {
        const scalarField& epsf = eps.boundaryField()[patchi];
        scalarField& lambda1f = lambda1.boundaryFieldRef()[patchi];
        scalarField& lambda2f = lambda2.boundaryFieldRef()[patchi];
        scalarField& lambda3f = lambda3.boundaryFieldRef()[patchi];
        forAll(epsf, facei)
        { scalar eps_ = 0.05*max(mag(lambda1f[facei]), mag(lambda3f[facei]));
            if (lambda1f[facei] < eps_)
            {
                lambda1f[facei] = (sqr(lambda1f[facei]) + sqr(eps_))/(2.0*eps_);
            }
            if (lambda2f[facei] < eps_)
            {
                lambda2f[facei] = (sqr(lambda2f[facei]) + sqr(eps_))/(2.0*eps_);
            }
            if (lambda3f[facei] < eps_)
            {
                lambda3f[facei] = (sqr(lambda3f[facei]) + sqr(eps_))/(2.0*eps_);
            }
        }
    }

    //// Step III: Compute the right eigenvectors

        // rho row:
    const scalar l1rho = 1.0;
    const scalar l2rho = 1.0;
    const scalar l3rho = 0.0;
    const scalar l4rho = 1.0;

        // first U column
    const surfaceVectorField l1U(uAvg() - cAvg()*n);

        // second U column
    const surfaceVectorField l2U(uAvg());

        // third U column
    const surfaceVectorField l3U(deltaU - deltaContrV*n);

        // fourth U column
    const surfaceVectorField l4U(uAvg() + cAvg()*n);

        // E row
    const surfaceScalarField l1e(HAvg() - cAvg()*uMagAvg());
    const surfaceScalarField l2e(0.5*(uAvg() & uAvg()));
    const surfaceScalarField l3e((uAvg() & deltaU) - uMagAvg()*deltaContrV);
    const surfaceScalarField l4e(HAvg() + cAvg()*uMagAvg());

    //// Step IV: Compute the wave strengths

    const surfaceScalarField r1((deltaP - rhoAvg()*cAvg()*deltaContrV)/(2.0*sqr(cAvg())));
    const surfaceScalarField r2(deltaRho - deltaP/sqr(cAvg()));
    const surfaceScalarField r3((deltaP + rhoAvg()*cAvg()*deltaContrV)/(2.0*sqr(cAvg())));

    //// Step V: Use all of the above quantities to compute F_{i+1/2}

        // Compute flux differences

            // Components of deltaF1
    const surfaceScalarField diffF11(lambda1*r1*l1rho);
    const surfaceVectorField diffF124(lambda1*r1*l1U);
    const surfaceScalarField diffF15(lambda1*r1*l1e);

            // Components of deltaF2
    const surfaceScalarField diffF21(lambda2*(r2*l2rho + rhoAvg()*l3rho));
    const surfaceVectorField diffF224(lambda2*(r2*l2U + rhoAvg()*l3U));
    const surfaceScalarField diffF25(lambda2*(r2*l2e + rhoAvg()*l3e));

            // Components of deltaF3
    const surfaceScalarField diffF31(lambda3*r3*l4rho);
    const surfaceVectorField diffF324(lambda3*r3*l4U);
    const surfaceScalarField diffF35(lambda3*r3*l4e);

        // Compute left and right fluxes

            // Left flux 5-vector
    const surfaceScalarField fluxLeft11(rho_l()*uMag_l());
    const surfaceVectorField fluxLeft124(U_l*fluxLeft11 + n*p_l());
    const surfaceScalarField fluxLeft15(H_l()*fluxLeft11);

            // Right flux 5-vector
    const surfaceScalarField fluxRight11(rho_r()*uMag_r());
    const surfaceVectorField fluxRight124(U_r*fluxRight11 + n*p_r());
    const surfaceScalarField fluxRight15(H_r()*fluxRight11);

        // Compute face flux 5-vector

    const surfaceScalarField flux1
        (0.5*(fluxLeft11 + fluxRight11 - (diffF11 + diffF21 + diffF31)));

    const surfaceVectorField flux24
        (0.5*(fluxLeft124 + fluxRight124 - (diffF124 + diffF224 + diffF324)));

    const surfaceScalarField flux5
        (0.5*(fluxLeft15 + fluxRight15 - (diffF15 + diffF25 + diffF35)));

            // Continuity

    phi =  (flux1)*magSf;

            // Momentum

    phiUp = (flux24)*magSf;

            // Energy

    phiEp = (flux5)*magSf;

            // Face velocity for sigmaDotU (Viscous flow)

    Up = linearInterpolate(U_)*magSf;

}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
