/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2019-2021 OpenCFD Ltd.
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

#include "SpalartAllmarasRCC.H"
#include "addToRunTimeSelectionTable.H"
#include "fvOptions.H"
#include "bound.H"
#include "wallDist.H"
//// S.A {01.12.2022}
#include "IOMRFZoneList.H"
//// S.A {01.12.2022}

#ifndef SpalartAllmarasRCC_C
#define SpalartAllmarasRCC_C

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace RASModels
{

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class BasicTurbulenceModel>
tmp<volScalarField> SpalartAllmarasRCC<BasicTurbulenceModel>::chi() const
{
    return nuTilda_/this->nu();
}


template<class BasicTurbulenceModel>
tmp<volScalarField> SpalartAllmarasRCC<BasicTurbulenceModel>::fv1
(
    const volScalarField& chi
) const
{
    const volScalarField chi3(pow3(chi));

    return chi3/(chi3 + pow3(Cv1_));
}


template<class BasicTurbulenceModel>
tmp<volScalarField::Internal> SpalartAllmarasRCC<BasicTurbulenceModel>::fv2
(
    const volScalarField::Internal& chi,
    const volScalarField::Internal& fv1
) const
{
    return scalar(1) - chi/(scalar(1) + chi*fv1);
}


template<class BasicTurbulenceModel>
tmp<volScalarField::Internal> SpalartAllmarasRCC<BasicTurbulenceModel>::Stilda()
const
{
    const volScalarField chi(this->chi());

    const volScalarField fv1(this->fv1(chi));

    const volScalarField::Internal Omega
    (
        ::sqrt(scalar(2))*mag(skew(fvc::grad(this->U_)().v()))
    );

    return
    (
        max
        (
            Omega + fv2(chi(), fv1())*nuTilda_()/sqr(kappa_*y_),
            Cs_*Omega
        )
    );
}


template<class BasicTurbulenceModel>
tmp<volScalarField::Internal> SpalartAllmarasRCC<BasicTurbulenceModel>::fw
(
    const volScalarField::Internal& Stilda
) const
{
    const volScalarField::Internal r
    (
        min
        (
            nuTilda_
           /(
               max
               (
                   Stilda,
                   dimensionedScalar(Stilda.dimensions(), SMALL)
               )
              *sqr(kappa_*y_)
            ),
            scalar(10)
        )
    );

    const volScalarField::Internal g(r + Cw2_*(pow6(r) - r));

    return
        g*pow
        (
            (scalar(1) + pow6(Cw3_))/(pow6(g) + pow6(Cw3_)),
            scalar(1)/scalar(6)
        );
}


template<class BasicTurbulenceModel>
void SpalartAllmarasRCC<BasicTurbulenceModel>::correctNut()
{
    this->nut_ = nuTilda_*this->fv1(this->chi());
    this->nut_.correctBoundaryConditions();
    fv::options::New(this->mesh_).correct(this->nut_);

    BasicTurbulenceModel::correctNut();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
SpalartAllmarasRCC<BasicTurbulenceModel>::SpalartAllmarasRCC
(
    const alphaField& alpha,
    const rhoField& rho,
    const volVectorField& U,
    const surfaceScalarField& alphaRhoPhi,
    const surfaceScalarField& phi,
    const transportModel& transport,
    const word& propertiesName,
    const word& type
)
:
    eddyViscosity<RASModel<BasicTurbulenceModel>>
    (
        type,
        alpha,
        rho,
        U,
        alphaRhoPhi,
        phi,
        transport,
        propertiesName
    ),

    sigmaNut_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "sigmaNut",
            this->coeffDict_,
            scalar(2)/scalar(3)
        )
    ),
    kappa_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "kappa",
            this->coeffDict_,
            0.41
        )
    ),

    Cb1_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "Cb1",
            this->coeffDict_,
            0.1355
        )
    ),
    Cb2_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "Cb2",
            this->coeffDict_,
            0.622
        )
    ),
    Cw1_(Cb1_/sqr(kappa_) + (scalar(1) + Cb2_)/sigmaNut_),
    Cw2_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "Cw2",
            this->coeffDict_,
            0.3
        )
    ),
    Cw3_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "Cw3",
            this->coeffDict_,
            2.0
        )
    ),
    Cv1_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "Cv1",
            this->coeffDict_,
            7.1
        )
    ),
    Cs_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "Cs",
            this->coeffDict_,
            0.3
        )
    ),
//// S.A {01.12.2022}
    cr1_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "cr1",
            this->coeffDict_,
            1.0
        )
    ),
    cr2_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "cr2",
            this->coeffDict_,
            12.0
        )
    ),
    cr3_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "cr3",
            this->coeffDict_,
            1.0
        )
    ),
    omga_(Zero),
//// S.A {01.12.2022}

    nuTilda_
    (
        IOobject
        (
            "nuTilda",
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_
    ),

    y_(wallDist::New(this->mesh_).y())
{
    if (type == typeName)
    {
        this->printCoeffs(type);
    }

//// S.A {01.12.2022}, TODO: include the dynamic mesh situation
    const auto* MRFZones =
        this->mesh_.objectRegistry::template cfindObject<IOMRFZoneList>("MRFProperties");
    const auto& mrf = MRFZones->MRFZoneList::getFromName("MRF1");
    omga_ = mrf.Omega();
    omga_.dimensions()[2] = -1; //correct the dimension to be 1/s
    //Info<<"omga_: "<<omga_<<endl;

//// S.A {01.12.2022}

}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
bool SpalartAllmarasRCC<BasicTurbulenceModel>::read()
{
    if (eddyViscosity<RASModel<BasicTurbulenceModel>>::read())
    {
        sigmaNut_.readIfPresent(this->coeffDict());
        kappa_.readIfPresent(this->coeffDict());

        Cb1_.readIfPresent(this->coeffDict());
        Cb2_.readIfPresent(this->coeffDict());
        Cw1_ = Cb1_/sqr(kappa_) + (scalar(1) + Cb2_)/sigmaNut_;
        Cw2_.readIfPresent(this->coeffDict());
        Cw3_.readIfPresent(this->coeffDict());
        Cv1_.readIfPresent(this->coeffDict());
        Cs_.readIfPresent(this->coeffDict());
//// S.A {01.12.2022}
	cr1_.readIfPresent(this->coeffDict());
        cr2_.readIfPresent(this->coeffDict());
        cr3_.readIfPresent(this->coeffDict());
//// S.A {01.12.2022}

        return true;
    }

    return false;
}


template<class BasicTurbulenceModel>
tmp<volScalarField> SpalartAllmarasRCC<BasicTurbulenceModel>::DnuTildaEff() const
{
    return tmp<volScalarField>::New
    (
        "DnuTildaEff",
        (nuTilda_ + this->nu())/sigmaNut_
    );
}


template<class BasicTurbulenceModel>
tmp<volScalarField> SpalartAllmarasRCC<BasicTurbulenceModel>::k() const
{
    // (B:Eq. 4.50)
    const scalar Cmu = 0.09;

    return tmp<volScalarField>::New
    (
        IOobject
        (
            IOobject::groupName("k", this->alphaRhoPhi_.group()),
            this->runTime_.timeName(),
            this->mesh_
        ),
        cbrt(this->fv1(this->chi()))
        *nuTilda_
        *::sqrt(scalar(2)/Cmu)
        *mag(symm(fvc::grad(this->U_))),
        this->nut_.boundaryField().types()
    );
}


template<class BasicTurbulenceModel>
tmp<volScalarField> SpalartAllmarasRCC<BasicTurbulenceModel>::epsilon() const
{
    // (B:Eq. 4.50)
    const scalar Cmu = 0.09;
    const dimensionedScalar nutSMALL(sqr(dimLength)/dimTime, SMALL);

    return tmp<volScalarField>::New
    (
        IOobject
        (
            IOobject::groupName("epsilon", this->alphaRhoPhi_.group()),
            this->runTime_.timeName(),
            this->mesh_
        ),
        pow(this->fv1(this->chi()), 0.5)
        *pow(::sqrt(Cmu)*this->k(), 2)
        /(nuTilda_ + this->nut_ + nutSMALL),
        this->nut_.boundaryField().types()
    );
}


template<class BasicTurbulenceModel>
tmp<volScalarField> SpalartAllmarasRCC<BasicTurbulenceModel>::omega() const
{
    // (P:p. 384)
    const scalar betaStar = 0.09;
    const dimensionedScalar k0(sqr(dimLength/dimTime), SMALL);

    return tmp<volScalarField>::New
    (
        IOobject
        (
            IOobject::groupName("omega", this->alphaRhoPhi_.group()),
            this->runTime_.timeName(),
            this->mesh_
        ),
        this->epsilon()/(betaStar*(this->k() + k0)),
        this->nut_.boundaryField().types()
    );
}


template<class BasicTurbulenceModel>
void SpalartAllmarasRCC<BasicTurbulenceModel>::correct()
{
    if (!this->turbulence_)
    {
        return;
    }

    {
        // Construct local convenience references
        const alphaField& alpha = this->alpha_;
        const rhoField& rho = this->rho_;
        const surfaceScalarField& alphaRhoPhi = this->alphaRhoPhi_;
        fv::options& fvOptions(fv::options::New(this->mesh_));

//// S.A {01.12.2022}
    Info<<"rotation correction starts"<<endl;
    tmp<volTensorField> tgradU = fvc::grad(this->U_);
    tmp<volSymmTensorField> S_ij = symm(tgradU()); //mean strain tensor
    volScalarField Ssquare(2.0*S_ij() && S_ij()); //S^(2)
    tmp<volTensorField> Omega_ij = skew(tgradU()); //vorticity tensor

    // create the Levi-Civita tensor in three separate Rank-2 tensors
    tensor Eps0(0.0, 0.0,  0.0,  0.0, 0.0, 1.0, 0.0, -1.0, 0.0);
    tensor Eps1(0.0, 0.0, -1.0,  0.0, 0.0, 0.0, 1.0,  0.0, 0.0);
    tensor Eps2(0.0, 1.0,  0.0, -1.0, 0.0, 0.0, 0.0,  0.0, 0.0);
    //Info<<"Eps0: "<<Eps0<<endl;

    Omega_ij.ref() += 2.0*(Eps0*omga_[0] + Eps1*omga_[1] + Eps2*omga_[2]);

    volScalarField Omegasquare(2.0*Omega_ij() && Omega_ij()); //Omega^(2)
    volScalarField OmegaMag(sqrt(2.0*Omega_ij() && Omega_ij())); //Omega tensor magnitude
    volScalarField rStar(sqrt(Ssquare)/OmegaMag); //r^{*}

    volScalarField D(sqrt(0.5*(Ssquare + Omegasquare)));
    tmp<volSymmTensorField> DS =
    (
        fvc::ddt(alpha*rho, S_ij())
       +fvc::div
        (
            alphaRhoPhi, S_ij()
        )
    );

    DS.ref()[0][0] += 2.0*( S_ij()[0][2]*omga_[1].value() - S_ij()[0][1]*omga_[2].value());
    DS.ref()[0][1] += 2.0*( S_ij()[1][2]*omga_[1].value() - S_ij()[1][1]*omga_[2].value());
    DS.ref()[0][2] += 2.0*( S_ij()[2][2]*omga_[1].value() - S_ij()[2][1]*omga_[2].value());

    DS.ref()[1][0] += 2.0*(-S_ij()[0][2]*omga_[0].value() + S_ij()[0][0]*omga_[2].value());
    DS.ref()[1][1] += 2.0*(-S_ij()[1][2]*omga_[0].value() + S_ij()[1][0]*omga_[2].value());
    DS.ref()[1][2] += 2.0*(-S_ij()[2][2]*omga_[0].value() + S_ij()[2][0]*omga_[2].value());

    DS.ref()[0][0] += 2.0*( S_ij()[0][1]*omga_[0].value() - S_ij()[0][0]*omga_[1].value());
    DS.ref()[0][1] += 2.0*( S_ij()[1][1]*omga_[0].value() - S_ij()[1][0]*omga_[1].value());
    DS.ref()[0][2] += 2.0*( S_ij()[2][1]*omga_[0].value() - S_ij()[2][0]*omga_[1].value());

    //volScalarField rTilda((Omega_ij().T() & S_ij()) && DS()); //r^{~}
    tmp<volSymmTensorField> twoOmgeaS_ij = symm(Omega_ij().T() & S_ij());
    for(int i=0; i<3; i++)
    {
        for(int j=0; j<3; j++)
        {
            twoOmgeaS_ij.ref()[i][j] = 2.0*(Omega_ij()[i][0] * S_ij()[j][0] + Omega_ij()[i][1] * S_ij()[j][1] + Omega_ij()[i][2] * S_ij()[j][2]);
        }
    }
    volScalarField rTilda(twoOmgeaS_ij && DS()); //r^{~}
    rTilda /= D*D*D*D;
    rTilda.dimensions()[0] = 0; //correct the dimension to be dimenstionless //// Result.dimensions().reset(dimensionSet(0, 1, -1, 0, 0, 0, 0));
    rTilda.dimensions()[1] = 0; //-
    //Info<<"rTilda: "<<rTilda.dimensions()<<endl;

    volScalarField Fr1
    (
        (1.0+cr1_)*2.0*rStar/(1.0+rStar)*(1.0-cr3_*atan(cr2_*rTilda))-cr1_
    );

    Info<<"min(Fr1): "<<min(Fr1).value()<<" max(Fr1): "<<max(Fr1).value()<<" average(Fr1): "<<average(Fr1).value()<<endl;
    Info<<"rotation correction ends"<<endl;

//// S.A {01.12.2022}

        eddyViscosity<RASModel<BasicTurbulenceModel>>::correct();

        const volScalarField::Internal Stilda(this->Stilda());

        tmp<fvScalarMatrix> nuTildaEqn
        (
            fvm::ddt(alpha, rho, nuTilda_)
          + fvm::div(alphaRhoPhi, nuTilda_)
          - fvm::laplacian(alpha*rho*DnuTildaEff(), nuTilda_)
          - Cb2_/sigmaNut_*alpha*rho*magSqr(fvc::grad(nuTilda_))
         ==
            Cb1_*alpha()*rho()*Stilda*nuTilda_()*Fr1() //Cb1_*alpha()*rho()*Stilda*nuTilda_() //// S.A {01.12.2022}
          - fvm::Sp(Cw1_*alpha()*rho()*fw(Stilda)*nuTilda_()/sqr(y_), nuTilda_)
          + fvOptions(alpha, rho, nuTilda_)
        );

        nuTildaEqn.ref().relax();
        fvOptions.constrain(nuTildaEqn.ref());
        solve(nuTildaEqn);
        fvOptions.correct(nuTilda_);
        bound(nuTilda_, dimensionedScalar(nuTilda_.dimensions(), Zero));
        nuTilda_.correctBoundaryConditions();
    }

    // Update nut with latest available nuTilda
    correctNut();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace Foam

#endif

// ************************************************************************* //
