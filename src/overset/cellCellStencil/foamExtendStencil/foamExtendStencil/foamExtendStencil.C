/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2017-2019 OpenCFD Ltd.
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify i
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

#include "foamExtendStencil.H"
#include "addToRunTimeSelectionTable.H"
#include "OBJstream.H"
#include "Time.H"
#include "fvMeshSubset.H"

#include "globalIndex.H"
#include "oversetFvPatch.H"
#include "zeroGradientFvPatchFields.H"
#include "syncTools.H"
#include "regionSplit.H"
#include "dynamicOversetFvMesh.H"
// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace cellCellStencils
{
    defineTypeNameAndDebug(foamExtendStencil, 0);
}
}

// the ESI implementation
void Foam::cellCellStencils::foamExtendStencil::extendCellStencil
(
    const globalIndex& globalCells
)
{
    if (debug)
    {
        Info<<"\nExtend the basic cellStencil\n"<<endl;
    }
    const Time& runTime = mesh_.time();
    scalar timeBefore = runTime.elapsedCpuTime();

    // Send cell centre back to donor
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // The complication is that multiple acceptors need the same donor
    // (but with different weights obviously)
    // So we do multi-pass:
    // - send over cc of acceptor for which we want stencil.
    //   Consistently choose the acceptor with smallest magSqr in case of
    //   multiple acceptors for the containing cell/donor.
    // - find the cell-cells and weights for the donor
    // - send back together with the acceptor cc
    // - use the acceptor cc to see if it was 'me' that sent it. If so
    //   mark me as complete so it doesn't get involved in the next loop.
    // - loop until all resolved.

    // Special value for unused points
    const vector greatPoint(GREAT, GREAT, GREAT);

    boolList isValidDonor(mesh_.nCells(), true);
    forAll(cellTypes_, celli)
    {
        if (cellTypes_[celli] == HOLE)
        {
            isValidDonor[celli] = false;
        }
    }


    // Has acceptor been handled already?
    bitSet doneAcceptor(interpolationCells_.size());

    while (true)
    {
        pointField samples(cellInterpolationMap().constructSize(), greatPoint);

        // Fill remote slots (override old content). We'll find out later
        // on which one has won and mark this one in doneAcceptor.
        label nSamples = 0;
        forAll(interpolationCells_, i)
        {
            if (!doneAcceptor[i])
            {
                label cellI = interpolationCells_[i];
                const point& cc = mesh_.cellCentres()[cellI];
                const labelList& slots = cellStencil_[cellI];

                if (slots.size() != 1)
                {
                    FatalErrorInFunction<< "Problem:" << slots
                        << abort(FatalError);
                }

                forAll(slots, slotI)
                {
                    label elemI = slots[slotI];
                    //Pout<< "    acceptor:" << cellI
                    //    << " at:" << mesh_.cellCentres()[cellI]
                    //    << " global:" << globalCells.toGlobal(cellI)
                    //    << " found in donor:" << elemI << endl;
                    minMagSqrEqOp<point>()(samples[elemI], cc);
                }
                nSamples++;
            }
        }


        if (returnReduce(nSamples, sumOp<label>()) == 0)
        {
            break;
        }

        // Send back to donor. Make sure valid point takes priority
        mapDistributeBase::distribute<point, minMagSqrEqOp<point>, flipOp>
        (
            Pstream::commsTypes::nonBlocking,
            List<labelPair>::null(), //List<labelPair>(), //S.A {13.09.22}
            mesh_.nCells(),
            cellInterpolationMap().constructMap(),
            false,
            cellInterpolationMap().subMap(),
            false,
            samples,
            greatPoint,                             // nullValue  //S.A {13.09.22}
            minMagSqrEqOp<point>(),
            flipOp(),                               // negateOp
            UPstream::msgType(), //S.A {13.09.22}
            cellInterpolationMap().comm() //S.A {13.09.22}
        );

        // All the donor cells will now have a valid cell centre. Construct a
        // stencil for these.

        DynamicList<label> donorCells(mesh_.nCells());
        forAll(samples, cellI)
        {
            if (samples[cellI] != greatPoint)
            {
                donorCells.append(cellI);
            }
        }


        // Get neighbours (global cell and centre) of donorCells.
        labelListList donorCellCells(mesh_.nCells());
        pointListList donorCellCentres(mesh_.nCells());
        globalCellCells
        (
            globalCells,
            mesh_,
            isValidDonor,
            donorCells,
            donorCellCells,
            donorCellCentres
        );

        // Determine the weights.
        scalarListList donorWeights(mesh_.nCells());
        forAll(donorCells, i)
        {
            label cellI = donorCells[i];
            const pointList& donorCentres = donorCellCentres[cellI];
            stencilWeights
            (
                samples[cellI],
                donorCentres,
                donorWeights[cellI]
            );
        }

        // Transfer the information back to the acceptor:
        // - donorCellCells : stencil (with first element the original donor)
        // - donorWeights : weights for donorCellCells
        cellInterpolationMap().distribute(donorCellCells);
        cellInterpolationMap().distribute(donorWeights);
        cellInterpolationMap().distribute(samples);

        // Check which acceptor has won and transfer
        forAll(interpolationCells_, i)
        {
            if (!doneAcceptor[i])
            {
                label cellI = interpolationCells_[i];
                const labelList& slots = cellStencil_[cellI];

                if (slots.size() != 1)
                {
                    FatalErrorInFunction << "Problem:" << slots
                        << abort(FatalError);
                }

                label slotI = slots[0];

                // Important: check if the stencil is actually for this cell
                if (samples[slotI] == mesh_.cellCentres()[cellI])
                {
                    cellStencil_[cellI].transfer(donorCellCells[slotI]);
                    cellInterpolationWeights_[cellI].transfer
                    (
                        donorWeights[slotI]
                    );
                    // Mark cell as being done so it does not get sent over
                    // again.
                    doneAcceptor.set(i);
                }
            }
        }
    }

    // Re-do the mapDistribute
    List<Map<label>> compactMap;
    cellInterpolationMap_.reset
    (
        new mapDistribute
        (
            globalCells,
            cellStencil_,
            compactMap
        )
    );

    Info<< "Execution time for creating the extended cellStencil = "
        << runTime.elapsedCpuTime() - timeBefore
        << " s" << endl;
}



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::cellCellStencils::foamExtendStencil::foamExtendStencil
(
    const fvMesh& mesh,
    const dictionary& dict,
    const bool doUpdate
)
:
    cellCellStencil(mesh),
    dict_(dict),
    cellTypes_(labelList(mesh.nCells(), CALCULATED)),
    interpolationCells_(0),
    cellInterpolationMap_(),
    cellStencil_(0),
    cellInterpolationWeights_(0),
    cellInterpolationWeight_
    (
        IOobject
        (
            "cellInterpolationWeight",
            mesh_.facesInstance(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        mesh_,
        dimensionedScalar(dimless, Zero),
        zeroGradientFvPatchScalarField::typeName
    ),
    om_(oversetMesh::New(mesh)),
    restartCondition_(true)
{
    // Protect local fields from interpolation
    nonInterpolatedFields_.insert("cellInterpolationWeight");
    nonInterpolatedFields_.insert("cellTypes");
    nonInterpolatedFields_.insert("maxMagWeight");

    // For convenience also suppress frequently used displacement field
    nonInterpolatedFields_.insert("cellDisplacement");
    nonInterpolatedFields_.insert("grad(cellDisplacement)");
    const word w("snGradCorr(cellDisplacement)");
    const word d("((viscosity*faceDiffusivity)*magSf)");
    nonInterpolatedFields_.insert("surfaceIntegrate(("+d+"*"+w+"))");

    // Read zoneID
    this->zoneID();

    // Read old-time cellTypes
    IOobject io
    (
        "cellTypes",
        mesh_.time().timeName(),
        mesh_,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE,
        false
    );
    if (io.typeHeaderOk<volScalarField>(true))
    {
        if (debug)
        {
            Pout<< "Reading cellTypes from time " << mesh_.time().timeName()
                << endl;
        }

        const volScalarField volCellTypes(io, mesh_);
        forAll(volCellTypes, celli)
        {
            // Round to integer
            cellTypes_[celli] = volCellTypes[celli];
        }
    }

    if (doUpdate)
    {
        update();
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::cellCellStencils::foamExtendStencil::~foamExtendStencil()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::cellCellStencils::foamExtendStencil::update()
{

////////////////////////////////////////////////////////
    // This part allows the OGA to work at specified freq.
    // We need this when no-motion is envolved when
    // the steady-state solution is purpose of the run
    // or when we need to reduce the computational cost
    // of a dynamic mesh transient analysis
    // S.A {19.06.21}

    // Re-read dictionary. small amountof time compared
    // to the OGA. Also very useful to be able
    // to modify on-the-fly.
    dictionary oversetDict
    (
        IOdictionary
        (
            IOobject
            (
                "dynamicMeshDict",
                mesh_.time().constant(),
                mesh_,
                IOobject::MUST_READ_IF_MODIFIED,
                IOobject::NO_WRITE,
                false
            )
        ).optionalSubDict(typeName + "Coeffs")
    );
    const label oversetFreq = oversetDict.get<label>("oversetFrequency");

    if (oversetFreq < 1)
    {
        FatalErrorInFunction
            << "Illegal oversetFrequency " << oversetFreq << nl
            << "The oversetFrequency setting in the dynamicMeshDict should"
            << " be >= 1." << nl
            << exit(FatalError);
    }
    //if (mesh_.time().timeIndex() > 1 && mesh_.time().timeIndex() % oversetFreq != 0
    if (!restartCondition_ && mesh_.time().timeIndex() % oversetFreq != 0)
    {
        return true;
    }

    if(restartCondition_)
    {
        restartCondition_ = false;
    }

////////////////////////////////////////////////////////

    if (debug)
    {
        Info<<"foamExtendStencil::update()"<<endl;
    }

////////////////////////////////////////////////////////

    if (debug)
    {
        Info<<"\nCreate the basic cellStencil: only the master donors\n"<<endl;
    }
    const Time& runTime = mesh_.time();
    scalar timeBefore = runTime.elapsedCpuTime();

    om_.movePoints(); //update foam-Extend overset mesh

    const globalIndex globalCells(mesh_.nCells());
    cellStencil_.setSize(mesh_.nCells());
    cellInterpolationWeights_.setSize(mesh_.nCells());
    interpolationCells_.setSize(om_.acceptorCells().size()); //triggering overset addressing too
    label totAcc = 0;

    // Loop through all overset regions
    forAll (om_.regions(), regionI)
    {
        // Get acceptors for this region
        const donorAcceptorList& curAcceptors =
            om_.regions()[regionI].acceptors();
        // Loop through acceptors of this region
        forAll (curAcceptors, aI)
        {
            // Get current donor/acceptor pair
            const donorAcceptor& da = curAcceptors[aI];
            // Get acceptor cell number
            label accID = da.acceptorCell();
            // Get the cell number of the master donor
            label donMasID = da.donorCell(); 

            // Adding the main donor as the first item in the slot
            DynamicList<label> slot; label donorProc = da.donorProcNo();
            slot.append(globalCells.toGlobal(donorProc, donMasID));

            // Get the donors for this acceptor
            cellStencil_[accID] = slot;

            // Fill the Interpolation stencil by one from each donor
            scalarList slotWeight(slot.size(), 1.0);
            cellInterpolationWeights_[accID] = slotWeight;
 

            // Fill the interpolation factors for the acceptor cells
            // Here all considered as one, so the interpolated/acceptors cells
            // take there values fully from the interpolation and nothing from
            // the old field values: newPhi = f*oldPhi + f*s
            cellInterpolationWeight_[accID] = 1.0;


            // Fill the IDs of acceptor cells
            interpolationCells_[totAcc] = accID;
            totAcc++;
        }
    }

    Info<< "Execution time for creating the basic cellStencil = "
        << runTime.elapsedCpuTime() - timeBefore
        << " s" << endl;

////////////////////////////////////////////////////////

    // Fill the cellTypes according to ESI
    // CALCULATED = 0, INTERPOLATED = 1, HOLE = 2
    // From the extend cellTypes
    // ACTIVE = 0, DONOR = 1, ACCEPTOR = 2, HOLE = -1

    cellTypes_ = CALCULATED; // reset the cellTypes_
    const labelList& holes = om_.holeCells();
    forAll (holes, cellI)
    {
        cellTypes_[holes[cellI]] = HOLE;
    }

    const labelList& acceptors = om_.acceptorCells();
    forAll (acceptors, cellI)
    {
        cellTypes_[acceptors[cellI]] = INTERPOLATED;
    }

    // Re-do the mapDistribute
    List<Map<label>> compactMap;
    cellInterpolationMap_.reset
    (
        new mapDistribute
        (
            globalCells,
            cellStencil_,
            compactMap
        )
    );

    // Extend stencil to the neighbours
    extendCellStencil(globalCells);

    // checking
    if (debug)
    {
        /*
        label totDnrsMain = 0;
        label totDnrsExt = 0;
        forAll(cellStencil_, slotID)
        {
            if(cellStencil_[slotID].size() > 0)
            {
                totDnrsMain++;
                totDnrsExt += cellStencil_[slotID].size();
            }
        }
        Info<<"Total number of main donors: "<<returnReduce(totDnrsMain, sumOp<label>())<<" and should be: "<<returnReduce(om_.donorCells().size(), sumOp<label>())<<endl;
        Info<<"Total number of donors (main and it's neighbours): "<<returnReduce(totDnrsExt, sumOp<label>())<<endl;
        */

        labelList nCells(count(3, cellTypes_));
        label nLocal = 0;
        label nMixed = 0;
        label nRemote = 0;
        forAll(interpolationCells_, i)
        {
            label celli = interpolationCells_[i];
            const labelList& slots = cellStencil_[celli];

            bool hasLocal = false;
            bool hasRemote = false;

            forAll(slots, sloti)
            {
                if (slots[sloti] >= mesh_.nCells())
                {
                    hasRemote = true;
                }
                else
                {
                    hasLocal = true;
                }
            }

            if (hasRemote)
            {
                if (!hasLocal)
                {
                    nRemote++;
                }
                else
                {
                    nMixed++;
                }
            }
            else if (hasLocal)
            {
                nLocal++;
            }
        }
        reduce(nLocal, sumOp<label>());
        reduce(nMixed, sumOp<label>());
        reduce(nRemote, sumOp<label>());

        Info<< "Overset analysis : nCells : "
            << returnReduce(cellTypes_.size(), sumOp<label>()) << nl
            << incrIndent
            << indent << "calculated   : " << nCells[CALCULATED] << nl
            << indent << "interpolated : " << nCells[INTERPOLATED]
            << " (interpolated from local:" << nLocal
            << "  mixed local/remote:" << nMixed
            << "  remote:" << nRemote << ")" << nl
            << indent << "hole (final) : " << nCells[HOLE] << nl
            << decrIndent << endl;
    }

    if (debug&2)
    {
        tmp<volScalarField> tfld
        (
            createField(mesh_, "CellTypes", cellTypes_)
        );
        tfld().write();

        tmp<volScalarField> tfld2
        (
            createField(mesh_, "zoneID", zoneID())
        );
        tfld2().write();
    }


    return true;
}


// ************************************************************************* //
