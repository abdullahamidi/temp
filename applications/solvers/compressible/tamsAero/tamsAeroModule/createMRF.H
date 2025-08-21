

//// S.A {Sep. 2022}
MRF_.reset(new IOMRFZoneList(mesh));
if(MRF_->active())
{
    Info << "Creating field minusMRFFlux" << nl << endl;
    minusMRFFlux_.reset
    (
        new surfaceScalarField
        (
            IOobject
            (
                "minusMRFFlux",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionedScalar("zero", dimensionSet(0, 3, -1, 0, 0), 0.0)
        )
    );
    surfaceScalarField& minusMRFFlux = minusMRFFlux_();

    // Get the { -(omega X r).Sf } for all the internal faces and for all the EXCLUDED boundary patches.
    // We need this value at all the rotating walls not to be zero (to be used in the velocity BC),
    // thus, the rotating walls shall be classified as EXCLUDED by adding them to the nonRotatingPatches
    // in the MRFProperties dict. It is a bit confusing, however, the purpose was to use this functionality
    // from the OpenFOAM MRF class to avoid code duplication.
    MRF_->makeRelative(minusMRFFlux);
}
////


//// S.A {Sep. 2022}
// Get the { omega X r } for all the rotating walls.
// To be used with rotating frame Navier-Stokes simulations
if(MRF_->active() && !inviscid_)
{
    Info << "Creating field omegaCrossR" << nl << endl;
    omegaCrossR_.reset
    (
        new surfaceVectorField
        (
            IOobject
            (
                "omegaCrossR",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionedVector("zero", dimensionSet(0, 1, -1, 0, 0), vector (0.0, 0.0, 0.0))
        )
    );
    surfaceVectorField& omegaCrossR = omegaCrossR_();
    surfaceVectorField::Boundary& omegaCrossR_b = omegaCrossR.boundaryFieldRef();

    const MRFZone& mrfZone(MRF_->getFromName("MRF1")); // only one MRF Zone for now
    //for (const auto& mrfZone: *MRF_)
    //{
        const vector Omega = mrfZone.Omega();
        const vector origin = mrfZone.origin();
        forAll(omegaCrossR_b, patchi)
        {
            if (mesh.boundary()[patchi].type() == "wall")
            {
                const vectorField& patchC = mesh.Cf().boundaryField()[patchi];
                fvsPatchVectorField& omegaCrossR_b_prop = omegaCrossR_b[patchi];
                forAll(omegaCrossR_b_prop, facei)
                {
                    omegaCrossR_b_prop[facei] = (Omega ^ (patchC[facei] - origin));
                }
            }
        }
    //}
}
////

