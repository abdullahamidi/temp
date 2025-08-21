// Create fields on object registry if necessary
word timeName = mesh_.time().timeName(mesh_.time().startTime().value());
if (!mesh_.foundObject<volVectorField>("F"))
{
    (
        new volVectorField
        (
            IOobject
            (
                "F",
                timeName,
                mesh_,
                IOobject::READ_IF_PRESENT,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedVector("zero", dimForce/dimVol, vector::zero)
        )
    )->store();
}

if (!mesh_.foundObject<volVectorField>("camber"))
{
    (
        new volVectorField
        (
            IOobject
            (
                "camber",
                timeName,
                mesh_,
                IOobject::READ_IF_PRESENT,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedVector("zero", dimless, vector::zero)
        )
    )->store();
}

if (!mesh_.foundObject<volScalarField>("devAngle"))
{
    (
        new volScalarField
        (
            IOobject
            (
                "devAngle",
                timeName,
                mesh_,
                IOobject::READ_IF_PRESENT,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedScalar("zero", dimless, 0.0)
        )
    )->store();
}

if (!mesh_.foundObject<volScalarField>("blockage"))
{
    (
        new volScalarField
        (
            IOobject
            (
                "blockage",
                timeName,
                mesh_,
                IOobject::READ_IF_PRESENT,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedScalar("zero", dimless, 0.0)
        )
    )->store();
}

if (!mesh_.foundObject<volVectorField>("blockageGradient"))
{
    (
        new volVectorField
        (
            IOobject
            (
                "blockageGradient",
                timeName,
                mesh_,
                IOobject::READ_IF_PRESENT,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedVector("zero", dimless/dimLength, vector::zero)
        )
    )->store();
}


