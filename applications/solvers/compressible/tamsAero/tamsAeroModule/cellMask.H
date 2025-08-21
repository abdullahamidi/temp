////S.A {Dec. 2023}
if(overset_)
{
    Info<<"Update the cellType list (cellMask)"<<'\n'<<endl;
    const cellCellStencil& overlap = Stencil::New(mesh);
    const labelList& cellTypes = overlap.cellTypes();

    forAll(cellTypes, i) {
        (*cellTypes_)[i] = cellTypes[i];
    }

    //label nZero = 0;
    volScalarField& cellMask = cellMask_();
    forAll(cellMask, cellI)
    {
        if (cellTypes[cellI] == cellCellStencil::HOLE || cellTypes[cellI] == cellCellStencil::INTERPOLATED)
        {
            cellMask[cellI] = 0.0;
            //nZero++;
        }
        else
        {
            cellMask[cellI] = 1.0;
        }
    }
    //reduce(nZero, sumOp<label>());
    //Info<<"nZero: "<<nZero<<endl;
}
////
