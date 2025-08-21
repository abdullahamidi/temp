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

#include "tiogaStencil.H"
#include "addToRunTimeSelectionTable.H"
#include "fvCFD.H"
#include "donorAcceptorTiogaList.H"
#include "prismMatcher.H"
#include "pointMesh.H"
#include "pointFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace cellCellStencils
{
    defineTypeNameAndDebug(tiogaStencil, 0);
}
}


void tiogaRegisterGridData(const fvMesh& msh, tiogaMeshAtrributes* tmsh)
{
    const cellZoneMesh& cz = msh.cellZones();
    if(cz.size() <= 1)
    {
        FatalErrorInFunction
            << "There is only one (or no) cellzone in the mesh"
            << "Make sure that you have a cellzone for each of your meshes in the overset grid system"
            << exit(FatalError);
    }

    const cellList& cells = msh.cells();
    
    const cellShapeList& shapes = msh.cellShapes();
    const cellModel& hexModel = cellModel::ref(cellModel::HEX);
    const cellModel& prismModel = cellModel::ref(cellModel::PRISM);
    const cellModel& pyrModel = cellModel::ref(cellModel::PYR);
    const cellModel& tetModel = cellModel::ref(cellModel::TET);

    const faceList& faces = msh.faces();
    const pointList& points = msh.points();
    const polyBoundaryMesh& boundaryMesh = msh.boundaryMesh();
    bitSet isZonePoint(0);
    bitSet isBoundaryPoint(0);
    DynamicList<point> pointList(0);
    std::vector<double> xyz;
    std::vector<int> iblank;
    std::vector<int> iblank_cell;
    std::vector<int> wbcnode;
    std::vector<int> obcnode;
    std::vector<int> nv;
    std::vector<int> nc;
    std::vector<std::vector<int>> vconn;
    std::vector<uint64_t> cell_gid;
    std::vector<uint64_t> node_gid;
    int nwbcFaces;
    int nobcFaces;
    int nHex;
    int nPrism;
    int nPyr;
    int nTet;
    prismMatcher prism;
    std::vector<int> vconnHex; // connectivity of Hex cell type
    std::vector<int> vconnPrism; // connectivity of Prism cell type
    std::vector<int> vconnPyr; // connectivity of Pyr cell type
    std::vector<int> vconnTet; // connectivity of Tet cell type

    int totNumOfVert = 0;

    for (const cellZone& cZone : cz)
    {
        tiogaMeshAtrributes& tma = tmsh[cZone.index()];
        tma.btag_ = cZone.index()+1; // TIOGA expects the btag to start from 1
        isZonePoint.reset();
        isBoundaryPoint.reset();
        pointList.clear();
        xyz.resize(points.size()*3);
        iblank.clear();
        iblank_cell.clear();
        wbcnode.clear();
        obcnode.clear();
        nv.clear();
        nc.clear();
        vconn.clear();
        cell_gid.clear();
        node_gid.clear();
        nwbcFaces = 0;
        nobcFaces = 0;
        nHex = 0;
        nTet = 0;
        nPrism = 0;
        nPyr = 0;
        vconnHex.clear();
        vconnPrism.clear();
        vconnPyr.clear();
        vconnTet.clear();
        for (const label celli : cZone)
        {
            for (const label facei : cells[celli]) // loop over all the faces
            { 
                const face& f = faces[facei];
                for (const label verti : f)
                {
                    if (isZonePoint.set(verti))
                    {
                        pointList.append(points[verti]);
                        xyz[3*(verti - totNumOfVert)  ] = points[verti].x();
                        xyz[3*(verti - totNumOfVert)+1] = points[verti].y();
                        xyz[3*(verti - totNumOfVert)+2] = points[verti].z();
                        iblank.push_back(-(cZone.index()+10));
                        node_gid.push_back(verti);
                    }
                }
                if(!msh.isInternalFace(facei))
                {
                     const label& patchi = boundaryMesh.whichPatch(facei);  // patch index
                     if (!boundaryMesh[patchi].coupled() && boundaryMesh[patchi].type() == "wall") //exclude the intra-processors patches and keep only the walls
                     {
                         nwbcFaces++;
                         for (const label verti : f)
                         {
                             if (isBoundaryPoint.set(verti))
                             {
                                 wbcnode.push_back(verti - totNumOfVert);

                                 tma.nwbc_++;
                             }
                         }
                     }
                     else if (!boundaryMesh[patchi].coupled() && boundaryMesh[patchi].type() == "overset") //exclude the intra-processors patches and keep only the overset
                     {
                         nobcFaces++;
                         for (const label verti : f)
                         {
                             if (isBoundaryPoint.set(verti))
                             {
                                 obcnode.push_back(verti - totNumOfVert);
                                 tma.nobc_++;
                             }
                         }
                     }
                }            
            }
            // Cells shape recognition
            const labelList& cPoints = shapes[celli];
            cell_gid.push_back(celli+1);
            iblank_cell.push_back(-(cZone.index()+10));
            if (shapes[celli].model() == hexModel)
            {
                nHex++; tma.hexCellsIDs_.push_back(celli);
                for (int i=0; i<cPoints.size(); i++)
                {
                    vconnHex.push_back(cPoints[i]+1 - totNumOfVert);
                }
            }
            else if (shapes[celli].model() == prismModel)
            {
                nPrism++; tma.prismCellsIDs_.push_back(celli);
                for (int i=0; i<cPoints.size(); i++)
                {
                    vconnPrism.push_back(cPoints[i]+1 - totNumOfVert);
                }
            }
            else if (shapes[celli].model() == pyrModel)
            {
                nPyr++; tma.pyrCellsIDs_.push_back(celli);
                for (int i=0; i<cPoints.size(); i++)
                {
                    vconnPyr.push_back(cPoints[i]+1 - totNumOfVert);
                }
            }
            else if (shapes[celli].model() == tetModel)
            {
                nTet++; tma.tetCellsIDs_.push_back(celli);
                for (int i=0; i<cPoints.size(); i++)
                {
                    vconnTet.push_back(cPoints[i]+1 - totNumOfVert);
                }
            }
            else
            {
                FatalErrorInFunction
                    << "The TIOGA library supports only four cell types: hexahedra, tetrahedra, prisms and pyramids"
                    << exit(FatalError);
            }
        }
        tma.nnodes_ = pointList.size(); // all the points including the ones on the boundary
        totNumOfVert += pointList.size(); 
        tma.ncells_ = cell_gid.size();
        tma.cellIdsMap_.resize(tma.ncells_);
        
        xyz.resize(3*tma.nnodes_); xyz.shrink_to_fit();
        if(xyz.size() != static_cast<unsigned int>(3*tma.nnodes_))
        {
            FatalErrorInFunction
                << "You are not reading the mesh nodes correctly!"
                << exit(FatalError);
        }
        else
        {   
            tma.xyz_ = static_cast<double *>(malloc(sizeof(double)*3*tma.nnodes_ ));
            for(int i=0; i<(3*tma.nnodes_); i++) tma.xyz_[i] = xyz[i];
            xyz.clear();
        }
        tma.iblank_ = static_cast<int *>(malloc(sizeof(int)*( tma.nnodes_ )));
        for(int i=0; i<(tma.nnodes_); i++) tma.iblank_[i] = iblank[i];

        tma.iblank_cell_ = static_cast<int *>(malloc(sizeof(int)*( cell_gid.size() )));
        for(unsigned int i=0; i<(cell_gid.size()); i++) tma.iblank_cell_[i] = iblank_cell[i];

        tma.wbcnode_ = static_cast<int *>(malloc(sizeof(int)*( wbcnode.size() )));
        for(unsigned int i=0; i<(wbcnode.size()); i++) tma.wbcnode_[i] = wbcnode[i]+1;

        tma.obcnode_ = static_cast<int *>(malloc(sizeof(int)*( obcnode.size() )));
        for(unsigned int i=0; i<(obcnode.size()); i++) tma.obcnode_[i]  = obcnode[i]+1;

        if(nHex>0)
        {
            nv.push_back(8);
            nc.push_back(nHex);
            vconn.push_back(vconnHex);
        }
        if(nPrism>0)
        {
            nv.push_back(6);
            nc.push_back(nPrism);
            vconn.push_back(vconnPrism);
        }
        if(nPyr>0)
        {
            nv.push_back(5);
            nc.push_back(nPyr);
            vconn.push_back(vconnPyr);
        }
        if(nTet>0)
        {
            nv.push_back(4);
            nc.push_back(nTet);
            vconn.push_back(vconnTet);
        }
        tma.ntypes_ = nv.size();

        tma.nv_ = static_cast<int *>(malloc(sizeof(int)*( nv.size() )));
        for(unsigned int i=0; i<(nv.size()); i++) tma.nv_[i] = nv[i];

        tma.nc_ = static_cast<int *>(malloc(sizeof(int)*( nc.size() )));
        for(unsigned int i=0; i<(nc.size()); i++) tma.nc_[i] = nc[i];

        tma.vconn_ = static_cast<int **>(malloc(sizeof(int*)*( tma.ntypes_ )));
        for(int i=0; i<tma.ntypes_; i++) tma.vconn_[i] = static_cast<int *>(malloc(sizeof(int)*( vconn[i].size() )));

        for(int i=0; i<tma.ntypes_; i++)
        {
            for(unsigned int j=0; j<vconn[i].size(); j++) tma.vconn_[i][j] = vconn[i][j];
        }

        tma.cell_gid_ = static_cast<uint64_t *>(malloc(sizeof(uint64_t)*( cell_gid.size() )));
        for(uint64_t i=0; i<cell_gid.size(); i++) tma.cell_gid_[i] = cell_gid[i];

        tma.node_gid_ = static_cast<uint64_t *>(malloc(sizeof(uint64_t)*( node_gid.size() )));
        for(uint64_t i=0; i<node_gid.size(); i++) tma.node_gid_[i] = node_gid[i];

// some ouputs for debugging
/*
        
        {
        Info<< "Overall number of nodes for zone :"<<cZone.name()<<" is: " <<returnReduce(pointList.size(),sumOp<label>())<<endl;

        reduce(nHex,sumOp<label>()); reduce(nTet,sumOp<label>()); reduce(nPrism,sumOp<label>()); reduce(nPyr,sumOp<label>());
        Info<< "Overall number of cells of each type for zone :"<<cZone.name()<< nl
            << "    hexahedra:     " << nHex <<nl
            << "    prisms:        " << nPrism << nl
            << "    pyramids:      " << nPyr << nl
            << "    tetrahedra:    " << nTet << nl
            << "    ntypes:        " << tma.ntypes_ << nl
            << "    total:         " << nHex + nPrism + nPyr + nTet << endl;

        Info<< "Overall number of points included each type of cells for zone :"<<cZone.name()<< nl
            << "    hexahedra:     " << returnReduce(vconnHex.size(),sumOp<label>()) <<nl
            << "    prisms:        " << returnReduce(vconnPrism.size(),sumOp<label>()) << nl
            << "    pyramids:      " << returnReduce(vconnPyr.size(),sumOp<label>()) << nl
            << "    tetrahedra:    " << returnReduce(vconnTet.size(),sumOp<label>()) << nl
            << "    total:         " << returnReduce(vconnHex.size(),sumOp<label>()) + returnReduce(vconnPrism.size(),sumOp<label>()) +
                                        returnReduce(vconnPyr.size(),sumOp<label>()) + returnReduce(vconnTet.size(),sumOp<label>()) << endl;
        }
*/        
    }

}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::cellCellStencils::tiogaStencil::tiogaStencil
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
    tg_(new TIOGA::tioga[1]), //- tioga
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

    //
    if (debug)
    {
        Info<<"TIOGA Initialization"<<endl;
    }
        // Valgrind shows memory leakage caused by the following line. 
        // Same has been reported on stackoverflow.
        // Even for a very simple mpi program (I tested my self), Valgrind showed memory leakage
        // Usually, the Valgrind leakages/warnings related to mpi are ignored/suppressed
    UPstream::initNull(); 
    tg_->setCommunicator(MPI_COMM_WORLD,Pstream::myProcNo(),Pstream::nProcs());
    // 

    if (doUpdate)
    {
        update();
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::cellCellStencils::tiogaStencil::~tiogaStencil()
{
    //- tioga
    delete [] tg_;
    // 
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::cellCellStencils::tiogaStencil::update()
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
        Info<<"tiogaStencil::update()"<<endl;
    }

////////////////////////////////////////////////////////

    const Time& runTime = mesh_.time();
    scalar timeStart{}, timeTiogaProcesses{}, timeBeforeBuilding{}, timeBeforeFilling{};
    timeStart = runTime.elapsedCpuTime();
    //- tioga
    if (debug)
    {
        Info<<"TIOGA reading mesh"<<endl;
    }
    
    tiogaMeshAtrributes *tmesh = static_cast<tiogaMeshAtrributes *>(calloc(mesh_.cellZones().size(), sizeof(struct tiogaMeshAtrributes)));
    
    tiogaRegisterGridData(mesh_, tmesh); //Preparing the mesh data for TIOGA

    if (debug)
    {
        Info<< "Required time for preparing the mesh data for TIOGA = "
        << mesh_.time().elapsedCpuTime() - timeStart << " s\n" << endl;
        timeTiogaProcesses = runTime.elapsedCpuTime();
    }    

    List<label> zeroNodesRegMap(mesh_.cellZones().size());
    label regZones = 0;
    for (int i=0; i<mesh_.cellZones().size(); i++)
    {
        // the below if statement will allow any type of parallel decomposition to work with TIOGA, 
        // even if it produces a decoposition where a processor doesn't contain any node.
        if(tmesh[i].nnodes_ > 0)
        {
            tg_->registerGridData(tmesh[i].btag_, tmesh[i].nnodes_, tmesh[i].xyz_, tmesh[i].iblank_, tmesh[i].nwbc_,
                                tmesh[i].nobc_, tmesh[i].wbcnode_, tmesh[i].obcnode_, tmesh[i].ntypes_, tmesh[i].nv_,
                                tmesh[i].nc_, tmesh[i].vconn_);//, tmesh[i].cell_gid_, tmesh[i].node_gid_);

            // Indicating that the element/cell IBLANK information are required and to be returned by TIOGA
            //tg_->set_cell_iblank(tmesh[i].btag_, tmesh[i].iblank_cell_);
            zeroNodesRegMap[regZones] = i; regZones++;
        }
    }


    if (debug)
    {
        Info<<"TIOGA preprocess_grids "<<endl;
    }
    tg_->profile();

    if (debug)
    {
        Info<<"TIOGA performconnectivity "<<endl;
    }    
    tg_->performConnectivity();
    //tg_->outputHoleMap();

    const bool reduceFringes = oversetDict.getOrDefault<bool>("reduceFringes", false);
    if (reduceFringes)
    {
        tg_->reduce_fringes();
        if (debug)
        {
            Info<<"TIOGA reduce_fringes "<<endl;
        }
    }  
    
    if (debug)
    {
        Info<<"Collect the cells classification from TIOGA"<<endl;
    }
    
    // Read the iblank outputs from TIOGA

    pointScalarField iblank_points
    (
        IOobject
        (
            "iblank_points",
            runTime.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        pointMesh::New(mesh_),
        dimensionedScalar("zero", dimless, 0.0)
    );
    
    int totPoints = 0;
    for (int i=0; i<mesh_.cellZones().size(); i++)
    {
        label nn = tmesh[i].nnodes_;
        tmesh[i].vertexIdsMap_.resize(nn);
        for (int j=0; j<nn; j++)
        {
            tmesh[i].vertexIdsMap_[j] = j+totPoints; //tiogavertexID to OFvertexID map
            iblank_points[j+totPoints] = tmesh[i].iblank_[j];
        }
        totPoints += nn;
    }  


    int totCells = 0;

    for (int i = 0; i < mesh_.cellZones().size(); i++)
    {
        int cellid = 0;
        const int ncells = tmesh[i].ncells_;
        const int* nv = tmesh[i].nv_;
        const int* nc = tmesh[i].nc_;

        for (int j = 0; j < tmesh[i].ntypes_; j++)
        {
            int nvert = nv[j];

            for (int l = 0; l < nc[j]; l++)
            {
                int ofCellID{};

                switch (nvert)
                {
                    case 8:
                        ofCellID = tmesh[i].hexCellsIDs_[l];
                        break;
                    case 6:
                        ofCellID = tmesh[i].prismCellsIDs_[l];
                        break;
                    case 5:
                        ofCellID = tmesh[i].pyrCellsIDs_[l];
                        break;
                    case 4:
                        ofCellID = tmesh[i].tetCellsIDs_[l];
                        break;
                    default:
                        FatalErrorInFunction
                            << "The TIOGA library supports only four cell types: hexahedra, tetrahedra, prisms, and pyramids."
                            << exit(FatalError);
                }

                tmesh[i].cellIdsMap_[cellid] = ofCellID;
                cellid++;
            }
        }

        totCells += ncells;
    }
    
    cellTypes_ = CALCULATED; // reset the cellTypes_   
    
    const labelListList& pointCells = mesh_.pointCells();
    for (int i=0; i<pointCells.size(); i++)
    {
        if(iblank_points[i] == 0)   //hole node
        {
            const labelList& pCells = pointCells[i];
            for(int q=0; q<pCells.size(); q++)
            {
                cellTypes_[ pCells[q] ] = 2;                                                                   
            }        
        }  
    }
    

    if (debug)
    {
        Info<<"getDonorCount from TIOGA"<<endl;
    }
    for (int i=0; i<mesh_.cellZones().size(); i++)
    {    
        tg_->getDonorCount(tmesh[i].btag_, &tmesh[i].dcount_, &tmesh[i].fcount_);
        tmesh[i].receptorInfo_.resize(4 * tmesh[i].dcount_);
        tmesh[i].inode_.resize(tmesh[i].fcount_);
        tmesh[i].frac_.resize(tmesh[i].fcount_);
        if (debug&2)
        {        
            Info<<"btag: "<< tmesh[i].btag_
                <<" Number of main donors (vertices) dcount: "<< returnReduce(tmesh[i].dcount_,sumOp<label>())
                <<", Number of total donors (vertices) fcount: "<< returnReduce(tmesh[i].fcount_,sumOp<label>())<<endl;
        }
    }          

    if (debug)
    {
        Info<<"getDonorInfo from TIOGA"<<endl;
    }
    
    for (int i=0; i<mesh_.cellZones().size(); i++)
    {
        tg_->getDonorInfo(tmesh[i].btag_, tmesh[i].receptorInfo_.data(), tmesh[i].inode_.data(), tmesh[i].frac_.data(), &tmesh[i].dcount_);           
    }
       
// Print Acceptor/Donors data
/*
if (Pstream::myProcNo() == 0)
{  
    for (int h=0; h<mesh_.cellZones().size(); h++)
    {
     if(tmesh[h].nnodes_ > 0)
     {
        int m = 0;
        std::ofstream donorFile("region" + std::to_string(tmesh[h].btag_-1) + ".txt");
        if (donorFile.is_open())
        {
            donorFile <<"receptorInfo[0] "<<" receptorInfo[1] "<<" receptorInfo[2] "<<" nweights "<<std::endl;
            for(int i=0; i<tmesh[h].dcount_; i++) // loop over donors at this region and processor
            {
                donorFile <<tmesh[h].receptorInfo_[4*i]<<" "<<tmesh[h].receptorInfo_[4*i+1]<<" "
                          <<tmesh[h].receptorInfo_[4*i+2]<<" "<<tmesh[h].receptorInfo_[4*i+3]<<std::endl;
                          
                int k = tmesh[h].receptorInfo_[4*i+3] + 1;
                donorFile <<"For acceptor/receptor vertex: "<<tmesh[h].receptorInfo_[4*i+1]
                <<" we have "<<k-1<<" donor vertices which belong to donor cell "<<tmesh[h].inode_[m+k-1]<<" (tioga cellID) "
                <<" or "<<tmesh[h].cellIdsMap_[ tmesh[h].inode_[m+k-1] ]<<" (OpenFOAM cellID) "<<std::endl;
                for(int j=0; j<k-1; j++)
                {
                    donorFile <<"donor vertex "<<j<<" ID is "<<tmesh[h].inode_[j+m]<<std::endl;
                }
                m += k;          
            }

            donorFile.close();
        }
     }
    }
}    
*/
    if (debug)
    {
        Info<< "Required time for all tioga processes = "
            << mesh_.time().elapsedCpuTime() - timeTiogaProcesses << " s\n" << endl;
        timeBeforeBuilding = mesh_.time().elapsedCpuTime(); 
    }
    
    if (debug)
    {
        Info<<"Building stencil data"<<"\n"<<endl;
    }    
       
    donorAcceptorTiogaDynamicList daList;
    donorAcceptorTiogaDynamicList daListRemoteAccept;
    const pointList& CC = mesh_.C();  
    const pointList& points = mesh_.points();
    const labelListList& cellCells = mesh_.cellCells();
    typedef DynamicList<label, 10> DynamicLabelList;
    typedef DynamicList<point, 10> DynamicPointList;

    const globalIndex globalCells(mesh_.nCells());
    
        // Count how many to send
    labelList nSend(Pstream::nProcs(), Zero);
              

    for (int h=0; h<mesh_.cellZones().size(); h++)
    {
        int nnodes_ = tmesh[h].nnodes_;
        if(nnodes_ > 0)
        {
            int m = 0; 
            for(int i=0; i<tmesh[h].dcount_; i++) // loop over donors at this region and this processor
            {
                donorAcceptorTioga da;                                  
                int k = tmesh[h].receptorInfo_[4*i+3] + 1;
                int donorCell_ = tmesh[h].cellIdsMap_[ tmesh[h].inode_[m+k-1] ]; // OpenFOAM local cellID for the main donor cell
                da.donorCell() = donorCell_;
                int donorProcNo_ = Pstream::myProcNo();
                da.donorProcNo() = donorProcNo_;
                da.donorPoint() = CC[donorCell_];
                labelList extendedDonorCells_ = cellCells[donorCell_];
                da.extendedDonorCells() = extendedDonorCells_;  // OpenFOAM local cellIDs for the extended donor cells
                DynamicLabelList& extendedDonorCellsMutable_ = da.extendedDonorCells();
                for(int q=0; q<da.extendedDonorCells().size(); q++)
                {
                    da.extendedDonorPoints().append( CC[extendedDonorCells_[q]] );
                    
                    // change the local IDs to global IDs for the OpenFOAM cellStencil
                    extendedDonorCellsMutable_[q] = globalCells.toGlobal(donorProcNo_, extendedDonorCells_[q]);
                }
                da.donorRegionID() = h;
                DynamicLabelList& mainDonorCellVertices = da.mainDonorCellVertices();
                for(int j=0; j<k-1; j++)
                {
                    label cellVertices = tmesh[h].vertexIdsMap_[tmesh[h].inode_[j+m]];
                    mainDonorCellVertices.append(cellVertices);  // OpenFOAM local pointID for the main donor vertex/node IDs                    
                    DynamicPointList& mainDonorCellPoints = da.mainDonorCellPoints();
                    mainDonorCellPoints.append( points[cellVertices] ); // coordinates for the main donor vertices/nodes                    
                }
                m += k; 
               
                int acceptorProcID = tmesh[h].receptorInfo_[4*i];
                da.acceptorProcNo() = acceptorProcID;
                int acceptorVertexID = tmesh[h].receptorInfo_[4*i+1];
                int acceptorRegionID = tmesh[h].receptorInfo_[4*i+2];
                if(acceptorProcID == Pstream::myProcNo())
                {
                    acceptorRegionID = zeroNodesRegMap[acceptorRegionID]; //zeroNodesRegMapLists[acceptorProcID][acceptorRegionID];
                    da.acceptorRegionID() = acceptorRegionID;
                    da.acceptorVertexID() =  tmesh[acceptorRegionID].vertexIdsMap_[acceptorVertexID];
                    const labelList& pCells = pointCells[da.acceptorVertexID()];
                    da.acceptorCandidateCells() = pCells.size();
                    for(int q=0; q<da.acceptorCandidateCells(); q++)
                    {
                        da.acceptorCandidateCellsIDs().append( pCells[q] );                                             
                        da.acceptorCandidateCellsPoints().append( CC[pCells[q]] );                       
                    }
                    daList.append(std::move(da));
                }               
                else
                {
                    da.setAcceptorRemote(true); 
                    da.acceptorVertexID() = acceptorVertexID;
                    da.acceptorRegionID() = acceptorRegionID;
                    // Count how many to send
                    nSend[acceptorProcID]++;
                    daListRemoteAccept.append(std::move(da));
                }
                //scalar d;
                //da.acceptorCell() = da.bestAcceptorCandidate(d);
                //checkDistances.append(d);         
                //daList.append(std::move(da));
            }
        }
    }
 
    //Info<<"gSum(nSend): "<<gSum(nSend)<<endl;

    // Collect items to be sent
    labelListList sendMap(Pstream::nProcs());
    forAll(sendMap, procI)
    {
        sendMap[procI].setSize(nSend[procI]);
    }
    nSend = 0;
    forAll(daListRemoteAccept, i)
    {
        const label procI = daListRemoteAccept[i].acceptorProcNo();
        sendMap[procI][nSend[procI]++] = i;
    }
    
    // Sync how many to send
    labelList nRecv;
    Pstream::exchangeSizes(sendMap, nRecv);
    
    // Collect items to be received
    labelListList recvMap(Pstream::nProcs());
    forAll(recvMap, procI)
    {
        recvMap[procI].setSize(nRecv[procI]);
    }
    
    label constructSize = 0;
    // Construct with my own elements first
    forAll(recvMap[Pstream::myProcNo()], i)
    {
        recvMap[Pstream::myProcNo()][i] = constructSize++;
    }
    // Construct from other processors
    forAll(recvMap, procI)
    {
        if (procI != Pstream::myProcNo())
        {
            forAll(recvMap[procI], i)
            {
                recvMap[procI][i] = constructSize++;
            }
        }
    }
  
    // Construct distribute map (destructively)
    mapDistribute map(constructSize, std::move(sendMap), std::move(recvMap));

    // Distribute daListRemoteAccept
    map.distribute(daListRemoteAccept);    
    //Info<<"daListRemoteAccept.size(): "<<returnReduce(daListRemoteAccept.size(),sumOp<label>())<<endl;

    // Fill the remote acceptor data
    for (int i=0; i<daListRemoteAccept.size(); i++)
    {
        donorAcceptorTioga da = daListRemoteAccept[i];
     
        int acceptorRegionID = zeroNodesRegMap[da.acceptorRegionID()];
        da.acceptorRegionID() = acceptorRegionID;
        da.acceptorVertexID() =  tmesh[acceptorRegionID].vertexIdsMap_[da.acceptorVertexID()/4]; // after mapping the indices are multiplied by 4 !! 
        const labelList& pCells = pointCells[da.acceptorVertexID()]; 
        da.acceptorCandidateCells() = pCells.size();
        for(int q=0; q<da.acceptorCandidateCells(); q++)
        {
            da.acceptorCandidateCellsIDs().append( pCells[q] );                                             
            da.acceptorCandidateCellsPoints().append( CC[pCells[q]] );                       
        }
        daList.append(std::move(da));        
    }      

    if (debug)
    {
        Info<< "Building stencil: total time = "
            << mesh_.time().elapsedCpuTime() - timeBeforeBuilding << " s\n" << endl;      
        timeBeforeFilling = mesh_.time().elapsedCpuTime();
    }


    cellStencil_.setSize(mesh_.nCells());
    cellInterpolationWeights_.setSize(mesh_.nCells());
    bitSet isLocalAccept(0);
    bitSet isRemoteAccept(0);

    for (int i=0; i<daList.size(); i++)
    {
        donorAcceptorTioga da = daList[i];
        
        DynamicLabelList acceptorCandidateCellsIDs = da.acceptorCandidateCellsIDs();
        DynamicPointList acceptorCandidateCellsPoints = da.acceptorCandidateCellsPoints();
        DynamicPointList extendedDonorPoints = da.extendedDonorPoints();
        extendedDonorPoints.append(da.donorPoint()); // add the main donor
        DynamicLabelList extendedDonorCells = da.extendedDonorCells();
        extendedDonorCells.append(globalCells.toGlobal(da.donorProcNo(), da.donorCell()));  // add the main donor
        
        for (int j=0; j<da.acceptorCandidateCells(); j++)
        {
            int acceptCellID = acceptorCandidateCellsIDs[j];
            if (cellTypes_[acceptCellID] != 2 && cellTypes_[acceptCellID] != 1)  //exclude hole and calculated cells
            {
                cellTypes_[acceptCellID] = 1;
                cellStencil_[acceptCellID] = extendedDonorCells;
                stencilWeights(acceptorCandidateCellsPoints[j], extendedDonorPoints, cellInterpolationWeights_[acceptCellID]);
                    
                if(da.acceptorRemote())
                {
                    isRemoteAccept.set(acceptCellID);
                }
                else
                {
                    isLocalAccept.set(acceptCellID);             
                }         
            }
        }
    }

    if (debug)
    {    
        Info<< "Required time for filling the stencil data = "
            << mesh_.time().elapsedCpuTime() - timeBeforeFilling << " s\n" << endl;
    }    
   
    if (debug)
    { 
        Info<<"Total local acceptors: "<< returnReduce(isLocalAccept.count(),sumOp<label>())<<endl;
        Info<<"Total remote acceptors: "<< returnReduce(isRemoteAccept.count(),sumOp<label>())<<"\n"<<endl;
    }    
      

    interpolationCells_.setSize(mesh_.nCells());
    label totAcc = 0;
    label totCalc = 0;
    label totHole = 0;
    cellInterpolationWeight_.setSize(mesh_.nCells());
    for (int i=0; i<mesh_.nCells(); i++)
    {
        if(cellTypes_[i] == 1)
        {
            interpolationCells_[totAcc] = i;
            totAcc++;
            cellInterpolationWeight_[i] = 1.0;
        }
        else if(cellTypes_[i] == 0)
        {
            totCalc++;
        }
        else if(cellTypes_[i] == 2)
        {
            totHole++;
        }
    }
    
    interpolationCells_.resize(totAcc);

    if (debug)
    {
        Info<< " DCI Summary:"<< nl
            << " Number of calculated cells: " << returnReduce(totCalc,sumOp<label>()) <<nl
            << " Number of interpolated cells: " << returnReduce(totAcc,sumOp<label>()) << nl
            << " Number of hole cells: " << returnReduce(totHole,sumOp<label>()) << nl
            << endl;
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
             
    //- tioga    
    // 
    if (debug)
    {
        Info<<"TIOGA release the local allocated resources"<<endl;
    }

    // release the local allocated resources
    for (int i=0; i<mesh_.cellZones().size(); i++) tmesh[i].clear();

    free(tmesh);
    
    Info<< "Total TIOGA stencil time = "
        << mesh_.time().elapsedCpuTime() - timeStart << " s\n" << endl;    

    return true;
}


// ************************************************************************* //
