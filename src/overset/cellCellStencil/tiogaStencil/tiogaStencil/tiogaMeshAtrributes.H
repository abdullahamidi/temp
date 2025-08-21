

#ifndef tiogaMeshAtrributes_H
#define tiogaMeshAtrributes_H

namespace Foam
{

// Forward declaration of classes
class Istream;
class Ostream;

// Forward declaration of friend functions and operators
class tiogaMeshAtrributes;

inline bool operator==(const tiogaMeshAtrributes&, const tiogaMeshAtrributes&);
inline bool operator!=(const tiogaMeshAtrributes&, const tiogaMeshAtrributes&);

Istream& operator>>(Istream&, tiogaMeshAtrributes&);
Ostream& operator<<(Ostream&, const tiogaMeshAtrributes&);


struct tiogaMeshAtrributes
{   
    int btag_ = 0; // mesh tag
    int nnodes_ = 0; // number of nodes 
    double* xyz_ = NULL; // points coordinates array {size = 3*nnodes}
    int* iblank_; // {size = nnodes} {Fring/Interpolated = -1, Hole = 0, Field/Calculated = 1}
    int* iblank_cell_; // {size = number of cells} {Fring/Interpolated = -1, Hole = 0, Field/Calculated = 1}
    int nwbc_ = 0; // number of wall boundary points
    int nobc_ = 0; // number of open/overset boundary points
    int* wbcnode_ = NULL; // {size = nwbc} 
    int* obcnode_ = NULL; // {size = nobc} 
    int ntypes_ = 0; // number of cell types
    int* nv_ = NULL; // number of vertices for each cell type
    int* nc_ = NULL; // number of cells for each cell type
    int** vconn_ = NULL; // connectivity of each kind of cell
    uint64_t* cell_gid_ = NULL; // Global cell IDs
    uint64_t* node_gid_ = NULL; // Global node IDs
    int ncells_ = 0; //number of cells
    int dcount_ = 0; // number of main donors (vertices)
    int fcount_ = 0; // number of total donors (vertices)
    std::vector<int> receptorInfo_; //receptorInfo = {receptorProcess Id, receptor Index, receptor Block id, number of fractions}*dcount
    std::vector<int> inode_; //inode = indices for each receptor one group after the other
    std::vector<double> frac_; //frac  = weights for each receptor one group after the other
    
    // for getting the cell ID right after TIOGA's calls
    std::vector<int> hexCellsIDs_;
    std::vector<int> prismCellsIDs_;
    std::vector<int> pyrCellsIDs_;
    std::vector<int> tetCellsIDs_;
    
    std::vector<int> cellIdsMap_;
    DynamicList<label> vertexIdsMap_;

    void clear()
    {
        free(xyz_);
        free(iblank_);
        free(iblank_cell_);
        free(wbcnode_);
        free(obcnode_);
        free(nv_);
        free(nc_);

        for(int i=0; i<ntypes_; i++) free(vconn_[i]);
        free(vconn_);

        free(cell_gid_);
        free(node_gid_);
        
        receptorInfo_.clear(); receptorInfo_.shrink_to_fit();
        inode_.clear(); inode_.shrink_to_fit();
        frac_.clear(); frac_.shrink_to_fit();
        hexCellsIDs_.clear(); hexCellsIDs_.shrink_to_fit();
        prismCellsIDs_.clear(); prismCellsIDs_.shrink_to_fit();
        pyrCellsIDs_.clear(); pyrCellsIDs_.shrink_to_fit();
        tetCellsIDs_.clear(); tetCellsIDs_.shrink_to_fit();
        cellIdsMap_.clear(); cellIdsMap_.shrink_to_fit();
        vertexIdsMap_.clear(); vertexIdsMap_.setCapacity(0);
    }
    
             //- Default constructor
    inline tiogaMeshAtrributes()
    {        
    }
    
             //- Copy assignement constructor
    inline tiogaMeshAtrributes& operator=(const tiogaMeshAtrributes& tmesh)
    {
        if (this == &tmesh)
        {
            return *this;
        }
        
            nnodes_ = tmesh.nnodes_;
            vertexIdsMap_ = tmesh.vertexIdsMap_;      
            return *this;
    }
    
             //- assignement constructor
    inline tiogaMeshAtrributes(const tiogaMeshAtrributes& tmesh)
    {        
            nnodes_ = tmesh.nnodes_;
            vertexIdsMap_ = tmesh.vertexIdsMap_;              
    }
    
    friend bool operator==
    (
    const tiogaMeshAtrributes& a,
    const tiogaMeshAtrributes& b
    )
    {
        return a.btag_ == b.btag_;
    }

    friend bool operator!=
    (
        const tiogaMeshAtrributes& a,
        const tiogaMeshAtrributes& b
    )
    {
        return !(a == b);
    }
    
            //- Istream operator
    friend Istream& operator>>(Istream& is, tiogaMeshAtrributes& tmesh)
    {
            is  >> tmesh.nnodes_
                >> tmesh.vertexIdsMap_;
            return is;
    }
    
    friend Ostream& operator<<(Ostream& os, const tiogaMeshAtrributes& tmesh)
    {
            os  << tmesh.nnodes_ << token::SPACE
                << tmesh.vertexIdsMap_ << token::SPACE;                                
            return os;
    }       
};

} // End namespace Foam

#endif

