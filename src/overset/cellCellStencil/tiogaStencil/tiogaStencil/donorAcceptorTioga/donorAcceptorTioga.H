/*---------------------------------------------------------------------------*\
Class
    Foam::donorAcceptorTioga

Description
    Class holds donor and acceptor data.  Used for searching and communications

SourceFiles

Author
    Saleh Abuhanieh, Modified version from the foam-extend one.

\*---------------------------------------------------------------------------*/

#ifndef donorAcceptorTioga_H
#define donorAcceptorTioga_H

#include "label.H"
#include "point.H"
#include "DynamicList.H"
#include "boolList.H"
#include "vectorField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class Istream;
class Ostream;

// Forward declaration of friend functions and operators
class donorAcceptorTioga;

inline bool operator==(const donorAcceptorTioga&, const donorAcceptorTioga&);
inline bool operator!=(const donorAcceptorTioga&, const donorAcceptorTioga&);

Istream& operator>>(Istream&, donorAcceptorTioga&);
Ostream& operator<<(Ostream&, const donorAcceptorTioga&);


/*---------------------------------------------------------------------------*\
                        Class donorAcceptorTioga Declaration
\*---------------------------------------------------------------------------*/

class donorAcceptorTioga
{
public:

    // Public data types
    typedef DynamicList<label, 10> DynamicLabelList;
    typedef DynamicList<point, 10> DynamicPointList;


private:

    // Private data

        // Acceptor side

            //- Acceptor cell number
            label acceptorCell_;

            //- Acceptor processor number
            label acceptorProcNo_;

            //- Acceptor location (cell center)
            point acceptorPoint_;
            
            //- Acceptor region ID            
            label acceptorRegionID_;
            
            //- Acceptor vertex ID            
            label acceptorVertexID_;          
            
            //- Acceptor remote?            
            bool acceptorRemote_;

            //- Acceptor vertex ID is final            
            bool acceptorVertexIDisFinal_; 

            //- Acceptor Candidate cells                
            label acceptorCandidateCells_;
            
            //- Acceptor Candidate cells IDs
            DynamicLabelList acceptorCandidateCellsIDs_;

            //- Acceptor Candidate cells center coordinates
            DynamicPointList acceptorCandidateCellsPoints_;         


        // Donor side

            //- Donor cell number (main donor cell ID)
            label donorCell_;

            //- Donor processor number
            label donorProcNo_;

            //- Donor location (cell center)
            point donorPoint_;
            
            //- Main donor cell vertex IDs
            DynamicLabelList mainDonorCellVertices_;
            
            //- Main donor cell vertex points (coordinates)
            DynamicPointList mainDonorCellPoints_;

            //- Donor region ID            
            label donorRegionID_;

            //- Extended donor cell numbers
            DynamicLabelList extendedDonorCells_;

            //- Extended donor points, used to construct weights in different
            //  overset interpolation schemes (cell center)
            DynamicPointList extendedDonorPoints_;

            // Note that the processor number is the same for all extended
            // donors

            //- Whether donor (or acceptor) is within bounding box of acceptor
            //  (or donor)
            bool withinBB_;


public:

    // Constructors

        //- Construct null
        inline donorAcceptorTioga()
        :
            acceptorCell_(-1),
            acceptorProcNo_(-1),
            acceptorPoint_(vector::zero),
            acceptorRegionID_(-1),
            acceptorVertexID_(-1), 
            acceptorRemote_(false),
            acceptorVertexIDisFinal_(true), 
            acceptorCandidateCells_(-1),
            acceptorCandidateCellsIDs_(),
            acceptorCandidateCellsPoints_(), 
            donorCell_(-1),
            donorProcNo_(-1),
            donorPoint_(vector::zero),
            mainDonorCellVertices_(),
            mainDonorCellPoints_(),
            donorRegionID_(-1),
            extendedDonorCells_(),
            extendedDonorPoints_(),
            withinBB_(false)
        {}

        //- Construct from acceptor data
        inline donorAcceptorTioga
        (
            const label acceptorCell,
            const label acceptorProcNo,
            const point& acceptorPoint
        )
        :
            acceptorCell_(acceptorCell),
            acceptorProcNo_(acceptorProcNo),
            acceptorPoint_(acceptorPoint),
            acceptorRegionID_(-1),
            acceptorVertexID_(-1),
            acceptorRemote_(false),
            acceptorVertexIDisFinal_(true),
            acceptorCandidateCells_(-1),
            acceptorCandidateCellsIDs_(),
            acceptorCandidateCellsPoints_(),             
            donorCell_(-1),
            donorProcNo_(-1),
            donorPoint_(vector::zero),
            mainDonorCellVertices_(),
            mainDonorCellPoints_(),
            donorRegionID_(-1),
            extendedDonorCells_(),
            extendedDonorPoints_(),
            withinBB_(false)
        {}

        //- Copy constructor
        inline donorAcceptorTioga(const donorAcceptorTioga& da)
        {
            if (this == &da)
            {
                return;
            }
             
            acceptorCell_ = da.acceptorCell_;
            acceptorProcNo_ = da.acceptorProcNo_;
            acceptorPoint_ = da.acceptorPoint_;
            acceptorRegionID_ = da.acceptorRegionID_;
            acceptorVertexID_ = da.acceptorVertexID_;
            acceptorRemote_ = da.acceptorRemote_;
            acceptorVertexIDisFinal_ = da.acceptorVertexIDisFinal_;
            acceptorCandidateCells_ = da.acceptorCandidateCells_;
            acceptorCandidateCellsIDs_ = da.acceptorCandidateCellsIDs_;
            acceptorCandidateCellsPoints_ = da.acceptorCandidateCellsPoints_;
            donorCell_ = da.donorCell_;
            donorProcNo_ = da.donorProcNo_;
            donorPoint_ = da.donorPoint_;
            mainDonorCellVertices_ = da.mainDonorCellVertices_;
            mainDonorCellPoints_ = da.mainDonorCellPoints_;
            donorRegionID_ = da.donorRegionID_;
            extendedDonorCells_ = da.extendedDonorCells_;
            extendedDonorPoints_ = da.extendedDonorPoints_;
            withinBB_ = da.withinBB_;
        }
        
         //- Copy assignement
        inline donorAcceptorTioga& operator=(const donorAcceptorTioga& da)
        {
            if (this == &da)
            {
                return *this;
            }
             
            acceptorCell_ = da.acceptorCell_;
            acceptorProcNo_ = da.acceptorProcNo_;
            acceptorPoint_ = da.acceptorPoint_;
            acceptorRegionID_ = da.acceptorRegionID_;
            acceptorVertexID_ = da.acceptorVertexID_*2;
            acceptorRemote_ = da.acceptorRemote_;
            acceptorVertexIDisFinal_ = da.acceptorVertexIDisFinal_;
            acceptorCandidateCells_ = da.acceptorCandidateCells_;
            acceptorCandidateCellsIDs_ = da.acceptorCandidateCellsIDs_;
            acceptorCandidateCellsPoints_ = da.acceptorCandidateCellsPoints_;            
            donorCell_ = da.donorCell_;
            donorProcNo_ = da.donorProcNo_;
            donorPoint_ = da.donorPoint_;
            mainDonorCellVertices_ = da.mainDonorCellVertices_;
            mainDonorCellPoints_ = da.mainDonorCellPoints_;
            donorRegionID_ = da.donorRegionID_;
            extendedDonorCells_ = da.extendedDonorCells_;
            extendedDonorPoints_ = da.extendedDonorPoints_;
            withinBB_ = da.withinBB_;
            
            return *this;
        }
        
    inline donorAcceptorTioga& operator=(donorAcceptorTioga&& da) noexcept
    {
        if (this == &da)
        {
            return *this;
        }

        // Move data from da to *this
        acceptorCell_ = std::move(da.acceptorCell_);
        acceptorProcNo_ = std::move(da.acceptorProcNo_);
        acceptorPoint_ = std::move(da.acceptorPoint_);
        acceptorRegionID_ = std::move(da.acceptorRegionID_);
        acceptorVertexID_ = std::move(da.acceptorVertexID_);
        acceptorRemote_ = std::move(da.acceptorRemote_);
        acceptorVertexIDisFinal_ = std::move(da.acceptorVertexIDisFinal_);
        acceptorCandidateCells_ = std::move(da.acceptorCandidateCells_);
        acceptorCandidateCellsIDs_ = std::move(da.acceptorCandidateCellsIDs_);
        acceptorCandidateCellsPoints_ = std::move(da.acceptorCandidateCellsPoints_);       
        donorCell_ = std::move(da.donorCell_);
        donorProcNo_ = std::move(da.donorProcNo_);
        donorPoint_ = std::move(da.donorPoint_);
        mainDonorCellVertices_ = std::move(da.mainDonorCellVertices_);
        mainDonorCellPoints_ = std::move(da.mainDonorCellPoints_);
        donorRegionID_ = std::move(da.donorRegionID_);
        extendedDonorCells_ = std::move(da.extendedDonorCells_);
        extendedDonorPoints_ = std::move(da.extendedDonorPoints_);
        withinBB_ = std::move(da.withinBB_);

        return *this;
    }


    //- Destructor
    ~donorAcceptorTioga() {}


    // Member Functions

        // Acceptor side

            //- Return acceptor cell number
            label acceptorCell() const
            {
                return acceptorCell_;
            }

            //- Return access to acceptor cell number
            label& acceptorCell()
            {
                return acceptorCell_;
            }

            //- Return acceptor processor number
            label acceptorProcNo() const
            {
                return acceptorProcNo_;
            }

            //- Return access to acceptor processor number
            label& acceptorProcNo()
            {
                return acceptorProcNo_;
            }

            //- Return acceptor location
            const point& acceptorPoint() const
            {
                return acceptorPoint_;
            }

            //- Return access to acceptor location
            point& acceptorPoint()
            {
                return acceptorPoint_;
            }
            
////
            //- Return acceptor region ID
            label acceptorRegionID() const
            {
                return acceptorRegionID_;
            }

            //- Return access to acceptor region ID
            label& acceptorRegionID()
            {
                return acceptorRegionID_;
            }
            
            //- Return acceptor vertex ID
            label acceptorVertexID() const
            {
                return acceptorVertexID_;
            }

            //- Return access to vertex ID
            label& acceptorVertexID()
            {
                return acceptorVertexID_;
            }
            
            //- Return acceptor is at remote processor
            bool acceptorRemote() const
            {
                return acceptorRemote_;
            }
            
            //- Set access acceptor remote processor
            void setAcceptorRemote(const bool acceptorRemote)
            {
                acceptorRemote_ = acceptorRemote;
            }

            //- Return acceptor vertex ID is final
            bool acceptorVertexIDisFinal() const
            {
                return acceptorVertexIDisFinal_;
            }
            
            //- Set access acceptor vertex ID is final
            void SetAcceptorVertexIDisFinal(const bool acceptorVertexIDisFinal)
            {
                acceptorVertexIDisFinal_ = acceptorVertexIDisFinal;
            } 
            
            //- Return acceptor vertex ID
            label acceptorCandidateCells() const
            {
                return acceptorCandidateCells_;
            }

            //- Return access to vertex ID
            label& acceptorCandidateCells()
            {
                return acceptorCandidateCells_;
            }            
        
             //- Return acceptor Candidate cells IDs
            const DynamicLabelList& acceptorCandidateCellsIDs() const
            {
                return acceptorCandidateCellsIDs_;
            }

            //- Return access to acceptor Candidate cells IDs
            DynamicLabelList& acceptorCandidateCellsIDs()
            {
                return acceptorCandidateCellsIDs_;
            }

            //- Return main acceptor Candidate cells center coordinates
            const DynamicPointList& acceptorCandidateCellsPoints() const
            {
                return acceptorCandidateCellsPoints_;
            }

            //- Return access to acceptor cell vertices coordinates
            DynamicPointList& acceptorCandidateCellsPoints()
            {
                return acceptorCandidateCellsPoints_;
            }
            
            label bestAcceptorCandidate(scalar& ds) const
            {
                std::vector<double> d(acceptorCandidateCells_);
                for(int q=0; q<acceptorCandidateCells_; q++)
                {
                    d[q] = mag(acceptorCandidateCellsPoints_[q] - donorPoint_);
                }
                                                                        
                auto minElementIter = std::min_element(d.begin(), d.end());
                size_t minIndex = std::distance(d.begin(), minElementIter);
                
                ds = d[minIndex];
                
                //Info<<"acceptorRemote(): "<<acceptorRemote_<<endl;
                //Info<<"smallest distance: "<<d[minIndex]<<endl;
                //Perr<< "smallest distance" << d[minIndex] <<endl;
                
                //const bool withinBB =  mesh_.pointInCellBB
                //(
                    //curP,
                    //curDonors[donorCandidateIndex]
                //);
                
                return acceptorCandidateCellsIDs_[minIndex]; 
            }                     
////            


        // Donor side

            //- Has a donor been found?
            bool donorFound() const
            {
                return donorCell_ > -1;
            }

            //- Return donor cell number
            label donorCell() const
            {
                return donorCell_;
            }

            //- Return access to donor cell number
            label& donorCell()
            {
                return donorCell_;
            }

            //- Return donor processor number
            label donorProcNo() const
            {
                return donorProcNo_;
            }

            //- Return access to donor processor number
            label& donorProcNo()
            {
                return donorProcNo_;
            }

/*
            //- Return donor point
            const point& donorPoint() const
            {
                if (donorFound())
                {
                    return donorPoint_;
                }
                else
                {
                    return vector::max;
                }
            }
*/
            //- Return access to donor point
            point& donorPoint()
            {
                return donorPoint_;
            }
////
            //- Return main donor cell vertex IDs
            const DynamicLabelList& mainDonorCellVertices() const
            {
                return mainDonorCellVertices_;
            }

            //- Return access to main donor cell vertex IDs
            DynamicLabelList& mainDonorCellVertices()
            {
                return mainDonorCellVertices_;
            }

            //- Return main donor cell vertices coordinates
            const DynamicPointList& mainDonorCellPoints() const
            {
                return mainDonorCellPoints_;
            }

            //- Return access to donor cell vertices coordinates
            DynamicPointList& mainDonorCellPoints()
            {
                return mainDonorCellPoints_;
            }
            
            //- Return donor region ID
            label donorRegionID() const
            {
                return donorRegionID_;
            }

            //- Return access to region ID
            label& donorRegionID()
            {
                return donorRegionID_;
            }
            
////

            //- Return extended donor cell numbers
            const DynamicLabelList& extendedDonorCells() const
            {
                return extendedDonorCells_;
            }

            //- Return access to extended donor cell numbers
            DynamicLabelList& extendedDonorCells()
            {
                return extendedDonorCells_;
            }

            //- Return extended donor cell centres
            const DynamicPointList& extendedDonorPoints() const
            {
                return extendedDonorPoints_;
            }

            //- Return access to extended donor cell centres
            DynamicPointList& extendedDonorPoints()
            {
                return extendedDonorPoints_;
            }

    // Friend Operators

        friend bool operator==
        (
            const donorAcceptorTioga& a,
            const donorAcceptorTioga& b
        )
        {
            return
                a.acceptorCell_ == b.acceptorCell_
             && a.acceptorProcNo_ == b.acceptorProcNo_
             && a.acceptorPoint_ == b.acceptorPoint_
             && a.donorCell_ == b.donorCell_
             && a.donorProcNo_ == b.donorProcNo_
             && a.withinBB_ == b.withinBB_;

             // Note: do not check whether extended neighbours are the same, we
             // assume they will be if donor data is the same
        }

        friend bool operator!=
        (
            const donorAcceptorTioga& a,
            const donorAcceptorTioga& b
        )
        {
            return !(a == b);
        }


    // IOstream Operators

        //- Istream operator
        friend Istream& operator>>(Istream& is, donorAcceptorTioga& rd)
        {
            is  >> rd.acceptorCell_
                >> rd.acceptorProcNo_
                >> rd.acceptorPoint_
                >> rd.acceptorRegionID_
                >> rd.acceptorVertexID_
                >> rd.acceptorRemote_
                >> rd.acceptorVertexIDisFinal_
                >> rd.acceptorCandidateCells_
                >> rd.acceptorCandidateCellsIDs_
                >> rd.acceptorCandidateCellsPoints_  
                >> rd.donorCell_
                >> rd.donorProcNo_
                >> rd.donorPoint_
                >> rd.mainDonorCellVertices_
                >> rd.mainDonorCellPoints_
                >> rd.donorRegionID_
                >> rd.extendedDonorCells_
                >> rd.extendedDonorPoints_
                >> rd.withinBB_;
            return is;
        }

        //- Ostream operator
        friend Ostream& operator<<(Ostream& os, const donorAcceptorTioga& rd)
        {
            os  << rd.acceptorCell_ << token::SPACE
                << rd.acceptorProcNo_ << token::SPACE
                << rd.acceptorPoint_ << token::SPACE
                << rd.acceptorRegionID_ << token::SPACE
                << rd.acceptorVertexID_ << token::SPACE
                << rd.acceptorRemote_ << token::SPACE
                << rd.acceptorVertexIDisFinal_ << token::SPACE
                << rd.acceptorCandidateCells_ << token::SPACE
                << rd.acceptorCandidateCellsIDs_ << token::SPACE
                << rd.acceptorCandidateCellsPoints_ << token::SPACE         
                << rd.donorCell_ << token::SPACE
                << rd.donorProcNo_ << token::SPACE
                << rd.donorPoint_ << token::SPACE
                << rd.mainDonorCellVertices_ << token::SPACE
                << rd.mainDonorCellPoints_ << token::SPACE
                << rd.donorRegionID_ << token::SPACE
                << rd.extendedDonorCells_ << token::SPACE
                << rd.extendedDonorPoints_ << token::SPACE
                << rd.withinBB_ << token::SPACE;                               
            return os;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
