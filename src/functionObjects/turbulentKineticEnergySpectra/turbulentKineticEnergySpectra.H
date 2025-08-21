/*---------------------------------------------------------------------------*\

    OSCFD Group

-------------------------------------------------------------------------------
Class
    turbulentKineticEnergySpectra

Description

Authors
    Barış Cumhur
        Turkish Aerospace Industries, Turkey

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_turbulentKineticEnergySpectra_H
#define functionObjects_turbulentKineticEnergySpectra_H

#include "fvMeshFunctionObject.H"
#include "writeFile.H"
#include "Vector.H"
#include "vectorField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam {
namespace functionObjects {

/*---------------------------------------------------------------------------*\
                       Class turbulentKineticEnergySpectra Declaration
\*---------------------------------------------------------------------------*/

class turbulentKineticEnergySpectra
    : public fvMeshFunctionObject, writeFile {
protected:
    word UProbeName_;

    scalar begin_time_;

    label window_count_;

    scalar window_overlap_;

    Field<vectorField> vdata_;

    scalar delta_t_;

    label data_count_;

    fileName probeDir_;

    //- Output file header information
    virtual void writeFileHeader(Ostream& os);

    //- No copy construct
    turbulentKineticEnergySpectra(const turbulentKineticEnergySpectra&) = delete;

    //- No copy assignment
    void operator=(const turbulentKineticEnergySpectra&) = delete;

    void read_probe();
public:

    //- Runtime type information
    TypeName("turbulentKineticEnergySpectra");


    // Constructors

    //- Construct from Time and dictionary
    turbulentKineticEnergySpectra
    (
        const word& name,
        const Time& runTime,
        const dictionary& dict
    );


    //- Destructor
    virtual ~turbulentKineticEnergySpectra() = default;


    // Member Functions

    //- Read the necessary data
    virtual bool read(const dictionary&);

    //- Execute, currently does nothing
    virtual bool execute();

    //- Write the turbulentKineticEnergySpectra
    virtual bool write();

    static void do_calculate(const Field<vectorField>& vdata, label data_count, scalar delta_t, label window_count, scalar window_overlap, fileName output);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
