/*---------------------------------------------------------------------------*\

    OSCFD Group

-------------------------------------------------------------------------------
Class
    soundPressureLevel

Description

Authors
    Barış Cumhur
        Turkish Aerospace Industries, Turkey

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_soundPressureLevel_H
#define functionObjects_soundPressureLevel_H

#include "fvMeshFunctionObject.H"
#include "writeFile.H"
#include "Vector.H"
#include "vectorField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam {
namespace functionObjects {

/*---------------------------------------------------------------------------*\
                       Class soundPressureLevel Declaration
\*---------------------------------------------------------------------------*/

class soundPressureLevel
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
    soundPressureLevel(const soundPressureLevel&) = delete;

    //- No copy assignment
    void operator=(const soundPressureLevel&) = delete;

    void read_probe();
public:

    //- Runtime type information
    TypeName("soundPressureLevel");


    // Constructors

    //- Construct from Time and dictionary
    soundPressureLevel
    (
        const word& name,
        const Time& runTime,
        const dictionary& dict
    );


    //- Destructor
    virtual ~soundPressureLevel() = default;


    // Member Functions

    //- Read the necessary data
    virtual bool read(const dictionary&);

    //- Execute, currently does nothing
    virtual bool execute();

    //- Write the soundPressureLevel
    virtual bool write();

    static void do_calculate(const Field<vectorField>& vdata, label data_count, scalar delta_t, label window_count, scalar window_overlap, fileName output);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
