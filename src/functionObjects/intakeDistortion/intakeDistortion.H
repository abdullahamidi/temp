/*---------------------------------------------------------------------------*\

    OSCFD Group

-------------------------------------------------------------------------------
Class
    intakeDistortion

Description

Authors
    Barış Cumhur
        Turkish Aerospace Industries, Turkey

Description
    Calculates the distortion coefficient that is described by

    DC\left(\alpha\right) = \frac{P_{t, mean} - P_{t, mean, \left(\alpha\right)}}{q}

    where $q$ is the mean dynamic pressure

    Calculation steps are listed below
    \begin{itemize}
    1) Calculate patch centroid, normal, tangent and binormal
    2) Calculate the mean total pressure and mean dynamic pressure over patch
    3) Select a plane using the binormal and centroid of the patch as 'reference' plane.
    4) Rotate the reference plane by $\alpha$ to define the 'interval' plane.
    5) For each angular position incremented by a scecified stride
        a) Rotate reference and interval planes to describe current angular position and interval
        b) Cut the patch with reference and interval planes
        c) Calculate area weighted average on cutted patch
        d) Store the results
    \end{itemize}


\*---------------------------------------------------------------------------*/

#ifndef functionObjects_intakeDistortion_H
#define functionObjects_intakeDistortion_H

#include "fvMeshFunctionObject.H"
#include "writeFile.H"
#include "Vector.H"
#include "vectorField.H"
#include "surfaceMesh.H"
#include "fvMesh.H"
#include "volFields.H"
#include "fvPatchFields.H"
#include "fvsPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam {
namespace functionObjects {

/*---------------------------------------------------------------------------*\
                       Class intakeDistortion Declaration
\*---------------------------------------------------------------------------*/

class intakeDistortion
    : public fvMeshFunctionObject, writeFile {

    //- Region type enumeration
    enum regionTypes {
        stFaceZone = 0x01,   //!< Calculate with faceZone(s)
        stPatch = 0x02,   //!< Calculate with patch(es)
    };

    enum avgTypes {
        areaAverage = 0x01,
        massAverage = 0x02,
    };

    static const Enum<regionTypes> region_type_names_;
    static const Enum<avgTypes> average_type_names_;

    regionTypes region_type_;
    avgTypes average_type_;
    word region_name_;
    wordRes selection_names_;
    bool calculate_intersections_;

    label n_faces_;
    labelList face_ids_;
    labelList face_patch_ids_;
    boolList face_flip_map_;

    scalar angle_interval_;
    scalar angle_stride_;
    label angle_count_;

    vector patch_surface_;
    scalar patch_area_;
    vector patch_normal_;
    vector patch_centroid_;
    vector patch_tangent_;
    vector patch_binormal_;

    scalar sum_mag_mass_flow_rate_;
    scalar avg_total_pressure_;
    scalar avg_dynamic_pressure_;

    scalarField avg_total_pressure_by_angle_;
    scalarField area_by_angle_;
    scalarField angle_;

    //- Output file header information
    virtual void writeFileHeader(Ostream& os);

    //- No copy construct
    intakeDistortion(const intakeDistortion&) = delete;

    //- No copy assignment
    void operator=(const intakeDistortion&) = delete;

    void setPatchFaces();

    void setFaceZoneFaces();

    scalar vAngle(const vector& ref, const vector& v);

    bool inside(const vector& ref, const vector& intervalRef, const vector& point);

    bool sideFlip(const vector& ref, const vector& v1, const vector& v2);

    vector project(const vector& v);

    vector intersect(const vector& refNormal, const vector& begin, const vector& end);

    vector surface(const vectorField& vertices);

    const face& selectedFace(label localIndex) const;
    
    template<class Type>
    Foam::Field<Type> filterField
    (
        const GeometricField<Type, fvPatchField, volMesh>& field
    ) const {
        const labelList& own = field.mesh().faceOwner();
        const labelList& nei = field.mesh().faceNeighbour();

        auto values = Field<Type>(face_ids_.size());

        forAll(values, i) {
            const label facei = face_ids_[i];
            const label patchi = face_patch_ids_[i];

            if (patchi >= 0) {
                // Boundary face - face id is the patch-local face id
                values[i] = field.boundaryField()[patchi][facei];
            } else {
                // Internal face
                values[i] = 0.5 * (field[own[facei]] + field[nei[facei]]);
            }
        }

        // No need to flip values - all boundary faces point outwards

        return values;
    }

    template<class Type>
    Foam::Field<Type> filterField
    (
        const GeometricField<Type, fvsPatchField, surfaceMesh>& field
    ) const {

        Field<Type> values(face_ids_.size());

        forAll(values, i) {
            const label facei = face_ids_[i];
            const label patchi = face_patch_ids_[i];

            if (patchi >= 0) {
                values[i] = field.boundaryField()[patchi][facei];
            } else {
                values[i] = field[facei];
            }
        }

        if (field.oriented()()) {
            forAll(values, i) {
                if (face_flip_map_[i]) {
                    values[i] *= -1;
                }
            }
        }

        return values;
    }

public:

    //- Runtime type information
    TypeName("intakeDistortion");


    // Constructors

    //- Construct from Time and dictionary
    intakeDistortion
    (
        const word& name,
        const Time& runTime,
        const dictionary& dict
    );


    //- Destructor
    virtual ~intakeDistortion() = default;

    // Member Functions

    //- Read the necessary data
    virtual bool read(const dictionary&);

    //- Execute, currently does nothing
    virtual bool execute();

    //- Write the intakeDistortion
    virtual bool write();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
