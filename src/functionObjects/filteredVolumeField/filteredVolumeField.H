/*---------------------------------------------------------------------------*\

    V7320 - Özgün Mühendislik Yazılımları

-------------------------------------------------------------------------------
Class
    FilteredVolumeField

Description

Authors
    Barış Cumhur
        Turkish Aerospace Industries, Turkey

Description
    Exports filtered field values from a given region in the mesh. Filtering
    can be done based on bounding box, wall distance, nut, nutilda, k, omega.
    The filtered values are written to a CSV file.

\*---------------------------------------------------------------------------*/

#ifndef FILTERED_FIELD_VALUE_H
#define FILTERED_FIELD_VALUE_H

#include "Field.H"
#include "GeometricField.H"
#include "IOobject.H"
#include "Pstream.H"
#include "autoPtr.H"
#include "boundBox.H"
#include "fvMeshFunctionObject.H"
#include "fvPatchField.H"
#include "labelIOList.H"
#include "volFieldsFwd.H"
#include "volMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam {
namespace functionObjects {

/*---------------------------------------------------------------------------*\
                       Class intakeDistortion Declaration
\*---------------------------------------------------------------------------*/

class filteredVolumeField : public fvMeshFunctionObject {

  word field_;

  bool has_bounding_box_;
  bool pending_bounding_box_filter_;
  boundBox bounding_box_;

  bool has_y_filter_;
  bool pending_y_filter_;
  scalar y_min_;
  scalar y_max_;

  bool has_nut_filter_;
  scalar nut_min_;
  scalar nut_max_;

  bool has_nutilda_filter_;
  scalar nutilda_min_;
  scalar nutilda_max_;

  bool has_k_filter_;
  scalar k_min_;
  scalar k_max_;

  bool has_omega_filter_;
  scalar omega_min_;
  scalar omega_max_;

  volScalarField y_;
  labelIOList cell_proc_addressing_;
  labelIOList label_list_;

  //- No copy construct
  filteredVolumeField(filteredVolumeField const &) = delete;

  //- No copy assignment
  void operator=(filteredVolumeField const &) = delete;

  bool has_static_filter() const;

  bool has_dynamic_filter() const;

  bool has_filter() const;

  bool requires_static_update() const;

  void compute_cells();

public:
  //- Runtime type information
  TypeName("filteredVolumeField");

  // Constructors

  //- Construct from Time and dictionary
  filteredVolumeField(word const &name, Time const &runTime,
                      dictionary const &dict);

  //- Destructor
  virtual ~filteredVolumeField();

  // Member Functions

  //- Read the necessary data
  virtual bool read(dictionary const &) override;

  //- Execute, currently does nothing
  virtual bool execute() override;

  //- Write the intakeDistortion
  virtual bool write() override;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
