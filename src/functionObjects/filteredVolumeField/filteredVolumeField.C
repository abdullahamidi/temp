/*---------------------------------------------------------------------------*\

    V7320 - Özgün Mühendislik Yazılımları

\*---------------------------------------------------------------------------*/

#include "IOList.H"
#include "IOobject.H"
#include "ISstream.H"
#include "addToRunTimeSelectionTable.H"
#include "boundBox.H"
#include "csv.H"
#include "filteredVolumeField.H"
#include "fvMesh.H"
#include "globalIndex.H"
#include "labelIOList.H"
#include "messageStream.H"
#include "scalarFwd.H"
#include "stdFoam.H"
#include "vectorIOField.H"
#include "volFieldsFwd.H"
#include "wallDist.H"
#include <utility>

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam {
namespace functionObjects {

defineTypeNameAndDebug(filteredVolumeField, 0);
addToRunTimeSelectionTable(functionObject, filteredVolumeField, dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

filteredVolumeField::filteredVolumeField(word const &name, Time const &runTime,
                                         dictionary const &dict)
    : fvMeshFunctionObject(name, runTime, dict),
      y_(name + "_y", wallDist(mesh_).y()),
      cell_proc_addressing_(IOobject(mesh_.meshDir() + "/cellProcAddressing",
                                     mesh_.facesInstance(), mesh_,
                                     IOobject::MUST_READ, IOobject::NO_WRITE)),
      label_list_(IOobject(name + "_local_indices", mesh_.time().timeName(),
                           mesh_, IOobject::NO_READ, IOobject::NO_WRITE)) {

  Info << "filteredFieldValue " << name << nl;

  read(dict);
}

filteredVolumeField::~filteredVolumeField() {}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool filteredVolumeField::read(dictionary const &dict) {
  if (fvMeshFunctionObject::read(dict)) {
    field_ = dict.getOrDefault<word>("field", "p");

    Info << "    Bounding box:  ";
    has_bounding_box_ = dict.getOrDefault<bool>("enableBoundingBox", false);
    pending_bounding_box_filter_ = true;
    if (has_bounding_box_) {
      Info << " Enabled";
      vector bounding_box_lower_left_ =
          dict.get<vector>("boundingBoxLowerLeft");
      Info << ", Min: " << bounding_box_lower_left_;
      vector bounding_box_upper_right_ =
          dict.get<vector>("boundingBoxUpperRight");
      Info << ", Max: " << bounding_box_upper_right_;
      bounding_box_ =
          boundBox(bounding_box_lower_left_, bounding_box_upper_right_);
    } else {
      Info << " Disabled";
    }
    Info << endl;

    Info << "    y filter:      ";
    has_y_filter_ = dict.getOrDefault<bool>("enableYFilter", false);
    pending_y_filter_ = true;
    if (has_y_filter_) {
      Info << " Enabled";
      y_min_ = dict.get<scalar>("minY");
      Info << ", Min: " << y_min_;
      y_max_ = dict.get<scalar>("maxY");
      Info << ", Max: " << y_max_;
    } else {
      Info << " Disabled";
    }
    Info << endl;

    Info << "    nut filter:    ";
    has_nut_filter_ = dict.getOrDefault<bool>("enableNutFilter", false);
    if (has_nut_filter_) {
      Info << " Enabled";
      nut_min_ = dict.get<scalar>("minNut");
      Info << ", Min: " << nut_min_;
      nut_max_ = dict.get<scalar>("maxNut");
      Info << ", Max: " << nut_max_;
    } else {
      Info << " Disabled";
    }
    Info << endl;

    Info << "    nuTilda filter:";
    has_nutilda_filter_ = dict.getOrDefault<bool>("enableNuTildaFilter", false);
    if (has_nutilda_filter_) {
      Info << " Enabled";
      nutilda_min_ = dict.get<scalar>("minNuTilda");
      Info << ", Min:  " << nutilda_min_;
      nutilda_max_ = dict.get<scalar>("maxNuTilda");
      Info << ", Max:  " << nutilda_max_;
    } else {
      Info << " Disabled";
    }
    Info << endl;

    Info << "    k filter:      ";
    has_k_filter_ = dict.getOrDefault("enableKFilter", false);
    if (has_k_filter_) {
      Info << " Enabled";
      k_min_ = dict.get<scalar>("minK");
      Info << ", Min:  " << k_min_;
      k_max_ = dict.get<scalar>("maxK");
      Info << ", Max:  " << k_max_;
    } else {
      Info << " Disabled";
    }
    Info << endl;

    Info << "    omega filter:  ";
    has_omega_filter_ = dict.getOrDefault("enableOmegaFilter", false);
    if (has_omega_filter_) {
      Info << " Enabled";
      omega_min_ = dict.get<scalar>("minOmega");
      Info << ", Min:  " << omega_min_;
      omega_max_ = dict.get<scalar>("maxOmega");
      Info << ", Max:  " << omega_max_;
    } else {
      Info << " Disabled";
    }
    Info << endl << endl;

    return true;
  }

  return false;
}

bool filteredVolumeField::has_static_filter() const {
  return has_bounding_box_ || has_y_filter_;
}

bool filteredVolumeField::has_dynamic_filter() const {
  return has_nut_filter_ || has_k_filter_ || has_omega_filter_;
}

bool filteredVolumeField::has_filter() const {
  return has_static_filter() || has_dynamic_filter();
}

bool filteredVolumeField::requires_static_update() const {
  return has_static_filter() &&
         (pending_bounding_box_filter_ || pending_y_filter_);
}

void filteredVolumeField::compute_cells() {
  label_list_.clear();

  vector const *c_begin = nullptr;
  vector const *c_end = c_begin - 1;
  scalar const *y_begin = nullptr;
  scalar const *y_end = y_begin - 1;
  scalar const *nut_begin = nullptr;
  scalar const *nut_end = nut_begin - 1;
  scalar const *nutilda_begin = nullptr;
  scalar const *nutilda_end = nutilda_begin - 1;
  scalar const *k_begin = nullptr;
  scalar const *k_end = k_begin - 1;
  scalar const *omega_begin = nullptr;
  scalar const *omega_end = omega_begin - 1;

  if (has_bounding_box_) {
    c_begin = mesh_.C().begin();
    c_end = mesh_.C().end();
  }

  if (has_y_filter_) {
    y_begin = y_.begin();
    y_end = y_.end();
  }
  if (has_nut_filter_) {
    auto &field = lookupObjectRef<volScalarField>("nut");
    nut_begin = field.begin();
    nut_end = field.end();
  }

  if (has_nutilda_filter_) {
    auto &field = lookupObjectRef<volScalarField>("nuTilda");
    nutilda_begin = field.begin();
    nutilda_end = field.end();
  }

  if (has_k_filter_) {
    auto &field = lookupObjectRef<volScalarField>("k");
    nutilda_begin = field.begin();
    nutilda_end = field.end();
  }

  if (has_omega_filter_) {
    auto &field = lookupObjectRef<volScalarField>("omega");
    nutilda_begin = field.begin();
    nutilda_end = field.end();
  }

  label idx = 0;
  vector const *c = c_begin;
  double const *y = y_begin;
  double const *nut = nut_begin;
  double const *nutilda = nutilda_begin;
  double const *k = k_begin;
  double const *omega = omega_begin;

  for (; c != c_end && y != y_end && nut != nut_end && nutilda != nutilda_end &&
         k != k_end && omega != omega_end && idx < mesh_.nCells();
       ++c, ++y, ++nut, ++nutilda, ++k, ++omega, ++idx) {

    bool includeCell = true;

    if (includeCell && has_bounding_box_) {
      vector value = *c;
      includeCell = (bounding_box_.min().x() <= value.x() &&
                     value.x() <= bounding_box_.max().x() &&
                     bounding_box_.min().y() <= value.y() &&
                     value.y() <= bounding_box_.max().y() &&
                     bounding_box_.min().z() <= value.z() &&
                     value.z() <= bounding_box_.max().z());
    }

    if (includeCell && has_y_filter_) {
      scalar value = *y;
      includeCell = y_min_ <= value && value <= y_max_;
    }

    if (includeCell && has_nut_filter_) {
      scalar value = *nut;
      includeCell = nut_min_ <= value && value <= nut_max_;
    }

    if (includeCell && has_k_filter_) {
      scalar value = *k;
      includeCell = k_min_ <= value && value <= k_max_;
    }

    if (includeCell && has_omega_filter_) {
      scalar value = *omega;
      includeCell = omega_min_ <= value && value <= omega_max_;
    }

    if (includeCell) {
      label_list_.append(idx);
    }
  }

  pending_bounding_box_filter_ = false;
  pending_y_filter_ = false;
}

bool filteredVolumeField::execute() {
  if (!has_filter()) {
    Warning << " " << name()
            << ": No filter applied, keeping last computed cell label list "
               "unchanged"
            << endl;
    return true;
  }

  if (!has_dynamic_filter() && !requires_static_update()) {
    return true;
  }

  compute_cells();

  return true;
}

bool filteredVolumeField::write() {
  if (foundObject<volScalarField>(field_) ||
      foundObject<volVectorField>(field_)) {
    Foam::writeCsv(IOobject(name() + ".csv", mesh_.time().timeName(), mesh_,
                            IOobject::NO_READ, IOobject::NO_WRITE),
                   label_list_, false, std::make_pair("index", &label_list_),
                   std::make_pair("global_index", &cell_proc_addressing_),
                   std::make_pair("centroid", &mesh_.C()),
                   std::make_pair("volume", &mesh_.V()),
                   std::make_pair(field_.c_str(),
                                  foundObject<volScalarField>(field_)
                                      ? &lookupObjectRef<volScalarField>(field_)
                                      : nullptr),
                   std::make_pair(field_.c_str(),
                                  foundObject<volVectorField>(field_)
                                      ? &lookupObjectRef<volVectorField>(field_)
                                      : nullptr));
  } else {
    Info << " Can't find field: " << field_ << endl;
    return false;
  }

  return true;
}
} // namespace functionObjects
} // namespace Foam

// ************************************************************************* //
