/*---------------------------------------------------------------------------*\

    OSCFD Group

\*---------------------------------------------------------------------------*/

#include "IFstream.H"
#include "IOmanip.H"
#include "OFstream.H"
#include "addToRunTimeSelectionTable.H"
#include "fvMesh.H"
#include "fvPatchFields.H"
#include "fvsPatchFields.H"
#include "intakeDistortion.H"
#include "mathematicalConstants.H"
#include "surfaceMesh.H"
#include "volFields.H"
#include "writeCSV.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam {
namespace functionObjects {

defineTypeNameAndDebug(intakeDistortion, 0);
addToRunTimeSelectionTable(functionObject, intakeDistortion, dictionary);

const Foam::Enum<intakeDistortion::regionTypes>
    intakeDistortion::region_type_names_({
        {regionTypes::stFaceZone, "faceZone"},
        {regionTypes::stPatch, "patch"},
    });

const Foam::Enum<intakeDistortion::avgTypes>
    intakeDistortion::average_type_names_({
        {avgTypes::areaAverage, "areaAverage"},
        {avgTypes::massAverage, "massAverage"},
    });

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void intakeDistortion::writeFileHeader(Ostream &os) {
  writeHeader(os, "intake Pressure Distortion Calculator");
  writeCommented(os, "Distorsion coefficient by angle");

  os << endl;
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

intakeDistortion::intakeDistortion(const word &name, const Time &runTime,
                                   const dictionary &dict)
    : fvMeshFunctionObject(name, runTime, dict),
      writeFile(mesh_, name, typeName, dict), n_faces_(0) {

  Info << "intakeDistortion " << name << nl;
  read(dict);

  writeFileHeader(file());
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool intakeDistortion::read(const dictionary &dict) {
  fvMeshFunctionObject::read(dict);
  writeFile::read(dict);

  face_ids_.clear();
  face_patch_ids_.clear();
  face_flip_map_.clear();
  region_name_.clear();
  selection_names_.clear();

  angle_stride_ = dict.getOrDefault<scalar>("angularStride", 1) / 180.0 *
                  constant::mathematical::pi;
  angle_interval_ = dict.getOrDefault<scalar>("angularInterval", 60) / 180.0 *
                    constant::mathematical::pi;
  calculate_intersections_ =
      dict.getOrDefault<bool>("calculateIntersections", true);

  angle_count_ = label(constant::mathematical::twoPi / angle_stride_);
  angle_stride_ = constant::mathematical::twoPi / angle_count_;

  Info << "    Distorsion coefficients going to be for " << nl << endl;
  Info << "    stride :" << angle_stride_ << nl << endl;
  Info << "    interval :" << angle_interval_ << nl << endl;

  region_type_ = region_type_names_.get("regionType", dict);
  average_type_ = average_type_names_.get("averageType", dict);

  {
    dict.readIfPresent("names", selection_names_);

    for (const auto &item : selection_names_) {
      if (item.isLiteral()) {
        region_name_ = item;
        break;
      }
    }

    // Mandatory if we didn't pick up a value from selectionNames_
    dict.readEntry("name", region_name_, keyType::LITERAL,
                   region_name_.empty()
                       ? IOobjectOption::readOption::MUST_READ
                       : IOobjectOption::readOption::READ_IF_PRESENT);

    // Ensure there is always content for selectionNames_
    if (selection_names_.empty()) {
      selection_names_.resize(1);
      selection_names_.first() = region_name_;
    }
  }

  switch (region_type_) {
  case stFaceZone: {
    setFaceZoneFaces();
    break;
  }
  case stPatch: {
    setPatchFaces();
    break;
  }
    // Compiler warning if we forgot an enumeration
  }

  if (n_faces_ == 0) {
    FatalErrorInFunction << type() << ' ' << name() << ": "
                         << region_type_names_[region_type_] << '('
                         << region_name_ << "):" << nl
                         << "    Region has no faces" << exit(FatalError);
  }

  Info << nl << endl;

  return true;
}

void intakeDistortion::setFaceZoneFaces() {
  // Indices for all matches, already sorted
  const labelList zoneIds(mesh_.faceZones().indices(selection_names_));

  // Total number of faces selected
  label numFaces = 0;
  for (const label zoneId : zoneIds) {
    numFaces += mesh_.faceZones()[zoneId].size();
  }

  if (zoneIds.empty()) {
    FatalErrorInFunction << type() << ' ' << name() << ": "
                         << region_type_names_[region_type_] << '('
                         << region_name_ << "):" << nl
                         << "    No matching face zone(s): "
                         << flatOutput(selection_names_) << nl
                         << "    Known face zones: "
                         << flatOutput(mesh_.faceZones().names()) << nl
                         << exit(FatalError);
  }

  face_ids_.resize(numFaces);
  face_patch_ids_.resize(numFaces);
  face_flip_map_.resize(numFaces);

  numFaces = 0;

  for (const label zoneId : zoneIds) {
    const faceZone &fZone = mesh_.faceZones()[zoneId];

    forAll(fZone, i) {
      const label meshFacei = fZone[i];
      const bool isFlip = fZone.flipMap()[i];

      // Internal faces
      label faceId = meshFacei;
      label facePatchId = -1;

      // Boundary faces
      if (!mesh_.isInternalFace(meshFacei)) {
        facePatchId = mesh_.boundaryMesh().whichPatch(meshFacei);
        const polyPatch &pp = mesh_.boundaryMesh()[facePatchId];
        const auto *cpp = isA<coupledPolyPatch>(pp);

        if (cpp) {
          faceId = (cpp->owner() ? pp.whichFace(meshFacei) : -1);
        } else if (!isA<emptyPolyPatch>(pp)) {
          faceId = pp.whichFace(meshFacei);
        } else {
          faceId = -1;
          facePatchId = -1;
        }
      }

      if (faceId >= 0) {
        face_ids_[numFaces] = faceId;
        face_patch_ids_[numFaces] = facePatchId;
        face_flip_map_[numFaces] = isFlip;

        ++numFaces;
      }
    }
  }

  // Shrink to size used
  face_ids_.resize(numFaces);
  face_patch_ids_.resize(numFaces);
  face_flip_map_.resize(numFaces);
  n_faces_ = returnReduce(face_ids_.size(), sumOp<label>());
}

void intakeDistortion::setPatchFaces() {
  // Patch indices for all matches
  labelList patchIds;

  // Total number of faces selected
  label numFaces = 0;

  labelList selected(mesh_.boundaryMesh()
                         .patchSet(selection_names_,
                                   false // warnNotFound - we do that ourselves
                                   )
                         .sortedToc());

  DynamicList<label> bad;
  for (const label patchi : selected) {
    const polyPatch &pp = mesh_.boundaryMesh()[patchi];

    if (isA<emptyPolyPatch>(pp)) {
      bad.append(patchi);
    } else {
      numFaces += pp.size();
    }
  }

  if (bad.size()) {
    label nGood = (selected.size() - bad.size());

    auto &os = (nGood > 0 ? WarningInFunction : FatalErrorInFunction);

    os << "Cannot sample an empty patch" << nl;

    for (const label patchi : bad) {
      os << "    " << mesh_.boundaryMesh()[patchi].name() << nl;
    }

    if (nGood) {
      os << "No non-empty patches selected" << endl << exit(FatalError);
    } else {
      os << "Selected " << nGood << " non-empty patches" << nl;
    }

    patchIds.resize(nGood);
    nGood = 0;
    for (const label patchi : selected) {
      if (!bad.found(patchi)) {
        patchIds[nGood] = patchi;
        ++nGood;
      }
    }
  } else {
    patchIds = std::move(selected);
  }

  if (patchIds.empty()) {
    FatalErrorInFunction << type() << ' ' << name() << ": "
                         << region_type_names_[region_type_] << '('
                         << region_name_ << "):" << nl
                         << "    No matching patch name(s): "
                         << flatOutput(selection_names_) << nl
                         << "    Known patch names:" << nl
                         << mesh_.boundaryMesh().names() << nl
                         << exit(FatalError);
  }

  face_ids_.resize(numFaces);
  face_patch_ids_.resize(numFaces);
  face_flip_map_.resize(numFaces);
  n_faces_ = returnReduce(face_ids_.size(), sumOp<label>());

  numFaces = 0;
  for (const label patchi : patchIds) {
    const polyPatch &pp = mesh_.boundaryMesh()[patchi];
    const label len = pp.size();

    SubList<label>(face_ids_, len, numFaces) = identity(len);
    SubList<label>(face_patch_ids_, len, numFaces) = patchi;
    SubList<bool>(face_flip_map_, len, numFaces) = false;

    numFaces += len;
  }
}

scalar intakeDistortion::vAngle(const vector &ref, const vector &v) {
  scalar angleSide = sign((v ^ ref) & patch_normal_);
  scalar cosAngle = (v / mag(v)) & (ref / mag(ref));
  cosAngle = cosAngle > 0 ? cosAngle - SMALL : cosAngle + SMALL;
  if (cosAngle > 1 || cosAngle < -1)
    Warning << "Cosine can't be greater then 1 or less then -1. Value was "
            << cosAngle << "." << endl;
  return fmod((1 + angleSide) * constant::mathematical::pi -
                  angleSide * acos(cosAngle),
              2 * constant::mathematical::pi);
}

bool intakeDistortion::inside(const vector &refNormal,
                              const vector &intervalRefNormal,
                              const vector &point) {
  return (point & refNormal) >= 0 && (point & intervalRefNormal) <= 0;
}

bool intakeDistortion::sideFlip(const vector &ref, const vector &v1,
                                const vector &v2) {
  return (v1 & ref) * (v2 & ref) < 0;
}

vector intakeDistortion::intersect(const vector &refNormal, const vector &begin,
                                   const vector &end) {
  const vector l = end - begin;
  const scalar denom = stabilise(l & refNormal, VSMALL);
  const scalar t = -(begin & refNormal) / denom;

  if (t < 0 || t > 1) {
    Warning << "line intersects the plane with its extension" << endl;
  }

  return begin + l * t;
}

vector intakeDistortion::surface(const vectorField &vertices) {
  const vector center = sum(vertices) / vertices.size();
  vector s = {0, 0, 0};

  for (label i = 0; i < vertices.size() - 1; ++i)
    s += (vertices[i + 1] - vertices[i]) ^ (center - vertices[i + 1]);

  s += (vertices[0] - vertices[vertices.size() - 1]) ^ (center - vertices[0]);

  return 0.5 * s;
}

const face &intakeDistortion::selectedFace(label localIndex) const {
  const label facei = face_ids_[localIndex];
  const label patchi = face_patch_ids_[localIndex];

  if (patchi >= 0) {
    // Boundary face - face id is the patch-local face id
    return mesh_.boundaryMesh()[patchi][facei];
  }

  return mesh_.faces()[facei];
}

bool intakeDistortion::execute() {
  Info << "intakeDistortion " << name() << nl;
  Info << "   nFace = " << n_faces_ << endl;

  vectorField Sf = filterField(mesh_.Sf());
  patch_surface_ = sum(Sf);
  reduce(patch_surface_, sumOp<vector>());
  Info << "   surface = " << patch_surface_ << endl;
  patch_area_ = mag(patch_surface_);
  Info << "   area = " << patch_area_ << endl;

  patch_normal_ = patch_surface_ / patch_area_;
  patch_tangent_ = {patch_normal_[1], patch_normal_[0], patch_normal_[2]};
  patch_binormal_ = patch_normal_ ^ patch_tangent_;
  patch_binormal_ = patch_binormal_ / mag(patch_binormal_);

  scalarField area(Sf & patch_normal_);
  vectorField Cf = filterField(mesh_.Cf());

  patch_centroid_ = sum(area * Cf);
  reduce(patch_centroid_, sumOp<vector>());
  patch_centroid_ /= patch_area_;

  scalarField p = filterField(lookupObjectRef<volScalarField>("p"));
  vectorField u = filterField(lookupObjectRef<volVectorField>("U"));
  scalarField rho = filterField(lookupObjectRef<volScalarField>("rho"));
  scalarField T = filterField(lookupObjectRef<volScalarField>("T"));
  scalarField magNormU(mag(patch_normal_ & u));
  scalarField pd(0.5 * rho * (u & u));
  scalarField pt(p + pd);

  sum_mag_mass_flow_rate_ = sum(area * rho * magNormU);
  switch (average_type_) {
  case areaAverage:
    avg_total_pressure_ = sum(pt * area);
    avg_dynamic_pressure_ = sum(pd * area);
    break;

  case massAverage:
    avg_total_pressure_ = sum(pt * area * rho * magNormU);
    // avg_dynamic_pressure_ = sum(pd * area * rho * magNormU);
    avg_dynamic_pressure_ = sum(pd * area);
    break;
  }

  reduce(avg_total_pressure_, sumOp<scalar>());
  reduce(avg_dynamic_pressure_, sumOp<scalar>());
  reduce(sum_mag_mass_flow_rate_, sumOp<scalar>());

  switch (average_type_) {
  case areaAverage:
    avg_total_pressure_ = avg_total_pressure_ / patch_area_;
    avg_dynamic_pressure_ = avg_dynamic_pressure_ / patch_area_;
    break;

  case massAverage:

    Info << "   MFR = " << sum_mag_mass_flow_rate_ << endl;
    Info << "   avgTotalPressure_1 = " << avg_total_pressure_ << endl;
    Info << "   avgDynamicPressure_1 = " << avg_dynamic_pressure_ << endl;
    avg_total_pressure_ = avg_total_pressure_ / sum_mag_mass_flow_rate_;
    // avg_dynamic_pressure_ = avg_dynamic_pressure_ / sum_mag_mass_flow_rate_;
    avg_dynamic_pressure_ = avg_dynamic_pressure_ / patch_area_;
    break;
  }

  Info << "   centroid = " << patch_centroid_ << endl;
  Info << "   normal = " << patch_normal_ << endl;
  Info << "   tangent = " << patch_tangent_ << endl;
  Info << "   avgTotalPressure = " << avg_total_pressure_ << endl;
  Info << "   avgDynamicPressure = " << avg_dynamic_pressure_ << endl;

  const vectorField points(mesh_.points() - patch_centroid_);

  vector ref = patch_tangent_;
  vector refNormal = patch_binormal_;
  quaternion rot(patch_normal_, angle_stride_);
  quaternion intervalRot(patch_normal_, angle_interval_);
  vector intervalRefNormal = transform(intervalRot.R(), refNormal).normalise();
  avg_total_pressure_by_angle_ = scalarField(angle_count_);
  area_by_angle_ = scalarField(angle_count_);
  angle_ = scalarField(angle_count_);

  scalarField weighting_mult(area.size());
  switch (average_type_) {
  case areaAverage:
    weighting_mult = 1;
    break;

  case massAverage:
    weighting_mult = rho * magNormU;
    break;
  }

  for (label n = 0; n < angle_count_;
       ++n, ref = transform(rot.R(), ref).normalise(),
             refNormal = transform(rot.R(), refNormal).normalise(),
             intervalRefNormal =
                 transform(intervalRot.R(), refNormal).normalise()) {

    scalar sum_total_pressure = 0;
    scalar sum_area = 0;
    scalar sum_weighting = 0;

    forAll(face_ids_, i) {
      const face &f = selectedFace(i);
      label insideCount = 0;
      forAll(f, j) {
        if (inside(refNormal, intervalRefNormal, points[f[j]])) {
          insideCount++;
        }
      }

      if (insideCount == f.size()) {
        sum_total_pressure =
            sum_total_pressure + weighting_mult[i] * area[i] * pt[i];
        sum_weighting = sum_weighting + weighting_mult[i] * area[i];
        sum_area = sum_area + area[i];
      } else if (calculate_intersections_ && insideCount > 0) {
        vectorField vertices(f.size() + 2);

        vector prevV = points[f[f.size() - 1]];
        label vertexIndex = 0;

        forAll(f, j) {
          const vector &v = points[f[j]];
          if (sideFlip(refNormal, prevV, v)) { // ref
            vertices[vertexIndex++] = intersect(refNormal, prevV, v);
          }
          if (sideFlip(intervalRefNormal, prevV, v)) { // ref + interval
            vertices[vertexIndex++] = intersect(intervalRefNormal, prevV, v);
          }
          if (inside(refNormal, intervalRefNormal, v)) {
            vertices[vertexIndex++] = v;
          }

          prevV = v;
        }

        vertices.resize(vertexIndex);

        if (vertexIndex < 3)
          Warning << "Face cutting resulted with " << vertexIndex
                  << " vertices." << endl;

        vector pSurface = surface(vertices);
        scalar pArea = pSurface & patch_normal_;

        if (pArea > area[i])
          Warning << "Face cutting resulted with greater face."
                  << "ratio = " << (pArea / area[i]) << " > 1" << endl;

        if ((pSurface & Sf[i]) < 0)
          Warning << "Face cutting resulted opposite face normals." << endl;

        sum_total_pressure =
            sum_total_pressure + weighting_mult[i] * pArea * pt[i];
        sum_weighting = sum_weighting + weighting_mult[i] * pArea;
        sum_area = sum_area + pArea;
      }
    }

    reduce(sum_total_pressure, sumOp<scalar>());
    reduce(sum_area, sumOp<scalar>());
    reduce(sum_weighting, sumOp<scalar>());

    avg_total_pressure_by_angle_[n] = sum_total_pressure / sum_weighting;
    area_by_angle_[n] = sum_area;
    angle_[n] = n * angle_stride_;
  }

  return true;
}

bool intakeDistortion::write() {
  if (Pstream::master()) {
    Log << type() << " " << name() << " write:" << nl
        << "    writing distortion" << endl;

    auto &stream = file();

    scalarField distortion_((avg_total_pressure_ - avg_total_pressure_by_angle_) /
        avg_dynamic_pressure_);

    writeCommented(stream, "Timestep : ");
    writeCurrentTime(stream);
    stream << endl;
    writeCommented(stream, "Surface Vector : ");
    stream << patch_surface_ << endl;
    writeCommented(stream, "Surface Area : ");
    stream << patch_area_ << endl;
    writeCommented(stream, "Surface Normal : ");
    stream << patch_normal_ << endl;
    writeCommented(stream, "Surface Tangent : ");
    stream << patch_tangent_ << endl;
    writeCommented(stream, "Surface Binormal : ");
    stream << patch_binormal_ << endl;
    writeCommented(stream, "Surface Centroid : ");
    stream << patch_centroid_ << endl;
    writeCommented(stream, "Min. Total Pressure : ");
    stream << min(avg_total_pressure_by_angle_) << endl;
    writeCommented(stream, "Max. Total Pressure : ");
    stream << max(avg_total_pressure_by_angle_) << endl;
    writeCommented(stream, "Min. Distortion : ");
    stream << min(distortion_) << endl;
    writeCommented(stream, "Max. Distortion : ");
    stream << max(distortion_) << endl;
    writeCommented(stream, "Avg. Total Pressure : ");
    stream << avg_total_pressure_ << endl;
    writeCommented(stream, "Avg. Dynamic Pressure : ");
    stream << avg_dynamic_pressure_ << endl;

    writeTabbed(stream, "angle");
    writeTabbed(stream, "partialAvgTotalPressure");
    writeTabbed(stream, "distorsion");
    writeTabbed(stream, "area");
    stream << endl;

    for (label i = 0; i < angle_count_; ++i) {
      stream << angle_[i] << token::TAB << avg_total_pressure_by_angle_[i]
             << token::TAB << distortion_[i] << token::TAB << area_by_angle_[i]
             << token::TAB << endl;
    }

    Log << "    Min. Distortion :  " << min(distortion_) << endl
        << "    Max. Distortion :  " << max(distortion_) << endl;
  }

  return true;
}

} // namespace functionObjects
} // namespace Foam

// ************************************************************************* //
