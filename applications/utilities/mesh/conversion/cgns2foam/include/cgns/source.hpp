#pragma once
#include <cgnslib.h>
#include <algorithm>
#include <limits>
#include "base.hpp"
#include "boundary.hpp"
#include "cell.hpp"
#include "element.hpp"
#include "face.hpp"
#include "node.hpp"
#include "openfoam/source.hpp"

#include <stdexcept>
#include <vector>
#include <string>

namespace CGNS {

class OpenfoamSource final
    : public openfoam::Source
    , public virtual impl::LoaderBase
    , public virtual impl::Nodes
    , public virtual impl::Elements
    , public virtual impl::Boundaries
    , public virtual impl::Faces
    , public virtual impl::Cells {

public:
    OpenfoamSource(std::string const &file_path,
                   std::string const &base_name = {},
                   std::string const &zone_name = {})
        : LoaderBase(file_path, base_name, zone_name)
        , Nodes(file_path, base_name, zone_name)
        , Elements(file_path, base_name, zone_name)
        , Boundaries(file_path, base_name, zone_name)
        , Faces(file_path, base_name, zone_name)
        , Cells(file_path, base_name, zone_name) {}

    void load() {

        auto [owner, neighbor] = get_owner_and_neighbor();

        enforce_owner_neighbor_relation(owner, neighbor);

        cgsize_t invalid_face_count =
          count_invalid_upper_triangular_faces(owner, neighbor);

        if (invalid_face_count > 0) {
            enforce_upper_triangular_face_ordering(owner, neighbor);

            BOOST_LOG_TRIVIAL(info)
              << "Corrected " << invalid_face_count
              << " faces to be in upper triangular order";
        }


        BOOST_LOG_TRIVIAL(info)
          << "Total " << get_face_nodes().size() << " faces";
        BOOST_LOG_TRIVIAL(info)
          << "Total " << get_cell_faces().size() << " cells";
    }

    std::vector<std::array<double, 3>> get_nodes() override {
        return Nodes::get_nodes();
    }

    //! Return the face definitions (each face is a list of point indices).
    std::vector<std::vector<long>> get_face_nodes() override {
        return Faces::get_face_nodes();
    }

    //! Return owner and neighbor info (format as you implement).
    std::tuple<std::vector<long>, std::vector<long>>
    get_owner_and_neighbor() override {
        return Cells::get_owner_and_neighbor();
    }

    //! Return info about mesh boundaries.
    virtual std::vector<openfoam::boundary_info> get_boundary_info() override {

        auto bd = Boundaries::get_boundaries();

        std::vector<openfoam::boundary_info> result(bd.size());

        cgsize_t offset = internal_face_count();

        std::transform(bd.cbegin(), bd.cend(), result.begin(), [&](auto &info) {
            openfoam::BOUNDARY_TYPE type = openfoam::patch;

            switch (info.type) {
                case BCType_t::BCWall: type = openfoam::wall; break;
                case BCType_t::BCSymmetryPlane:
                    type = openfoam::symmetry;
                    break;
                default: type = openfoam::patch; break;
            }

            if (info.remap.first != info.remap.second) {
                return openfoam::boundary_info {
                    .name = info.name,
                    .type = type,
                    .start_face = size_t(offset + info.remap.first),
                    .n_faces = size_t(info.remap.second - info.remap.first)
                };
            }

            throw std::runtime_error("Unable to handle boundary "
                                     "defination");
        });

        return result;
    }
};
}
