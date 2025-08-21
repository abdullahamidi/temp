#pragma once

#include <vector>
#include <string>

namespace openfoam {

enum BOUNDARY_TYPE {
    wall,
    patch,
    symmetry
};

static constexpr inline std::string to_string(BOUNDARY_TYPE type) {
    switch (type) {
        case wall: return "wall";
        case patch: return "patch";
        case symmetry: return "symmetry";
    }

    return "patch";
}

struct boundary_info {
    std::string name;
    BOUNDARY_TYPE type;
    size_t start_face;
    size_t n_faces;
};

class Source {
public:
    virtual ~Source() noexcept(false) = default;

    //! Return the mesh points (each point is a 3D coordinate).
    virtual std::vector<std::array<double, 3>> get_nodes() = 0;

    //! Return the face definitions (each face is a list of point indices).
    virtual std::vector<std::vector<long>> get_face_nodes() = 0;

    //! Return owner and neighbor info (format as you implement).
    virtual std::tuple<std::vector<long>, std::vector<long>>
    get_owner_and_neighbor() = 0;

    //! Return info about mesh boundaries.
    virtual std::vector<openfoam::boundary_info> get_boundary_info() = 0;
};
}
