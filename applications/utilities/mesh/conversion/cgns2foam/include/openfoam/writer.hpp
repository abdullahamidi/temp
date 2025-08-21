#pragma once

#include <string>
#include <vector>
#include <fstream>
#include <filesystem>
#include <iostream>
#include <format>
#include <boost/log/trivial.hpp>
#include "openfoam/source.hpp"

namespace openfoam {

namespace {
/**
 * @brief OpenFOAM file header template
 * 
 * Standard header format for OpenFOAM polyMesh files including
 * version information, format specification, and placeholders
 * for class and object names.
 */
constexpr char const *header =
  R"(/*--------------------------------*- C++ -*----------------------------------*\
| =========                |
| \\      /  F ield        | OpenFOAM: The Open Source CFD Toolbox 
|  \\    /   O peration    | Version: v2312
|   \\  /    A nd          | Website: www.openfoam.com
|    \\/     M anipulation |
\*---------------------------------------------------------------------------*/
FoamFile
{{
    version     2.0;
    format      ascii;
    class       {};
    arch        "LSB;label=32;scalar=64";
    location    "constant/polyMesh";
    object      {};
}}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
)";
}

/**
 * @brief OpenFOAM polyMesh writer class
 * 
 * Handles the conversion and writing of mesh data from CGNS format
 * to OpenFOAM polyMesh format. Creates all necessary files including
 * points, faces, owner, neighbour, and boundary definitions.
 */
class Writer {
    /** @brief Reference to the mesh data source */
    Source &_source;

    /** @brief Output directory path */
    std::string _output_dir;

    /** @brief polyMesh subdirectory path */
    std::string _polymesh_dir;

    /**
     * @brief Ensure directory exists, creating it if necessary
     * 
     * @param directory Path to the directory to create
     */
    static inline void ensure_directory_exists(std::string const &directory) {
        std::filesystem::create_directories(directory);
    }

    /**
     * @brief Format a 3D vector for OpenFOAM output
     * 
     * Converts a 3D coordinate array to OpenFOAM vector format "(x y z)".
     * 
     * @param vec 3D coordinate array
     * @return Formatted string representation
     */
    static inline std::string format_vector(std::array<double, 3> const &vec) {
        std::ostringstream ss;
        ss << "(" << vec[0] << " " << vec[1] << " " << vec[2] << ")";
        return ss.str();
    }

    /**
     * @brief Format a list of integers for OpenFOAM output
     * 
     * Converts a vector of integers to OpenFOAM list format "n(v1 v2 ... vn)".
     * 
     * @param lst Vector of integers to format
     * @return Formatted string representation
     */
    static inline std::string format_list(std::vector<long> const &lst) {
        std::ostringstream ss;
        ss << lst.size() << "(";
        for (size_t i = 0; i < lst.size(); ++i) {
            ss << lst[i];
            if (i + 1 != lst.size()) { ss << " "; }
        }
        ss << ")";
        return ss.str();
    }

    /**
     * @brief Create OpenFOAM file header
     * 
     * Generates a properly formatted OpenFOAM file header with the specified
     * class and object names.
     * 
     * @param class_name OpenFOAM class type for the file
     * @param object_name Object name for the file
     * @return Formatted header string
     */
    static inline std::string
    create_header(std::string const &class_name,
                  std::string const &object_name) {
        return std::vformat(header,
                            std::make_format_args(class_name, object_name));
    }


public:
    /**
     * @brief Constructor for OpenFOAM writer
     * 
     * Initializes the writer with a mesh data source and output directory.
     * Creates the necessary directory structure for polyMesh files.
     * 
     * @param src Reference to the mesh data source
     * @param out Output directory path
     */
    Writer(Source &src, std::string const &out)
        : _source(src)
        , _output_dir(out)
        , _polymesh_dir(out + "/constant/polyMesh") {
        ensure_directory_exists(_polymesh_dir);
        BOOST_LOG_TRIVIAL(info)
          << "Initialized OpenFOAMWriter with output directory: "
          << _output_dir;
    }

    /**
     * @brief Write points file
     * 
     * Writes the mesh node coordinates to the OpenFOAM points file.
     * 
     * @param points Vector of 3D coordinate points
     */
    void write_points(std::vector<std::array<double, 3>> const &points) {
        std::ofstream f(_polymesh_dir + "/points");
        f << create_header("vectorField", "Points");
        f << points.size() << "\n(\n";
        for (auto const &pt: points) { f << format_vector(pt) << "\n"; }
        f << ")\n";
        BOOST_LOG_TRIVIAL(info) << "Wrote " << points.size() << " points";
    }

    /**
     * @brief Write faces file
     * 
     * Writes the mesh face connectivity to the OpenFOAM faces file.
     * 
     * @param faces Vector of face node connectivity lists
     */
    void write_faces(std::vector<std::vector<long>> const &faces) {
        std::ofstream f(_polymesh_dir + "/faces");
        f << create_header("faceList", "Faces");
        f << faces.size() << "\n(\n";
        for (auto const &face: faces) { f << format_list(face) << "\n"; }
        f << ")\n";
        BOOST_LOG_TRIVIAL(info) << "Wrote " << faces.size() << " faces";
    }

    /**
     * @brief Write owner file
     * 
     * Writes the face owner information to the OpenFOAM owner file.
     * Each entry specifies which cell owns the corresponding face.
     * 
     * @param owners Vector of cell indices that own each face
     */
    void write_owner(std::vector<long> const &owners) {
        std::ofstream f(_polymesh_dir + "/owner");
        f << create_header("labelList", "Owner");
        f << owners.size() << "\n(\n";
        for (int o: owners) { f << o << "\n"; }
        f << ")\n";
        BOOST_LOG_TRIVIAL(info) << "Wrote " << owners.size() << " face owners";
    }

    /**
     * @brief Write neighbour file
     * 
     * Writes the face neighbour information to the OpenFOAM neighbour file.
     * Each entry specifies which cell is the neighbour of the corresponding face.
     * Only internal faces have neighbours.
     * 
     * @param neighbours Vector of cell indices that are neighbours of each internal face
     */
    void write_neighbour(std::vector<long> const &neighbours) {
        std::ofstream f(_polymesh_dir + "/neighbour");
        f << create_header("labelList", "Neighbour");
        f << neighbours.size() << "\n(\n";
        for (int n: neighbours) { f << n << "\n"; }
        f << ")\n";
        BOOST_LOG_TRIVIAL(info)
          << "Wrote " << neighbours.size() << " face neighbours";
    }

    /**
     * @brief Write boundary file
     * 
     * Writes the boundary patch information to the OpenFOAM boundary file.
     * Defines boundary patches with their types, face counts, and starting face indices.
     * 
     * @param boundaries Vector of boundary patch information
     */
    void write_boundary(std::vector<boundary_info> const &boundaries) {
        std::ofstream f(_polymesh_dir + "/boundary");
        f << create_header("polyBoundaryMesh", "Boundary");
        f << boundaries.size() << "\n(\n";
        for (auto const &b: boundaries) {
            f << " \"" << b.name << "\"\n {\n"
              << " type " << to_string(b.type) << ";\n"
              << " nFaces " << b.n_faces << ";\n"
              << " startFace " << b.start_face << ";\n"
              << " }\n";
        }
        f << ")\n";
        BOOST_LOG_TRIVIAL(info)
          << "Wrote " << boundaries.size() << " boundaries";
    }

    /**
     * @brief Write complete polyMesh
     * 
     * Writes all polyMesh files (points, faces, owner, neighbour, boundary)
     * by extracting data from the source and calling individual write methods.
     */
    void write() {
        write_points(_source.get_nodes());
        write_faces(_source.get_face_nodes());
        auto [owner, neighbour] = _source.get_owner_and_neighbor();
        write_owner(owner);
        write_neighbour(neighbour);
        write_boundary(_source.get_boundary_info());
    }
};
}
