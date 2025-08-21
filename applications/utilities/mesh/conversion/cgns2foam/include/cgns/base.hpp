#pragma once

#include <array>
#include <set>
#include <stdexcept>
#include <string>
#include <vector>
#include "error.hpp"
#include "macro.hpp"

#include <boost/log/trivial.hpp>


#include <cgnslib.h>
#include <cgns_io.h>
#include <cgnstypes.h>

namespace CGNS {

using coord = std::array<double, 3>;

/**
 * @brief 3D coordinate representation using double precision
 * 
 * Represents a single node coordinate in 3D space with x, y, z components.
 */
using node_coord = coord;

/**
 * @brief Collection of 3D node coordinates
 * 
 * Vector containing multiple node coordinates for mesh representation.
 */
using node_list = std::vector<node_coord>;

/**
 * @brief List of CGNS size type identifiers
 * 
 * Vector of cgsize_t values used for indexing and connectivity information.
 */
using id_list = std::vector<cgsize_t>;

/**
 * @brief Set of unique CGNS identifiers
 * 
 * Set container for storing unique cgsize_t values, typically used for
 * efficient lookup and duplicate elimination.
 */
using id_set = std::set<cgsize_t>;

/**
 * @brief Pair representing a connected element with its connectivity list
 * 
 * First element is the identifier, second is the list of connected nodes/elements.
 */
using connected_pair = std::pair<cgsize_t, id_list>;

/**
 * @brief Connectivity mapping structure
 * 
 * Vector of id_list elements representing connectivity relationships
 * between mesh entities (nodes, faces, cells).
 */
using connectivity_map = std::vector<id_list>;

namespace impl {
/**
 * @brief Base class for loading CGNS files
 * 
 * Provides fundamental functionality for opening CGNS files and accessing
 * base and zone information. This class serves as the foundation for
 * specialized loaders that handle specific mesh data types.
 */
class LoaderBase {

    /** @brief Path to the CGNS file */
    std::string _file_path;

    /** @brief CGNS file handle */
    int _file;

    /** @brief Index of the selected base in the CGNS file */
    int _base_index;

    /** @brief Index of the selected zone in the CGNS file */
    int _zone_index;

    /** @brief Cell dimension (1D, 2D, or 3D) */
    int _cell_dim;

    /** @brief Physical dimension of the coordinate space */
    int _phys_dim;

    /** @brief Zone size information array */
    std::array<cgsize_t, CGIO_MAX_DIMENSIONS> _zone_sizes;

    /** @brief Type of the zone (Structured or Unstructured) */
    ZoneType_t _zone_type;

    /**
     * @brief Load and open the CGNS file
     * 
     * Opens the CGNS file in read mode and validates the file handle.
     * 
     * @throws std::runtime_error if file cannot be opened or is invalid
     */
    void load_file() {
        BOOST_LOG_TRIVIAL(info) << "Loading CGNS file: " << _file_path;

        CG_CHECK(cg_open, (_file_path.c_str(), CG_MODE_READ, &_file));

        DYN_ASSERT(_file != 0)
    }

    /**
     * @brief Find and select a base in the CGNS file
     * 
     * Searches for a base with the specified name. If base_name is empty,
     * selects the first available base. Reads base information including
     * cell and physical dimensions.
     * 
     * @param base_name Name of the base to find (empty for first base)
     * @throws std::runtime_error if no suitable base is found
     */
    void find_base(std::string const &base_name) {
        int nbases;

        CG_CHECK(cg_nbases, (_file, &nbases));

        for (cgsize_t ibase = 1; ibase <= nbases; ++ibase) {
            char basename[CGIO_MAX_NAME_LENGTH + 1];

            CG_CHECK(cg_base_read,
                     (_file, ibase, basename, &_cell_dim, &_phys_dim));

            if (base_name.empty() || base_name == basename) {
                _base_index = ibase;

                BOOST_LOG_TRIVIAL(info)
                  << "Base " << ibase << ": " << basename
                  << ", Cell Dim: " << _cell_dim << ", Phys Dim: " << _phys_dim;

                return;
            }
        }

        DYN_ASSERT(_base_index > 0)
    }

    /**
     * @brief Find and select a zone in the current base
     * 
     * Searches for a zone with the specified name. If zone_name is empty,
     * selects the first available zone. Reads zone information including
     * type and size data.
     * 
     * @param zone_name Name of the zone to find (empty for first zone)
     * @throws std::runtime_error if no suitable zone is found
     */
    void find_zone(std::string const &zone_name) {
        int nzones;

        CG_CHECK(cg_nzones, (_file, _base_index, &nzones));

        for (cgsize_t izone = 1; izone <= nzones; ++izone) {
            char zonename[CGIO_MAX_NAME_LENGTH];

            CG_CHECK(cg_zone_type, (_file, _base_index, izone, &_zone_type));

            CG_CHECK(cg_zone_read,
                     (_file, _base_index, izone, zonename, _zone_sizes.data()));

            if (zone_name.empty() || zone_name == zonename) {
                _zone_index = izone;

                BOOST_LOG_TRIVIAL(info) << "Zone " << izone << ": " << zonename;
            }
        }

        DYN_ASSERT(_zone_index > 0)
    }


protected:
    /** @brief Get the CGNS file handle @return File handle */
    int get_file() const { return _file; }

    /** @brief Get the base index @return Base index */
    int get_base_index() const { return _base_index; }

    /** @brief Get the zone index @return Zone index */
    int get_zone_index() const { return _zone_index; }

    /** @brief Get the cell dimension @return Cell dimension (1, 2, or 3) */
    int get_cell_dim() const { return _cell_dim; }

    /** @brief Get the physical dimension @return Physical dimension (2 or 3) */
    int get_phys_dim() const { return _phys_dim; }

    /** @brief Get the zone sizes array @return Reference to zone sizes */
    std::array<cgsize_t, CGIO_MAX_DIMENSIONS> const &get_zone_sizes() const {
        return _zone_sizes;
    }

    /** @brief Get the zone type @return Zone type (Structured or Unstructured) */
    ZoneType_t get_zone_type() const { return _zone_type; }


public:
    /**
     * @brief Constructor for CGNS loader base
     * 
     * Initializes the loader by opening the specified CGNS file and
     * selecting the appropriate base and zone.
     * 
     * @param file_path Path to the CGNS file
     * @param base_name Name of the base to load (empty for first base)
     * @param zone_name Name of the zone to load (empty for first zone)
     * @throws std::runtime_error if file cannot be loaded or base/zone not found
     */
    LoaderBase(std::string const &file_path,
               std::string const &base_name,
               std::string const &zone_name)
        : _file_path(file_path)
        , _file(0)
        , _base_index(-1)
        , _zone_index(-1)
        , _zone_sizes({}) {
        load_file();
        find_base(base_name);
        find_zone(zone_name);
    }

    /** @brief Copy constructor (deleted) */
    LoaderBase(LoaderBase const &) = delete;

    /** @brief Move constructor (default) */
    LoaderBase(LoaderBase &&) noexcept = default;

    /** @brief Copy assignment operator (deleted) */
    LoaderBase &operator=(LoaderBase const &) = delete;

    /** @brief Move assignment operator (default) */
    LoaderBase &operator=(LoaderBase &&) noexcept = default;

    /**
     * @brief Virtual destructor
     * 
     * Closes the CGNS file if it's still open.
     * 
     * @throws std::runtime_error if file cannot be closed properly
     */
    virtual ~LoaderBase() noexcept(false) {
        if (_file >= 0) { CG_CHECK(cg_close, (_file)); }
    }

    /** @brief Check if mesh contains 3D volume cells @return true if cell dimension is 3 */
    bool containing_volume_cells() const { return _cell_dim == 3; }

    /** @brief Check if mesh contains 2D surface cells @return true if cell dimension is 2 */
    bool containing_surface_cells() const { return _cell_dim == 2; }

    /** @brief Check if mesh contains 1D line cells @return true if cell dimension is 1 */
    bool containing_line_cells() const { return _cell_dim == 1; }

    /** @brief Check if coordinate space is 3D @return true if physical dimension is 3 */
    bool is_3d() const { return _phys_dim == 3; }

    /** @brief Check if coordinate space is 2D @return true if physical dimension is 2 */
    bool is_2d() const { return _phys_dim == 2; }

    /**
     * @brief Get the total number of points in the zone
     * 
     * Calculates the point count based on zone type and dimensions.
     * For structured zones, multiplies the dimensional sizes.
     * For unstructured zones, uses the first zone size value.
     * 
     * @return Total number of points
     * @throws std::runtime_error if zone type is unhandled
     */
    cgsize_t get_point_count() const {
        if (is_3d()) {
            if (_zone_type == ZoneType_t::Structured) {
                return _zone_sizes[0] * _zone_sizes[1] * _zone_sizes[2];
            }

            if (_zone_type == ZoneType_t::Unstructured) {
                return _zone_sizes[0];
            }
        } else if (is_2d()) {

            if (_zone_type == ZoneType_t::Structured) {
                return _zone_sizes[0] * _zone_sizes[1];
            }

            if (_zone_type == ZoneType_t::Unstructured) {
                return _zone_sizes[0];
            }
        }

        BOOST_LOG_TRIVIAL(fatal) << "Unable to determine node count";
        throw std::runtime_error("unhandled zone type");
    }

    /**
     * @brief Get the total number of elements in the zone
     * 
     * Calculates the element count based on zone type and dimensions.
     * For 3D structured zones, multiplies the cell dimensional sizes.
     * For 3D unstructured zones, uses the second zone size value.
     * 
     * @return Total number of elements
     * @throws std::runtime_error if element count cannot be determined
     */
    cgsize_t get_element_count() const {
        if (is_3d()) {
            if (_zone_type == ZoneType_t::Structured) {
                return _zone_sizes[3] * _zone_sizes[4] * _zone_sizes[5];
            }

            if (_zone_type == ZoneType_t::Unstructured) {
                return _zone_sizes[1];
            }
        } else if (is_2d()) {

            if (_zone_type == ZoneType_t::Structured) {
                return _zone_sizes[3] * _zone_sizes[4];
            }

            if (_zone_type == ZoneType_t::Unstructured) {
                return _zone_sizes[1];
            }
        }

        BOOST_LOG_TRIVIAL(fatal) << "Unable to determine element count";
        throw std::runtime_error("unknown element count");
    }
};

}

}
