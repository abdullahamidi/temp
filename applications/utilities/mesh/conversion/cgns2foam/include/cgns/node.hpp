#pragma once


#include <cgns_io.h>
#include <cgnslib.h>
#include "cgns/base.hpp"

#include <algorithm>
#include <boost/iterator/iterator_categories.hpp>
#include <boost/iterator/iterator_traits.hpp>
#include <boost/iterator/transform_iterator.hpp>
#include <boost/iterator/zip_iterator.hpp>
#include <type_traits>
#include <vector>
#include <array>
#include <string>


#include "cgns/error.hpp"
#include "macro.hpp"

namespace CGNS {


namespace impl {
/**
 * @brief Class for loading and managing CGNS node coordinates
 * 
 * Handles reading coordinate data from CGNS files, supporting both
 * structured and unstructured zones in 2D and 3D configurations.
 * Supports single and double precision coordinate data.
 */
class Nodes : public virtual LoaderBase {

    /** @brief Storage for node coordinate data */
    node_list _data;

    /**
     * @brief Template method to read coordinates of specific data type
     * 
     * Reads coordinate data from the CGNS file for a specific coordinate
     * component (x, y, or z) and data type (single or double precision).
     * 
     * @tparam D Data type (RealSingle or RealDouble)
     * @param coordname Name of the coordinate component
     * @param index Index of the coordinate component (0=x, 1=y, 2=z)
     * @throws std::runtime_error if coordinate reading fails
     */
    template <DataType_t D>
    void read_coordinates(char const *const coordname, int const index) {
        using scalar_type = std::conditional_t<
          D == DataType_t::RealDouble,
          double,
          std::conditional_t<D == DataType_t::RealSingle, float, void>>;

        std::vector<scalar_type> result(_data.size());

        auto zone_type = get_zone_type();

        auto sizes = get_zone_sizes();

        if (is_3d() && zone_type == ZoneType_t::Structured) {
            cgsize_t rmin[] = { 1, 1, 1 };
            cgsize_t rmax[] = { sizes[0], sizes[1], sizes[2] };
            CG_CHECK(cg_coord_read,
                     (get_file(),
                      get_base_index(),
                      get_zone_index(),
                      coordname,
                      D,
                      rmin,
                      rmax,
                      result.data()));

        } else if (is_2d() && zone_type == ZoneType_t::Structured) {
            cgsize_t rmin[] = {
                1,
                1,
            };
            cgsize_t rmax[] = { sizes[0], sizes[1] };
            CG_CHECK(cg_coord_read,
                     (get_file(),
                      get_base_index(),
                      get_zone_index(),
                      coordname,
                      D,
                      rmin,
                      rmax,
                      result.data()));
        } else {
            cgsize_t rmin[] = { 1 };
            cgsize_t rmax[] = { sizes[0] };
            CG_CHECK(cg_coord_read,
                     (get_file(),
                      get_base_index(),
                      get_zone_index(),
                      coordname,
                      D,
                      rmin,
                      rmax,
                      result.data()));
        }


        auto begin = boost::make_zip_iterator(boost::make_tuple(result.begin(),
                                                                _data.begin()));
        auto end = boost::make_zip_iterator(boost::make_tuple(result.end(),
                                                              _data.end()));

        std::for_each(begin, end, [&](auto zip) {
            auto &[source, target] = zip;
            target[index] = source;
        });
    }

    /**
     * @brief Read all coordinate components from the CGNS file
     * 
     * Reads coordinate information for all available coordinate components
     * (typically x, y, z). Handles both single and double precision data
     * automatically based on the data type stored in the file.
     * 
     * @throws std::runtime_error if coordinate data cannot be read or data type is unsupported
     */
    void load_data() {

        _data.resize(get_point_count());

        BOOST_LOG_TRIVIAL(trace) << "Reading " << _data.size() << " points";

        int ncoords;
        CG_CHECK(cg_ncoords,
                 (get_file(), get_base_index(), get_zone_index(), &ncoords));

        for (int i = 1; i <= ncoords; ++i) {
            char coordname[CGIO_MAX_NAME_LENGTH];
            DataType_t datatype;

            CG_CHECK(cg_coord_info,
                     (get_file(),
                      get_base_index(),
                      get_zone_index(),
                      i,
                      &datatype,
                      coordname));

            if (datatype == DataType_t::RealSingle) {
                read_coordinates<DataType_t::RealSingle>(coordname, i - 1);
            } else if (datatype == DataType_t::RealDouble) {
                read_coordinates<DataType_t::RealDouble>(coordname, i - 1);
            } else {
                BOOST_LOG_TRIVIAL(fatal) << "Unknown coordinate data type";
                throw std::runtime_error("unhandled data type");
            }
        }


        BOOST_LOG_TRIVIAL(info) << "Read " << _data.size() << " points";
    }

protected:
    /**
     * @brief Get mutable reference to coordinate data
     * 
     * Provides access to the coordinate data for derived classes.
     * Triggers lazy loading if coordinates haven't been loaded yet.
     * 
     * @return Reference to the coordinate data vector
     */
    std::vector<std::array<double, 3>> &get_points() { return _data; }

public:
    /**
     * @brief Constructor for nodes loader
     * 
     * Initializes the nodes loader with the specified CGNS file and location.
     * 
     * @param file_path Path to the CGNS file
     * @param base_name Name of the base to load (empty for first base)
     * @param zone_name Name of the zone to load (empty for first zone)
     */
    Nodes(std::string const &file_path,
          std::string const &base_name,
          std::string const &zone_name)
        : LoaderBase(file_path, base_name, zone_name) {
        load_data();
    }

    /**
     * @brief Get read-only access to node coordinates
     * 
     * Provides const access to the loaded coordinate data.
     * Note: This method returns the internal data without triggering
     * lazy loading, so coordinates must be loaded separately.
     * 
     * @return Const reference to the coordinate data vector
     */
    std::vector<std::array<double, 3>> const &get_nodes() const {
        return _data;
    }
};
}
}
