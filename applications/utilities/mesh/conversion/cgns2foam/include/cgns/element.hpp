#pragma once

#include <cgnslib.h>
#include <cgnstypes.h>
#include <boost/iterator/zip_iterator.hpp>
#include <boost/iterator/counting_iterator.hpp>
#include <boost/iterator/transform_iterator.hpp>
#include <list>
#include <numeric>
#include <set>
#include <tuple>
#include <vector>
#include <string>
#include <algorithm>


#include "cgns/base.hpp"
#include "cgns/error.hpp"

namespace CGNS {

namespace impl {
/**
 * @brief Class for loading and managing CGNS element data
 * 
 * Handles reading element connectivity information from CGNS files,
 * classifying elements into face and cell categories, and providing
 * access to element data for mesh processing.
 */
class Elements : public virtual LoaderBase {

    /** @brief Indices of sections containing face elements */
    std::vector<int> _face_data;

    /** @brief Indices of sections containing cell elements */
    std::vector<int> _cell_data;

    /** @brief Total number of face elements */
    cgsize_t _provided_face_count;

    /** @brief Total number of cell elements */
    cgsize_t _provided_cell_count;

    bool _requires_internal_build = true;

    /**
     * @brief Classify element sections into face and cell categories
     * 
     * Reads all element sections from the CGNS file and classifies them
     * based on element type. Face elements include triangles, quadrilaterals,
     * and n-gons. Cell elements include tetrahedra, pyramids, prisms, 
     * hexahedra, and n-faces.
     * 
     * @throws std::runtime_error if element type cannot be classified
     */
    void load_data() {
        int nsections;

        CG_CHECK(cg_nsections,
                 (get_file(), get_base_index(), get_zone_index(), &nsections));

        BOOST_LOG_TRIVIAL(trace) << "Found " << nsections << " sections";

        std::vector<cgsize_t> face_starts;
        std::vector<int> face_data;
        std::vector<cgsize_t> cell_starts;
        std::vector<int> cell_data;

        auto add_node =
          [&](int index, ElementType_t type, cgsize_t start, cgsize_t end) {
              if (type == ElementType_t::TRI_3 ||
                  type == ElementType_t::QUAD_4 ||
                  type == ElementType_t::NGON_n) {
                  if (type == NGON_n) { _requires_internal_build = false; }
                  face_data.push_back(index);
                  face_starts.push_back(start);
                  _provided_face_count += end - start + 1;
              } else if (type == ElementType_t::TETRA_4 ||
                         type == ElementType_t::PYRA_5 ||
                         type == ElementType_t::PENTA_6 ||
                         type == ElementType_t::HEXA_8 ||
                         type == ElementType_t::NFACE_n) {
                  cell_data.push_back(index);
                  cell_starts.push_back(start);
                  _provided_cell_count += end - start + 1;
              } else {
                  BOOST_LOG_TRIVIAL(fatal)
                    << "Can't classify mixed zone "
                       "elements";

                  throw std::runtime_error("unhandled element type");
              }
          };

        for (int isec = 1; isec <= nsections; ++isec) {
            char secname[CGIO_MAX_NAME_LENGTH + 1];
            ElementType_t elem_type;
            cgsize_t start;
            cgsize_t end;
            int nbndry;
            int nparent;
            CG_CHECK(cg_section_read,
                     (get_file(),
                      get_base_index(),
                      get_zone_index(),
                      isec,
                      secname,
                      &elem_type,
                      &start,
                      &end,
                      &nbndry,
                      &nparent));

            cgsize_t conn_size;
            CG_CHECK(cg_ElementDataSize,
                     (get_file(),
                      get_base_index(),
                      get_zone_index(),
                      isec,
                      &conn_size));

            std::vector<cgsize_t> elem_conn(conn_size);

            if (elem_type == ElementType_t::MIXED) {
                std::vector<cgsize_t> offset;

                offset.resize(end - start + 2);

                CG_CHECK(cg_poly_elements_read,
                         (get_file(),
                          get_base_index(),
                          get_zone_index(),
                          isec,
                          elem_conn.data(),
                          offset.data(),
                          NULL));


                ElementType_t
                  first_type = static_cast<ElementType_t>(elem_conn[0]);
                add_node(isec, first_type, start, end);

            } else {
                add_node(isec, elem_type, start, end);
            }
        }

        std::vector<cgsize_t> face_sort_indices(face_starts.size());
        std::iota(face_sort_indices.begin(), face_sort_indices.end(), 0);
        std::sort(face_sort_indices.begin(),
                  face_sort_indices.end(),
                  [&](cgsize_t a, cgsize_t b) {
                      return face_starts[a] < face_starts[b];
                  });
        _face_data.resize(face_data.size());
        std::transform(face_sort_indices.begin(),
                       face_sort_indices.end(),
                       _face_data.begin(),
                       [&](cgsize_t index) { return face_data[index]; });


        std::vector<cgsize_t> cell_sort_indices(cell_starts.size());
        std::iota(cell_sort_indices.begin(), cell_sort_indices.end(), 0);
        std::sort(cell_sort_indices.begin(),
                  cell_sort_indices.end(),
                  [&](cgsize_t a, cgsize_t b) {
                      return cell_starts[a] < cell_starts[b];
                  });
        _cell_data.resize(cell_data.size());
        std::transform(cell_sort_indices.begin(),
                       cell_sort_indices.end(),
                       _cell_data.begin(),
                       [&](cgsize_t index) { return cell_data[index]; });
    }

protected:

    bool requires_internal_build() const { return _requires_internal_build; }

    /**
     * @brief Hash function for sets of identifiers
     * 
     * Provides a hash function for std::set containers used in
     * unordered containers for efficient face lookup operations.
     */
    struct id_set_hasher {
        /**
         * @brief Hash operator for sets
         * 
         * @tparam T Type of elements in the set
         * @param v Set to hash
         * @return Hash value
         */
        template <typename T>
        std::size_t operator()(std::set<T> const &v) const {
            std::size_t h = v.size();
            for (T const &x: v) {
                h ^= std::hash<T>()(x) + 0x9e3779b9 + (h << 6) + (h >> 2);
            }
            return h;
        }
    };

    /** @brief Get face element section indices @return Reference to face element indices */
    std::vector<int> const &get_face_elements_indices() const {
        return _face_data;
    }

    /** @brief Get cell element section indices @return Reference to cell element indices */
    std::vector<int> const &get_cell_elements_indices() const {
        return _cell_data;
    }

    /** @brief Get total cell count @return Number of cell elements */
    size_t get_provided_cell_count() const { return _provided_cell_count; }

    size_t get_provided_face_count() const { return _provided_face_count; }

    /**
     * @brief Structure containing element section data
     * 
     * Holds all information about an element section including
     * connectivity, offsets, and metadata.
     */
    struct element_data {
        /** @brief Element type */
        ElementType_t type;

        /** @brief Starting element index */
        cgsize_t start;

        /** @brief Ending element index */
        cgsize_t end;

        /** @brief Section name */
        std::string name;

        /** @brief Element connectivity data */
        std::vector<cgsize_t> connectivity;

        /** @brief Offset data for polygonal elements */
        std::vector<cgsize_t> offset;

        /** @brief Number of boundary elements */
        int nbndry;

        /** @brief Number of parent elements */
        int nparent;
    };

    /**
     * @brief Read element data from a specific section
     * 
     * Reads complete element information including connectivity data
     * for the specified element section.
     * 
     * @param ielements Section index to read
     * @return Complete element data structure
     * @throws std::runtime_error if element data cannot be read
     */
    element_data read_element_data(int ielements) const {
        char secname[CGIO_MAX_NAME_LENGTH + 1];
        ElementType_t elem_type;
        cgsize_t start;
        cgsize_t end;
        int nbndry;
        int nparent;

        CG_CHECK(cg_section_read,
                 (get_file(),
                  get_base_index(),
                  get_zone_index(),
                  ielements,
                  secname,
                  &elem_type,
                  &start,
                  &end,
                  &nbndry,
                  &nparent));

        cgsize_t conn_size;

        CG_CHECK(cg_ElementDataSize,
                 (get_file(),
                  get_base_index(),
                  get_zone_index(),
                  ielements,
                  &conn_size));

        std::vector<cgsize_t> connectivity(conn_size);
        std::vector<cgsize_t> offset;


        if (elem_type == ElementType_t::NGON_n ||
            elem_type == ElementType_t::NFACE_n ||
            elem_type == ElementType_t::MIXED) {
            offset.resize(end - start + 2);
            CG_CHECK(cg_poly_elements_read,
                     (get_file(),
                      get_base_index(),
                      get_zone_index(),
                      ielements,
                      connectivity.data(),
                      offset.data(),
                      NULL));
        } else {
            CG_CHECK(cg_elements_read,
                     (get_file(),
                      get_base_index(),
                      get_zone_index(),
                      ielements,
                      connectivity.data(),
                      NULL));
        }

        return element_data { elem_type,    start,  end,    secname,
                              connectivity, offset, nbndry, nparent };
    }

    /**
     * @brief Build connectivity lists from element data
     * 
     * Processes element data to extract individual element connectivity
     * and type information. Handles mixed, polygonal, and standard element types.
     * 
     * @param data Element data structure containing connectivity information
     * @return Tuple containing connectivity list and element types
     * @throws std::runtime_error if element processing fails
     */
    static std::tuple<connectivity_map, std::vector<ElementType_t>>
    build_connectivity(element_data const &data) {
        cgsize_t count = data.end - data.start + 1;

        connectivity_map result;
        std::vector<ElementType_t> types(count, data.type);

        result.resize(count);

        if (data.type == ElementType_t::MIXED) {
            size_t pos = 0;
            cgsize_t index = 0;
            while (pos < data.connectivity.size()) {
                int nnodes;
                ElementType_t
                  type = static_cast<ElementType_t>(data.connectivity[pos]);
                CG_CHECK(cg_npe,
                         (static_cast<ElementType_t>(data.connectivity[pos]),
                          &nnodes));


                result[index] = std::vector(data.connectivity.data() + pos + 1,
                                            data.connectivity.data() + pos +
                                              nnodes + 1);
                types[index++] = type;

                pos += nnodes + 1;
            }


        } else if (data.type == ElementType_t::NGON_n ||
                   data.type == ElementType_t::NFACE_n) {

            auto begin = boost::make_zip_iterator(
              boost::make_tuple(data.offset.begin(), ++data.offset.begin()));

            auto end = boost::make_zip_iterator(
              boost::make_tuple(--data.offset.end(), data.offset.end()));

            std::fill(types.begin(), types.end(), data.type);

            std::transform(begin, end, result.begin(), [&](auto zip) {
                auto [start, end] = zip;
                return std::vector<cgsize_t>(data.connectivity.data() + start,
                                             data.connectivity.data() + end);
            });

        } else {
            int nnodes;
            CG_CHECK(cg_npe, (static_cast<ElementType_t>(data.type), &nnodes));

            std::transform(boost::make_counting_iterator<cgsize_t>(0),
                           boost::make_counting_iterator<cgsize_t>(count),
                           result.begin(),
                           [&](cgsize_t x) {
                               return std::vector<cgsize_t>(
                                 data.connectivity.data() + x * nnodes,
                                 data.connectivity.data() + (x + 1) * nnodes);
                           });
        }

        return std::make_tuple(std::move(result), std::move(types));
    }

public:
    /**
     * @brief Constructor for elements loader
     * 
     * Initializes the elements loader and classifies all element sections
     * in the specified CGNS file into face and cell categories.
     * 
     * @param file_path Path to the CGNS file
     * @param base_name Name of the base to load (empty for first base)
     * @param zone_name Name of the zone to load (empty for first zone)
     */
    Elements(std::string const &file_path,
             std::string const &base_name,
             std::string const &zone_name)
        : LoaderBase(file_path, base_name, zone_name) {

        _provided_face_count = 0;
        _provided_cell_count = 0;

        load_data();
    }
};
}
}
