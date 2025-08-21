#pragma once

#include <cgnslib.h>
#include <cgnstypes.h>
#include <algorithm>
#include <boost/log/trivial.hpp>
#include <execution>
#include <functional>
#include <iostream>
#include <iterator>
#include <limits>
#include <stdexcept>
#include <vector>
#include <array>
#include <string>
#include <unordered_map>

#include <boost/iterator/counting_iterator.hpp>
#include <boost/iterator/transform_iterator.hpp>
#include <boost/iterator/zip_iterator.hpp>
#include <boost/iterator/filter_iterator.hpp>


#include "base.hpp"
#include "cgns/boundary.hpp"
#include "node.hpp"
#include "face.hpp"
#include "element.hpp"

/** @brief Reserved number of faces per cell */
#define _RESERVED_CELL_FACE_COUNT 6

/** @brief Reserved number of nodes per face */
#define _RESERVED_FACE_NODE_COUNT 4

namespace CGNS {

namespace impl {
/**
 * @brief Cell face connectivity lookup table
 * 
 * Maps CGNS element types to their face connectivity patterns.
 * Each entry contains the local node indices that form each face
 * of the cell type. Used for building cell-face connectivity from
 * cell-node connectivity.
 */
std::unordered_map<ElementType_t, connectivity_map> const
  cgns_cell_face_schema =
    // clang-format off
    {
    {
            ElementType_t::TETRA_4,
            {   
                { 0, 2, 1 }, 
                { 0, 1, 3 }, 
                { 1, 2, 3 }, 
                { 2, 0, 3 } 
            }
        },
    { 
            ElementType_t::PYRA_5,
            { 
                { 0, 3, 2, 1 }, 
                { 0, 1, 4 }, 
                { 1, 2, 4 }, 
                { 2, 3, 4 }, 
                { 3, 0, 4 } 
            } 
        },
    { 
            ElementType_t::PENTA_6,
            { 
                { 0, 1, 4, 3 },
                { 1, 2, 5, 4 },
                { 2, 0, 3, 5 },
                { 0, 2, 1 },
                { 3, 4, 5 } 
            } 
        }, 
    {
            ElementType_t::HEXA_8, 
            {
                { 0, 3, 2, 1 },
                { 0, 1, 5, 4 }, 
                { 1, 2, 6, 5 }, 
                { 2, 3, 7, 6 },
                { 0, 4, 7, 3 }, 
                { 4, 5, 6, 7 }
            }
        }
};
// clang-format on

/**
 * @brief Class for loading and managing CGNS cell data
 * 
 * Handles cell connectivity, face relationships, and owner-neighbor
 * information. Provides functionality for building cell-face connectivity
 * from various input formats and enforcing proper mesh topology.
 */
class Cells
    : public virtual Elements
    , public virtual Faces {

    /** @brief Cell-face connectivity data with signed face indices */
    connectivity_map _data;

    /**
     * @brief Build face node lists for each cell
     * 
     * Uses the cell face connectivity table to determine which nodes
     * form each face of each cell. This is used when building cell-face
     * connectivity from cell-node connectivity.
     * 
     * @param cell_nodes List of node connectivity for each cell
     * @param cell_types Vector of element types for each cell
     * @return List of face node lists for each cell
     * @throws std::runtime_error if unknown cell type is encountered
     */
    static std::vector<std::vector<id_list>>
    build_cells_faces_nodes(connectivity_map const &cell_nodes,
                            std::vector<ElementType_t> const &cell_types) {

        BOOST_LOG_TRIVIAL(info) << "Generating internal faces";

        std::vector<std::vector<id_list>> result;
        result.reserve(cell_nodes.size());

        auto begin = boost::make_zip_iterator(boost::make_tuple(
          boost::make_counting_iterator<size_t>(0),
          cell_nodes.begin(),
          cell_types.begin()));
        auto end = boost::make_zip_iterator(boost::make_tuple(
          boost::make_counting_iterator<size_t>(cell_nodes.size()),
          cell_nodes.end(),
          cell_types.end()));

        std::for_each(begin, end, [&](auto elem) {
            auto [cid, nodes, type] = elem;
            auto it = cgns_cell_face_schema.find(type);
            if (it == cgns_cell_face_schema.end()) {
                throw std::runtime_error("Unknown cell type: " +
                                         std::to_string(type));
            }

            auto const &face_templates = it->second;

            std::vector<id_list> cell_faces;
            cell_faces.reserve(_RESERVED_CELL_FACE_COUNT);
            for (auto const &face_template: face_templates) {
                id_list face_nodes;
                face_nodes.reserve(_RESERVED_FACE_NODE_COUNT);
                for (cgsize_t local_idx: face_template) {
                    face_nodes.push_back(nodes[local_idx] - 1);
                }
                cell_faces.push_back(face_nodes);
            }
            result.push_back(cell_faces);
        });

        return result;
    }

    /**
     * @brief Build cell-face connectivity from available data
     * 
     * Determines the appropriate method for building cell-face connectivity
     * based on the available data in the CGNS file. Uses direct connectivity
     * if available, otherwise builds from cell-node connectivity.
     * 
     * @throws std::runtime_error if connectivity cannot be built
     */
    void load_data() {
        connectivity_map connectivity;
        std::vector<ElementType_t> types;

        for (int index: get_cell_elements_indices()) {
            auto data = read_element_data(index);
            auto [conn, ty] = Elements::build_connectivity(data);

            std::copy(conn.begin(),
                      conn.end(),
                      std::back_inserter(connectivity));

            std::copy(ty.begin(), ty.end(), std::back_inserter(types));

            BOOST_LOG_TRIVIAL(info)
              << "Read " << data.end - data.start + 1 << " cells from node "
              << data.name;
        }

        if (requires_internal_build()) {

            auto
              cells_faces_nodes = build_cells_faces_nodes(connectivity, types);

            std::
              unordered_map<id_set, std::pair<cgsize_t, id_list>, id_set_hasher>
                face_registry = Faces::generate(cells_faces_nodes);

            BOOST_LOG_TRIVIAL(info) << "Generating cells";

            // Build cell-face connectivity with proper signs
            for (auto const &cell_faces: cells_faces_nodes) {
                id_list cell_face_indices;
                cell_face_indices.reserve(_RESERVED_CELL_FACE_COUNT);

                for (auto const &face_nodes: cell_faces) {
                    id_set face_key = Faces::key_transform(face_nodes);

                    auto it = face_registry.find(face_key);

                    if (it != face_registry.end()) {
                        auto [face_idx, registered_face_nodes] = it->second;

                        if (Faces::check_orientation(face_nodes,
                                                     registered_face_nodes)) {
                            cell_face_indices.push_back(face_idx);
                        } else {
                            cell_face_indices.push_back(-face_idx);
                        }
                    } else {
                        throw std::runtime_error("Face not found in registry.");
                    }
                }
                _data.push_back(cell_face_indices);
            }
        } else {
            cgsize_t offset = get_provided_face_count() -
                              get_defined_boundary_face_count();

            id_list boundary_map;
            boundary_map.resize(get_defined_boundary_face_count());
            cgsize_t reorder_count = 0;

            for (auto &info: Boundaries::get_boundaries()) {
                if (info.ps_type == PointSetType_t::ElementList ||
                    info.ps_type == PointSetType_t::PointList) {

                    for (cgsize_t i = 0; i < (cgsize_t) info.points.size();
                         ++i) {
                        boundary_map[info.points[i]] = i;
                        reorder_count++;
                    }
                }
            }

            if (reorder_count > 0) {

                if (reorder_count != (cgsize_t) boundary_map.size()) {
                    throw std::runtime_error("boundary face reordering map "
                                             "does not conform with boundary "
                                             "face count");
                }

                std::transform(
                  connectivity.begin(),
                  connectivity.end(),
                  connectivity.begin(),
                  [&](id_list faces) {
                      std::transform(
                        faces.begin(),
                        faces.end(),
                        faces.begin(),
                        [&](cgsize_t fid) {
                            if (std::abs(fid) - 1 >= offset) {
                                return boundary_map[fid - offset - 1] + offset +
                                       1;
                            }

                            return fid;
                        });
                      return faces;
                  });
            }


            _data.insert(_data.end(), connectivity.begin(), connectivity.end());
        }
    }

    /**
     * @brief Build cell-face connectivity from owner-neighbor data
     * 
     * Constructs cell-face connectivity using face owner and neighbor
     * information. Positive indices indicate ownership, negative indices
     * indicate neighbor relationship.
     * 
     * @param owner Vector of cell indices that own each face
     * @param neighbor Vector of cell indices that neighbor each internal face
     */
    void
    build_cell_faces_from_owner_neighbor(id_list const &owner,
                                         id_list const &neighbor) {

        connectivity_map result(_data.size());

        // Owner assignments: +1-based face index
        for (cgsize_t i = 0; i < (cgsize_t) owner.size(); ++i) {
            result[owner[i]].push_back(i + 1);
        }

        // Neighbor assignments: -1-based face index
        for (cgsize_t i = 0; i < (cgsize_t) neighbor.size(); ++i) {
            result[neighbor[i]].push_back(-i - 1);
        }

        _data = { result.begin(), result.end() };
    }

protected:
    /**
     * @brief Check if face orientation information is available
     * 
     * Verifies that the cell-face connectivity contains both positive
     * and negative face indices, indicating proper face orientation.
     * 
     * @return true if face flip map is present, false otherwise
     */
    bool has_face_flip_map() const {
        int fmin = INT_MAX;
        int fmax = INT_MIN;
        for (auto const &faces: _data) {
            for (int face: faces) {
                fmin = std::min(fmin, face);
                fmax = std::max(fmax, face);
            }
        }
        return (fmin < 0 && fmax > 0);
    }

    /**
     * @brief Get mutable cell-face connectivity
     * 
     * Provides access to cell-face connectivity data for derived classes.
     * Triggers lazy loading if data hasn't been loaded yet.
     * 
     * @return Reference to cell-face connectivity data
     */
    connectivity_map &get_cell_faces() { return _data; }

    cgsize_t internal_face_count() const {
        auto const [owner, neighbor] = get_owner_and_neighbor();
        return neighbor.size();
    }

    cgsize_t
    count_invalid_upper_triangular_faces(id_list const &owner,
                                         id_list const &neighbor) {

        cgsize_t internal_face_count = neighbor.size();
        cgsize_t result = 0;
        cgsize_t max_int = std::numeric_limits<cgsize_t>::max();

        std::for_each(
          boost::make_zip_iterator(boost::make_tuple(
            boost::make_counting_iterator<cgsize_t>(0),
            _data.begin())),
          boost::make_zip_iterator(boost::make_tuple(
            boost::make_counting_iterator<cgsize_t>(_data.size()),
            _data.end())),
          [&](auto const &zip) {
              auto [ci, faces] = zip;

              id_list cells(faces.size(), -1);

              // Masks for vectors
              for (size_t i = 0; i < faces.size(); ++i) {
                  if (faces[i] > 0 && faces[i] <= internal_face_count) {
                      cells[i] = neighbor[faces[i] - 1];
                  } else if (faces[i] < 0) {
                      cells[i] = owner[-faces[i] - 1];
                  }
              }

              // Non-owned faces to max_int
              for (size_t i = 0; i < cells.size(); ++i) {
                  if (ci > cells[i]) { cells[i] = max_int; }
              }

              // Order of sorted cells
              id_list order(cells.size());
              std::iota(order.begin(), order.end(), 0);
              std::sort(order.begin(),
                        order.end(),
                        [&cells](size_t a, size_t b) {
                            return cells[a] < cells[b];
                        });

              std::vector<bool> valid(cells.size(), false);
              id_list face_sorted(faces.size());
              for (size_t i = 0; i < faces.size(); ++i) {
                  face_sorted[i] = faces[order[i]];
                  valid[i] = (cells[order[i]] != max_int);
              }

              // Shifted arrays for comparison
              id_list face_shifted(faces.size());
              std::vector<bool> valid_shift(faces.size());
              face_shifted[0] = face_sorted.back();
              valid_shift[0] = valid.back();
              for (size_t i = 1; i < faces.size(); ++i) {
                  face_shifted[i] = face_sorted[i - 1];
                  valid_shift[i] = valid[i - 1];
              }

              // Compare faces
              std::vector<bool> mask(faces.size(), false);
              for (size_t i = 1; i < faces.size(); ++i) {// skip the first
                  if (face_sorted[i] < face_shifted[i] && valid[i] &&
                      valid_shift[i]) {
                      mask[i] = true;
                  }
              }
              result += std::count(mask.begin(), mask.end(), true);
          });

        return result;
    }

    bool check_boundary_defination(id_list &owner, id_list &neighbor) {
        cgsize_t total_face_count = get_face_nodes().size();
        cgsize_t
          defined_boundary_face_count = get_defined_boundary_face_count();
        cgsize_t
          internal_face_count = total_face_count - defined_boundary_face_count;

        return (((cgsize_t) owner.size()) == total_face_count &&
                ((cgsize_t) neighbor.size()) == internal_face_count);
    }

    void move_undefined_boundaries(id_list &owner, id_list &neighbor) {
        // TODO:
    }

    void enforce_owner_neighbor_relation(id_list &owner, id_list &neighbor) {
        if (! check_boundary_defination(owner, neighbor)) {
            BOOST_LOG_TRIVIAL(warning)
              << "Unconforming face counts possibly due to "
                 "missing boundary face data";
        }

        using zip_type = boost::
          tuple<cgsize_t, id_list &, cgsize_t &, cgsize_t &>;
        cgsize_t minsize = std::min(owner.size(), neighbor.size());

        auto begin = boost::make_zip_iterator(boost::make_tuple(
          boost::make_counting_iterator<cgsize_t>(0),
          get_face_nodes().begin(),
          owner.begin(),
          neighbor.begin()));

        auto end = boost::make_zip_iterator(boost::make_tuple(
          boost::make_counting_iterator<cgsize_t>(minsize),
          get_face_nodes().begin() + minsize,
          owner.begin() + minsize,
          neighbor.begin() + minsize));

        auto filter = [](zip_type const &zip) {
            auto const &[fid, nodes, oid, nid] = zip;
            return nid < oid;
        };

        cgsize_t count = 0;
        auto new_cell_faces = std::vector(_data.begin(), _data.end());

        std::for_each(boost::make_filter_iterator(filter, begin, end),
                      boost::make_filter_iterator(filter, end, end),
                      [&](zip_type zip) {
                          auto &[fid, nodes, oid, nid] = zip;

                          // Flip face orientation
                          std::reverse(nodes.begin(), nodes.end());
                          // Find the local indices in cell_faces_ arrays and flip sign
                          for (cgsize_t &cf: new_cell_faces[nid]) {
                              if (cf == -fid - 1) { cf *= -1; }
                          }
                          for (cgsize_t &cf: new_cell_faces[oid]) {
                              if (cf == fid + 1) { cf *= -1; }
                          }
                          std::swap(oid, nid);
                          count++;
                      });

        if (count > 0) {
            BOOST_LOG_TRIVIAL(info)
              << "Corrected " << count << " owner-neighbor inconsistency";
        }
    }

    void
    enforce_upper_triangular_face_ordering(id_list &owner, id_list &neighbor) {

        if (! check_boundary_defination(owner, neighbor)) {
            BOOST_LOG_TRIVIAL(warning)
              << "Owner - neighbor map doesn't conform "
                 "with boundary face count";
        }

        cgsize_t total_face_count = get_face_nodes().size();
        cgsize_t minsize = (cgsize_t) std::min(owner.size(), neighbor.size());
        id_list new_face_order(total_face_count);
        std::iota(new_face_order.begin(), new_face_order.end(), 0);
        cgsize_t last_new_face = -1;

        auto begin = boost::make_zip_iterator(boost::make_tuple(
          boost::make_counting_iterator<cgsize_t>(0),
          _data.begin()));
        auto end = boost::make_zip_iterator(boost::make_tuple(
          boost::make_counting_iterator<cgsize_t>(_data.size()),
          _data.end()));

        std::for_each(begin, end, [&](auto zip) {
            auto [ci, faces] = zip;

            id_list internal_face_ids;

            std::copy_if(faces.begin(),
                         faces.end(),
                         std::back_insert_iterator { internal_face_ids },
                         [&](auto id) { return std::abs(id) <= minsize; });

            std::transform(internal_face_ids.begin(),
                           internal_face_ids.end(),
                           internal_face_ids.begin(),
                           [](cgsize_t id) { return std::abs(id) - 1; });

            id_list neighbor_face_ids;
            std::copy_if(internal_face_ids.begin(),
                         internal_face_ids.end(),
                         std::back_insert_iterator { neighbor_face_ids },
                         [&](auto id) { return owner[id] == ci; });

            id_list neigbor_ids(neighbor_face_ids.size());
            std::transform(neighbor_face_ids.begin(),
                           neighbor_face_ids.end(),
                           neigbor_ids.begin(),
                           [&](cgsize_t fid) { return neighbor[fid]; });

            std::vector<size_t> order(neighbor_face_ids.size());
            std::iota(order.begin(), order.end(), 0);
            std::sort(order.begin(), order.end(), [&](size_t i, size_t j) {
                return neigbor_ids[i] < neigbor_ids[j];
            });

            std::for_each(order.begin(), order.end(), [&](size_t index) {
                new_face_order[++last_new_face] = neighbor_face_ids[index];
            });
        });

        if (last_new_face != minsize - 1) {
            BOOST_LOG_TRIVIAL(warning)
              << "Incomplete number of faces produced! Please check if "
                 "boundary data is correct and complete in provided mesh";
        }

        reorder_faces(new_face_order);

        // Update owner/neigh maps to the new order
        id_list new_owner(owner.size());
        id_list new_neighbor(neighbor.size());
        for (cgsize_t fi = 0; fi < total_face_count; ++fi) {
            new_owner[fi] = owner[new_face_order[fi]];
        }
        for (cgsize_t fi = 0; fi < minsize; ++fi) {
            new_neighbor[fi] = neighbor[new_face_order[fi]];
        }

        build_cell_faces_from_owner_neighbor(new_owner, new_neighbor);
    }

public:
    /**
     * @brief Constructor for cells loader
     * 
     * Initializes the cells loader by inheriting from multiple base classes
     * to provide access to loader, element, node, and face functionality.
     * 
     * @param file_path Path to the CGNS file
     * @param base_name Name of the base to load (empty for first base)
     * @param zone_name Name of the zone to load (empty for first zone)
     */
    Cells(std::string const &file_path,
          std::string const &base_name,
          std::string const &zone_name)
        : LoaderBase(file_path, base_name, zone_name)
        , Elements(file_path, base_name, zone_name)
        , Nodes(file_path, base_name, zone_name)
        , Boundaries(file_path, base_name, zone_name)
        , Faces(file_path, base_name, zone_name) {
        load_data();

        if (! has_face_flip_map()) {
            BOOST_LOG_TRIVIAL(info) << "File does not contain face flip map.";
            throw std::runtime_error("No face flip map!");
        }
    }

    /**
     * @brief Get read-only access to cell-face connectivity
     * 
     * Provides const access to the cell-face connectivity data.
     * Note: This method returns the internal data without triggering
     * lazy loading, so connectivity must be loaded separately.
     * 
     * @return Const reference to cell-face connectivity data
     */
    connectivity_map const &get_cell_faces() const { return _data; }

    /**
     * @brief Extract owner and neighbor information from cell-face connectivity
     * 
     * Analyzes the signed face indices in cell-face connectivity to determine
     * which cell owns each face and which cell is the neighbor for internal faces.
     * Positive face indices indicate ownership, negative indices indicate neighbor.
     * 
     * @return Tuple containing owner and neighbor vectors
     */
    std::tuple<id_list, id_list> get_owner_and_neighbor() const {

        // Find max face index (positive and negative!)
        cgsize_t positive_max = 0;
        cgsize_t negative_min = 0;
        for (auto const &faces: _data) {
            for (cgsize_t fi: faces) {
                positive_max = std::max(positive_max, fi);
                negative_min = std::min(negative_min, fi);
            }
        }

        id_list owner(positive_max, std::numeric_limits<cgsize_t>::min());
        id_list neighbor(-negative_min, std::numeric_limits<cgsize_t>::min());

        std::for_each(
          std::execution::par,
          boost::make_zip_iterator(boost::make_tuple(
            boost::make_counting_iterator(size_t(0)),
            _data.cbegin())),
          boost::make_zip_iterator(boost::make_tuple(
            boost::make_counting_iterator(_data.size()),
            _data.cend())),
          [&](boost::tuple<size_t const &, id_list const &> faces) {
              auto [cid, fl] = faces;
              for (cgsize_t face: fl) {

                  if (face > 0) {
                      owner[face - 1] = cid;
                  } else if (face < 0) {
                      neighbor[-face - 1] = cid;
                  }
              }
          });

        bool invalid_owner = std::transform_reduce(
          std::execution::par,
          owner.begin(),
          owner.end(),
          false,
          std::logical_or {},
          [](cgsize_t oid) { return oid < 0; });

        if (invalid_owner) {
            BOOST_LOG_TRIVIAL(warning)
              << "Invalid owner cell indices found possibly due to missing "
                 "boundary data.";
        }

        cgsize_t index = -1;
        bool invalid_neighbor = std::transform_reduce(
          // std::execution::par,
          neighbor.begin(),
          neighbor.end(),
          false,
          std::logical_or {},
          [&](cgsize_t nid) {
              ++index;
              if (nid < 0) {
                  std::cout
                    << "index: " << index << " owner " << owner[index]
                    << std::endl;
              }
              return nid < 0;
          });

        if (invalid_neighbor) {
            BOOST_LOG_TRIVIAL(warning)
              << "Invalid neighbor cell indices found possibly due to missing "
                 "boundary data.";
        }


        return std::make_tuple(std::move(owner), std::move(neighbor));
    }
};
}
}
