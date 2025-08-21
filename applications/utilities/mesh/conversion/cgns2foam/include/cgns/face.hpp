#pragma once
#include <cgnslib.h>
#include <unordered_set>
#include <vector>
#include <stdexcept>
#include <algorithm>
#include <iterator>

#include <boost/iterator/zip_iterator.hpp>
#include <boost/iterator/counting_iterator.hpp>
#include <boost/iterator/transform_iterator.hpp>

#include "cgns/base.hpp"
#include "cgns/boundary.hpp"
#include "cgns/element.hpp"
#include "cgns/node.hpp"

namespace CGNS {

namespace impl {
class Faces
    : public virtual Nodes
    , public virtual Boundaries
    , public virtual Elements {

    connectivity_map _data;

    cgsize_t _built_face_count = 0;

    void load_data() {
        std::list<id_list> face_nodes;

        for (int index: get_face_elements_indices()) {
            auto data = read_element_data(index);

            auto [connectivity, _] = Elements::build_connectivity(data);
            std::for_each(connectivity.begin(),
                          connectivity.end(),
                          [](auto &nodes) {
                              std::transform(nodes.begin(),
                                             nodes.end(),
                                             nodes.begin(),
                                             [](cgsize_t nid) {
                                                 return nid - 1;
                                             });
                          });

            std::copy(connectivity.begin(),
                      connectivity.end(),
                      std::back_inserter(face_nodes));

            BOOST_LOG_TRIVIAL(info)
              << "Read " << data.end - data.start + 1 << " faces from node "
              << data.name;
        }


        _data = connectivity_map(face_nodes.begin(), face_nodes.end());

        cgsize_t reordered_face_count = 0;
        id_list boundary_map;
        boundary_map.reserve(get_defined_boundary_face_count());

        for (auto &info: Boundaries::get_boundaries()) {
            if (info.ps_type == PointSetType_t::ElementList ||
                info.ps_type == PointSetType_t::PointList) {
                std::copy(info.points.begin(),
                          info.points.end(),
                          std::back_inserter(boundary_map));

                info.remap = {
                    reordered_face_count,
                    reordered_face_count + (cgsize_t) info.points.size()
                };
                reordered_face_count += info.points.size();
            } else if (info.ps_type == PointSetType_t::ElementRange ||
                       info.ps_type == PointSetType_t::PointRange) {
                info.remap = { info.points[0], info.points[1] + 1 };
            }
        }

        if (! boundary_map.empty()) {
            reorder_faces(boundary_map,
                          get_provided_face_count() -
                            get_defined_boundary_face_count());
        }
    }


protected:

    static id_set key_transform(id_list const &nodes) {

        return { nodes.begin(), nodes.end() };
    }

    std::unordered_set<id_set, id_set_hasher> get_face_keys() const {
        std::unordered_set<id_set, id_set_hasher> result;

        std::transform(_data.begin(),
                       _data.end(),
                       std::inserter(result, result.end()),
                       Faces::key_transform);

        return result;
    }

    connectivity_map &get_face_nodes() { return _data; }

    std::unordered_map<id_set, connected_pair, id_set_hasher>
    generate(std::vector<std::vector<id_list>> const &list) {
        std::unordered_set<id_set, id_set_hasher>
          existing_face_keys = get_face_keys();

        std::unordered_map<id_set, connected_pair, id_set_hasher> face_registry;

        std::list<id_list> generated_faces;

        cgsize_t gen_face_idx = 1;
        for (auto const &faces_nodes: list) {
            for (auto const &face_nodes: faces_nodes) {
                id_set face_key = Faces::key_transform(face_nodes);
                if (! existing_face_keys.contains(face_key) &&
                    ! face_registry.contains(face_key)) {
                    // New internal face
                    generated_faces.push_back(face_nodes);
                    face_registry[face_key] = { gen_face_idx++, face_nodes };
                }
            }
        }

        auto boundary_faces_begin = boost::make_zip_iterator(boost::make_tuple(
          boost::make_counting_iterator<cgsize_t>(gen_face_idx),
          _data.begin()));
        auto boundary_faces_end = boost::make_zip_iterator(boost::make_tuple(
          boost::make_counting_iterator<cgsize_t>(generated_faces.size() +
                                                  _data.size() + 1),
          get_face_nodes().end()));

        std::for_each(boundary_faces_begin,
                      boundary_faces_end,
                      [&](auto const &elem) {
                          auto [idx, nodes] = elem;

                          id_set face_key = Faces::key_transform(nodes);

                          face_registry[face_key] = { idx, nodes };
                      });

        connectivity_map new_data(generated_faces.size() + _data.size());

        std::copy(generated_faces.begin(),
                  generated_faces.end(),
                  new_data.begin());
        std::copy(_data.begin(),
                  _data.end(),
                  new_data.begin() + generated_faces.size());

        _data = std::move(new_data);
        _built_face_count += generated_faces.size();

        BOOST_LOG_TRIVIAL(info) << "Built " << _built_face_count << " faces";

        return face_registry;
    }

    static bool
    check_orientation(id_list const &face_nodes,
                      id_list const &reference_nodes) {
        if (face_nodes.size() != reference_nodes.size()) { return false; }

        id_list search = face_nodes;
        search.insert(search.end(), face_nodes.begin(), face_nodes.end());
        size_t window = face_nodes.size();

        for (size_t i = 0; i < face_nodes.size(); ++i) {
            if (std::equal(search.begin() + i,
                           search.begin() + i + window,
                           reference_nodes.begin())) {
                return true;
            }
        }
        return false;
    }

    void reorder_faces(id_list const &map, cgsize_t offset = 0) {
        if (map.size() + offset != _data.size()) {
            BOOST_LOG_TRIVIAL(fatal) << "invalid ordering map size";
            throw std::runtime_error("invalid ordering map size");
        }

        connectivity_map reordered(_data.size());
        std::transform(map.begin(),
                       map.end(),
                       reordered.begin(),
                       [&](size_t old_idx) {
                           auto it = _data.begin();
                           std::advance(it, old_idx + offset);
                           return *it;
                       });

        std::copy(reordered.begin(), reordered.end(), _data.begin() + offset);
    }


public:
    Faces(std::string const &file_path,
          std::string const &base_name,
          std::string const &zone_name)
        : LoaderBase(file_path, base_name, zone_name)
        , Nodes(file_path, base_name, zone_name)
        , Boundaries(file_path, base_name, zone_name)
        , Elements(file_path, base_name, zone_name) {
        load_data();
    }

    connectivity_map const &get_face_nodes() const { return _data; }

    cgsize_t built_face_count() const { return _built_face_count; }
};
}
}
