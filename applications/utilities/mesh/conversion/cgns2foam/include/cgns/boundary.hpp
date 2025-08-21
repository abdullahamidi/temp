
#pragma once

#include <cgns_io.h>
#include <cgnslib.h>
#include <cgnstypes.h>
#include <functional>
#include <numeric>
#include <vector>
#include <string>
#include <algorithm>
#include <stdexcept>
#include "base.hpp"

namespace CGNS {
namespace impl {
class Boundaries : public virtual LoaderBase {
    struct boundary_info {
        std::string name;

        BCType_t type;

        PointSetType_t ps_type;

        std::vector<cgsize_t> points;

        std::pair<cgsize_t, cgsize_t> remap;

        GridLocation_t location;

        int ndataset;
    };

    std::vector<boundary_info> _data;

    cgsize_t _boundary_count;

    void read_boundary_data() {
        int nbocos;

        cg_nbocos(get_file(), get_base_index(), get_zone_index(), &nbocos);

        _data.resize(nbocos);

        for (int iboco = 1; iboco <= nbocos; ++iboco) {
            cgsize_t npnts;
            boundary_info boco = {};
            int inormal;
            cgsize_t nnormal;
            DataType_t normal_type;
            char boco_name[CGIO_MAX_NAME_LENGTH + 1];

            CG_CHECK(
              cg_boco_info,
              (get_file(),
               get_base_index(),
               get_zone_index(),
               iboco,
               boco_name,
               &boco.type,
               &boco.ps_type,
               &npnts,
               &inormal,
               &nnormal,
               &normal_type,
               &boco.ndataset));
            boco.name = boco_name;

            CG_CHECK(cg_boco_gridlocation_read,
                     (get_file(),
                      get_base_index(),
                      get_zone_index(),
                      iboco,
                      &boco.location));

            if (boco.location != GridLocation_t::FaceCenter) {
                BOOST_LOG_TRIVIAL(fatal)
                  << "Boundaries required to be defined on face centers";
                throw std::runtime_error("Boundaries required to be defined on "
                                         "face centers");
            }

            boco.points.resize(npnts);

            CG_CHECK(cg_boco_read,
                     (get_file(),
                      get_base_index(),
                      get_zone_index(),
                      iboco,
                      boco.points.data(),
                      NULL));

            _data[iboco - 1] = std::move(boco);
        }

        std::sort(_data.begin(),
                  _data.end(),
                  [](boundary_info const &a, boundary_info const &b) {
                      return a.points[0] < b.points[0];
                  });

        cgsize_t first = std::transform_reduce(
          _data.begin(),
          _data.end(),
          std::numeric_limits<cgsize_t>::max(),
          [&](cgsize_t result, cgsize_t current) {
              return std::min(result, current);
          },
          [](auto const &bd) {
              return std::reduce(
                bd.points.begin(),
                bd.points.end(),
                std::numeric_limits<cgsize_t>::max(),
                [](cgsize_t a, cgsize_t b) { return std::min(a, b); });
          });

        std::for_each(_data.begin(), _data.end(), [&](auto &bd) {
            std::transform(bd.points.begin(),
                           bd.points.end(),
                           bd.points.begin(),
                           [&](cgsize_t p) { return p - first; });
        });


        BOOST_LOG_TRIVIAL(info) << "Read " << nbocos << " boundaries";
    }

protected:
    std::vector<boundary_info> &get_boundaries() { return _data; }

    cgsize_t get_defined_boundary_face_count() const {
        return std::transform_reduce(
          _data.begin(),
          _data.end(),
          cgsize_t {},
          std::plus {},
          [](auto const &info) {
              if (info.ps_type == PointSetType_t::ElementList ||
                  info.ps_type == PointSetType_t::PointList) {
                  return (cgsize_t) info.points.size();
              }
              if (info.ps_type == PointSetType_t::ElementRange ||
                  info.ps_type == PointSetType_t::PointRange) {
                  return info.points[1] - info.points[0] + 1;
              }
              return cgsize_t {};
          });
    }

public:
    Boundaries(std::string const &file_path,
               std::string const &base_name,
               std::string const &zone_name)
        : impl::LoaderBase(file_path, base_name, zone_name) {
        read_boundary_data();
    }

    std::vector<boundary_info> const &get_boundaries() const { return _data; }

    cgsize_t get_defined_boundary_count() const { return _boundary_count; }
};

}

}
