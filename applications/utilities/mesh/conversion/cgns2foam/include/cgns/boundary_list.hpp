#pragma once
#include <cgnslib.h>
#include <algorithm>
#include <boost/json/array.hpp>
#include <boost/json/object.hpp>
#include <boost/json/value.hpp>
#include <iterator>
#include <string>
#include <boost/json.hpp>

#include "base.hpp"
#include "boundary.hpp"

namespace CGNS {

class BoundaryList final
    : public virtual impl::LoaderBase
    , public virtual impl::Boundaries {
public:
    BoundaryList(std::string const &file_path,
                 std::string const &base_name = {},
                 std::string const &zone_name = {})
        : impl::LoaderBase(file_path, base_name, zone_name)
        , impl::Boundaries(file_path, base_name, zone_name) {}
};

boost::json::value to_json(BoundaryList const &obj) {
    boost::json::array result = {};

    std::transform(obj.get_boundaries().begin(),
                   obj.get_boundaries().end(),
                   std::inserter(result, result.end()),
                   [](auto &info) {
                       return boost::json::object {
                           { "name",             info.name },
                           { "type", BCTypeName[info.type] },
                       };
                   });

    return result;
}
}
