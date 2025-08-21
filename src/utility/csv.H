/*---------------------------------------------------------------------------*\

    V7320 - Özgün Mühendislik Yazılımları

-------------------------------------------------------------------------------
Description

Authors
    Barış Cumhur
        Turkish Aerospace Industries, Turkey
\*---------------------------------------------------------------------------*/

#ifndef CSV_H
#define CSV_H

#include <algorithm>
#include <fstream>
#include <iosfwd>
#include <string>
#include <type_traits>
#include "IOobject.H"
#include "ISstream.H"
#include "IOList.H"
#include "OSspecific.H"
#include "labelFwd.H"
#include "stdFoam.H"
#include "vector.H"
#include "List.H"
#include "OFstream.H"


#ifndef CSV_SEPERATOR
#    define CSV_SEPERATOR ';'
#endif


namespace {
using namespace Foam;

template <typename T, typename std::enable_if<std::is_same<T, label>::value, int>::type = 0>
void _append_header(OFstream &ofs, bool skip, char const *header) {
    if (! skip) { ofs << header << CSV_SEPERATOR; }
}

template <typename T, typename std::enable_if<std::is_same<T, scalar>::value, int>::type = 0>
void _append_header(OFstream &ofs, bool skip, char const *header) {
    if (! skip) { ofs << header << CSV_SEPERATOR; }
}

template <typename T, typename std::enable_if<std::is_same<T, vector>::value, int>::type = 0>
void _append_header(OFstream &ofs, bool skip, char const *header) {
    if (! skip) {
        ofs << (std::string(header) + " - x").c_str() << CSV_SEPERATOR;
        ofs << (std::string(header) + " - y").c_str() << CSV_SEPERATOR;
        ofs << (std::string(header) + " - z").c_str() << CSV_SEPERATOR;
    }
}

template <size_t J, size_t I = 0>
struct _append_headers {
    template <typename First, typename... Rest, typename... S>
    static void
    perform(OFstream &ofs, std::array<bool, J> const &skips, char const *first, S... rest) {
        _append_header<First>(ofs, skips[I], first);
        _append_headers<J, I + 1>::template perform<Rest...>(ofs, skips, rest...);
    }
};

template <size_t J>
struct _append_headers<J, J - 1> {
    template <typename First>
    static void perform(OFstream &ofs, std::array<bool, J> const &skips, char const *first) {
        _append_header<First>(ofs, skips[J - 1], first);
    }
};

template <typename... FieldType, typename... S>
void append_headers(OFstream &ofs, std::array<bool, sizeof...(S)> const &skips, S... headers) {
    _append_headers<sizeof...(FieldType)>::template perform<typename std::decay<
      decltype(std::declval<FieldType>()[std::declval<int>()])>::type...>(ofs, skips, headers...);
}

template <typename T, typename std::enable_if<std::is_same<T, label>::value, int>::type = 0>
void __append_data(OFstream &ofs, bool skip, T const &data) {
    if (! skip) { ofs << data << CSV_SEPERATOR; }
}

template <typename T, typename std::enable_if<std::is_same<T, scalar>::value, int>::type = 0>
void __append_data(OFstream &ofs, bool skip, T const &data) {
    if (! skip) { ofs << data << CSV_SEPERATOR; }
}

template <typename T, typename std::enable_if<std::is_same<T, vector>::value, int>::type = 0>
void __append_data(OFstream &ofs, bool skip, T const &data) {
    if (! skip) {
        ofs << data.x() << CSV_SEPERATOR;
        ofs << data.y() << CSV_SEPERATOR;
        ofs << data.z() << CSV_SEPERATOR;
    }
}

template <size_t J, size_t I = 0>
struct _append_data {
    template <typename First, typename... Rest>
    static void
    perform(OFstream &ofs,
            std::array<bool, J> const &skips,
            First const &first,
            Rest const &...rest) {
        __append_data<First>(ofs, skips[I], first);
        _append_data<J, I + 1>::template perform<Rest...>(ofs, skips, rest...);
    }
};

template <size_t J>
struct _append_data<J, J - 1> {
    template <typename T>
    static void perform(OFstream &ofs, std::array<bool, J> const &skips, T const &first) {
        __append_data<T>(ofs, skips[J - 1], first);
    }
};

template <typename... T>
void append_data(OFstream &ofs, std::array<bool, sizeof...(T)> const &skips, T const &...values) {
    _append_data<sizeof...(T)>::perform(ofs, skips, values...);
}
}

namespace Foam {

template <typename First, typename... FieldType>
void writeCsv(IOobject obj,
              labelList const &map,
              bool map_first,
              std::pair<char const *, First *> first,
              std::pair<char const *, FieldType *>... rest) {

    if (! isDir(obj.path())) { mkDir(obj.path()); }

    OFstream ofs(obj.objectPath());

    if (! ofs.opened()) {
        FatalErrorInFunction << "Cannot open file: " << obj.objectPath() << abort(FatalError);
    }

    append_headers<First, FieldType...>(
      ofs,
      { first.second == nullptr, rest.second == nullptr... },
      first.first,
      rest.first...);
    ofs << endl;

    forAll(map, i) {
        label data_index = map[i];
        if (first.first != nullptr) {
            if (map_first) {
                ofs << (*first.second)[data_index] << CSV_SEPERATOR;
            } else {
                ofs << (*first.second)[i] << CSV_SEPERATOR;
            }
        }
        append_data(ofs,
                    { rest.second == nullptr... },
                    (rest.second != nullptr
                       ? (*rest.second)[data_index]
                       : typename std::decay<decltype((*rest.second)[data_index])>::type {})...);
        ofs << endl;
    }

    ofs.flush();
}

template <typename... FieldType>
void writeCsv(IOobject obj, label data_count, std::pair<char const *, FieldType *>... fields) {
    labelList map(data_count);
    forAll(map, i) { map[i] = i; }
    writeCsv(obj, map, true, fields...);
}
}

namespace {
using namespace Foam;

template <typename T, typename std::enable_if<std::is_same<T, label>::value, int>::type = 0>
void _check_header(std::istringstream &iss, bool skip, char const *header) {
    if (skip) { return; }
    std::string read_header;
    std::getline(iss, read_header, CSV_SEPERATOR);
    if (read_header != header) {
        FatalErrorInFunction << "Header mismatch: expected \"" << header << "\", but found \""
                             << read_header << "\"" << exit(FatalError);
    }
}

template <typename T, typename std::enable_if<std::is_same<T, scalar>::value, int>::type = 0>
void _check_header(std::istringstream &iss, bool skip, char const *header) {
    if (skip) { return; }

    std::string read_header;
    std::getline(iss, read_header, CSV_SEPERATOR);
    if (read_header != header) {
        FatalErrorInFunction << "Header mismatch: expected \"" << header << "\", but found \""
                             << read_header << "\"" << exit(FatalError);
    }
}

template <typename T, typename std::enable_if<std::is_same<T, vector>::value, int>::type = 0>
void _check_header(std::istringstream &iss, bool skip, char const *header) {
    if (skip) { return; }

    std::string read_header_x;
    std::string read_header_y;
    std::string read_header_z;
    std::getline(iss, read_header_x, CSV_SEPERATOR);
    std::getline(iss, read_header_y, CSV_SEPERATOR);
    std::getline(iss, read_header_z, CSV_SEPERATOR);
    if (read_header_x != std ::string(header) + " - x" ||
        read_header_y != std ::string(header) + " - y" ||
        read_header_z != std ::string(header) + " - z") {
        FatalErrorInFunction
          << "Header mismatch: expected \""
          << "(" + std ::string(header) + " - x" + ", " + std ::string(header) + " - y" + ", " +
               std ::string(header) + " - z" + ")"
          << "\", but found \""
          << "(" + read_header_x + ", " + read_header_y + ", " + read_header_z + ")" << "\""
          << exit(FatalError);
    }
}

template <size_t J, size_t I = 0>
struct _check_headers {
    template <typename First, typename... Rest, typename... S>
    static void
    perform(std::istringstream &ofs, std::array<bool, J> skips, char const *first, S... rest) {
        _check_header<First>(ofs, skips[I], first);
        _check_headers<J, I + 1>::template perform<Rest...>(ofs, skips, rest...);
    }
};

template <size_t J>
struct _check_headers<J, J - 1> {
    template <typename First>
    static void
    perform(std::istringstream &ofs, std::array<bool, J> const &skips, char const *first) {
        _check_header<First>(ofs, skips[J - 1], first);
    }
};

template <typename... FieldType, typename... S>
void check_headers(std::istringstream &ofs,
                   std::array<bool, sizeof...(FieldType)> skips,
                   S... headers) {
    _check_headers<sizeof...(FieldType)>::template perform<typename std::decay<
      decltype(std::declval<FieldType>()[std::declval<int>()])>::type...>(ofs, skips, headers...);
}

template <typename T, typename std::enable_if<std::is_same<T, label>::value, int>::type = 0>
void __read_data(std::istringstream &iss, bool skip, T *data) {
    if (skip) { return; }
    std::string line;
    std::getline(iss, line, CSV_SEPERATOR);
    line.erase(std::remove_if(line.begin(), line.end(), [](char c) { return std::isspace(c); }),
               line.end());
    *data = static_cast<label>(std::stoll(line.c_str()));
}

template <typename T, typename std::enable_if<std::is_same<T, scalar>::value, int>::type = 0>
void __read_data(std::istringstream &iss, bool skip, T *data) {
    if (skip) { return; }
    std::string line;
    std::getline(iss, line, CSV_SEPERATOR);
    line.erase(std::remove_if(line.begin(), line.end(), [](char c) { return std::isspace(c); }),
               line.end());
    *data = static_cast<scalar>(std::stod(line.c_str()));
}

template <typename T, typename std::enable_if<std::is_same<T, vector>::value, int>::type = 0>
void __read_data(std::istringstream &iss, bool skip, T *data) {
    if (skip) { return; }
    std::string line_x, line_y, line_z;
    std::getline(iss, line_x, CSV_SEPERATOR);
    std::getline(iss, line_y, CSV_SEPERATOR);
    std::getline(iss, line_z, CSV_SEPERATOR);
    line_x
      .erase(std::remove_if(line_x.begin(), line_x.end(), [](char c) { return std::isspace(c); }),
             line_x.end());
    line_y
      .erase(std::remove_if(line_y.begin(), line_y.end(), [](char c) { return std::isspace(c); }),
             line_y.end());
    line_z
      .erase(std::remove_if(line_z.begin(), line_z.end(), [](char c) { return std::isspace(c); }),
             line_z.end());
    data->x() = static_cast<scalar>(std::stod(line_x.c_str()));
    data->y() = static_cast<scalar>(std::stod(line_y.c_str()));
    data->z() = static_cast<scalar>(std::stod(line_z.c_str()));
}

template <size_t J, size_t I = 0>
struct _read_util {
    template <typename First, typename... Rest>
    static void
    perform(std::istringstream &iss, std::array<bool, J> const &skips, First *first, Rest *...rest) {
        __read_data<First>(iss, skips[I], first);
        _read_util<J, I + 1>::template perform<Rest...>(iss, skips, rest...);
    }

    template <typename First, typename... Rest>
    static void
    resize(size_t count, std::array<bool, J> const &skips, First *first, Rest *...rest) {
        if (! skips[I]) { first->resize(count); }
        _read_util<J, I + 1>::template resize<Rest...>(count, skips, rest...);
    }
};

template <size_t J>
struct _read_util<J, J - 1> {
    template <typename T>
    static void perform(std::istringstream &iss, std::array<bool, J> const &skips, T *first) {
        __read_data<T>(iss, skips[J - 1], first);
    }

    template <typename First>
    static void resize(size_t count, std::array<bool, J> const &skips, First *first) {
        if (! skips[J - 1]) { first->resize(count); }
    }
};

template <typename... T>
void read_data(std::istringstream &iss, std::array<bool, sizeof...(T)> const &skips, T *...values) {
    _read_util<sizeof...(T)>::perform(iss, skips, values...);
}

template <typename... T>
void resize_fields(size_t count, std::array<bool, sizeof...(T)> const &skips, T *...fields) {
    _read_util<sizeof...(T)>::resize(count, skips, fields...);
}
}

namespace Foam {
template <typename... FieldType>
void readCsv(IOobject obj, std::pair<char const *, FieldType *>... fields) {

    if (! isDir(obj.path())) { mkDir(obj.path()); }

    std::fstream ifs(obj.objectPath().c_str());

    if (! ifs.is_open()) {
        FatalErrorInFunction << "Cannot open file: " << obj.objectPath() << abort(FatalError);
    }

    std::array<bool, sizeof...(FieldType)> const skips = { fields.second == nullptr... };

    string line;
    std::getline(ifs, line);
    std::istringstream iss(line);
    check_headers<FieldType...>(iss, skips, fields.first...);

    auto state = ifs.rdstate();
    std::streampos pos = ifs.tellg();
    size_t count = std::count(std::istreambuf_iterator<char>(ifs),
                              std::istreambuf_iterator<char>(),
                              '\n');
    ifs.seekg(pos);
    ifs.setstate(state);

    resize_fields(count, skips, fields.second...);


    for (size_t i = 0; i < count; ++i) {
        std::getline(ifs, line);
        iss = std::istringstream(line);
        read_data(iss, skips, fields.second != nullptr ? &((*fields.second)[i]) : nullptr...);
    }
}
}

#endif
