/*---------------------------------------------------------------------------*\

                        Open Source CFD Group

\*---------------------------------------------------------------------------*/

#ifndef WRITE_CSV_H
#define WRITE_CSV_H

#include "fileName.H"
#include "label.H"
#include "scalarField.H"
#include "vectorField.H"

#ifndef CSV_SEPERATOR
#define CSV_SEPERATOR "; "
#endif

namespace Foam {

template <typename FieldType>
void write_csv(fileName path, label data_count, std::string firstHeader,
               const FieldType &firstColumn,
               const std::vector<std::string> &headers,
               const Field<scalarField> &data) {
  OFstream ofs(path);

  ofs << firstHeader.c_str();
  for (std::size_t i = 0; i < headers.size(); ++i) {
    ofs << CSV_SEPERATOR << headers[i].c_str();
  }
  ofs << nl;

  label data_index = 0;
  for (; data_index < data_count; ++data_index) {
    ofs << firstColumn[data_index];

    for (auto &d : data) {
      ofs << CSV_SEPERATOR << d[data_index];
    }
    ofs << nl;
  }
}

template <typename FieldType>
void write_csv(fileName path, label data_count, std::string firstHeader,
               const FieldType &firstColumn,
               const std::vector<std::string> &headers,
               const Field<vectorField> &data) {
  OFstream ofs(path);

  ofs << firstHeader.c_str();
  for (std::size_t i = 0; i < headers.size(); ++i) {
    ofs << CSV_SEPERATOR << (headers[i] + " - x").c_str();
    ofs << CSV_SEPERATOR << (headers[i] + " - y").c_str();
    ofs << CSV_SEPERATOR << (headers[i] + " - z").c_str();
  }
  ofs << nl;

  label data_index = 0;
  for (; data_index < data_count; ++data_index) {
    ofs << firstColumn[data_index];

    for (auto &d : data) {
      ofs << CSV_SEPERATOR << d[data_index].x();
      ofs << CSV_SEPERATOR << d[data_index].y();
      ofs << CSV_SEPERATOR << d[data_index].z();
    }
    ofs << nl;
  }
}

} // namespace Foam

#ifdef NoRepository
#include "writeCSV.H"
#endif

#endif
