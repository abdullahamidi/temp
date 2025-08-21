#pragma once

#include "macro.hpp"

/**
 * @brief CGNS API error checking macro
 * 
 * This macro wraps CGNS API function calls and automatically checks for errors.
 * If the function returns a non-zero value (indicating an error), it logs
 * detailed error information including the function name, file location, line number,
 * and the specific CGNS error message, then throws a runtime exception.
 * 
 * @param fn The CGNS API function to call
 * @param args The arguments to pass to the function (in parentheses)
 * @throws std::runtime_error if the CGNS API call returns an error code
 * 
 * Usage example:
 * @code
 * CG_CHECK(cg_open, ("file.cgns", CG_MODE_READ, &file_id));
 * @endcode
 */
#define CG_CHECK(fn, args)                                                     \
    if (fn args) {                                                             \
        BOOST_LOG_TRIVIAL(fatal)                                               \
          << "Error while executing " << QUOTE(fn)                             \
          << ".\n\tfile: " << __FILE__ << ".\n\tline: " << __LINE__            \
          << "\n\terror: " << cg_get_error();                                  \
        throw std::runtime_error("CGNS API call returned error");              \
    }
