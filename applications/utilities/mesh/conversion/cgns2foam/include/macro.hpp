#pragma once

/**
 * @brief Stringification macro - converts argument to string literal
 * 
 * @param x The token to be converted to a string
 */
#define Q(x) #x

/**
 * @brief Double-level stringification macro for macro expansion
 * 
 * This macro first expands the argument, then converts it to a string.
 * Useful for stringifying macro values.
 * 
 * @param x The macro or token to be expanded and converted to string
 */
#define QUOTE(x) Q(x)

/**
 * @brief Dynamic assertion macro with logging support
 * 
 * Evaluates the given expression and throws a runtime error with logging
 * if the assertion fails. Uses BOOST_LOG_TRIVIAL for fatal error logging.
 * 
 * @param expr The expression to evaluate for assertion
 * @throws std::runtime_error if the assertion fails
 */
#define DYN_ASSERT(expr)                                                       \
    if (! (expr)) {                                                            \
        BOOST_LOG_TRIVIAL(fatal) << "Error while assertion " << QUOTE(expr);   \
        throw std::runtime_error("Unable to validate assertion");              \
    }
