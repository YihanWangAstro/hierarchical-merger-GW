// Writing results as Mathematica input, so that the notebooks' analysis cells
// can take them as they are: Get["code0_spacehub.m"].
#pragma once

#include <cmath>
#include <string>

#include "io.hpp"

namespace bhs::mma {

/// A machine real the way Mathematica writes one: 0.31735951552887776` or
/// 6.937361306167924`*^8.  (Without the backtick a 17-digit literal would be
/// read as an arbitrary-precision number.)
inline std::string real(double x) {
    if (std::isnan(x)) return "Indeterminate";
    if (std::isinf(x)) return x > 0 ? "Infinity" : "-Infinity";
    if (x == 0.0) return "0.";

    std::string scientific;   // shortest d.ddd...e+XX that reads back exactly
    for (int decimals = 14; decimals <= 16; ++decimals) {
        scientific = format("%.*e", decimals, x);
        if (std::strtod(scientific.c_str(), nullptr) == x) break;
    }
    const std::size_t e = scientific.find('e');
    const int power = std::stoi(scientific.substr(e + 1));
    std::string mantissa = scientific.substr(0, e);
    mantissa.erase(mantissa.find_last_not_of('0') + 1);   // trailing zeros; at least "d." is left

    if (power < -5 || power >= 6) return mantissa + "`*^" + std::to_string(power);
    const int digits = static_cast<int>(mantissa.size()) - (x < 0 ? 2 : 1);   // significant digits
    const int decimals = digits - 1 - power;
    return decimals > 0 ? format("%.*f`", decimals, x) : format("%.0f.`", x);
}

inline std::string integer(long n) { return std::to_string(n); }
inline std::string text(const std::string& s) { return '"' + s + '"'; }

/// "{f(x1),f(x2),...}" for the items of a range.
template <class Range, class F>
std::string list(const Range& items, F f) {
    std::string out = "{";
    for (const auto& item : items) {
        if (out.size() > 1) out += ',';
        out += f(item);
    }
    return out + "}";
}

/// "{f(0),f(1),...,f(count-1)}".
template <class F>
std::string table(int count, F f) {
    std::string out = "{";
    for (int i = 0; i < count; ++i) {
        if (i > 0) out += ',';
        out += f(i);
    }
    return out + "}";
}

/// "name=value;" on a line of its own.
inline std::string assign(const std::string& name, const std::string& value) { return name + "=" + value + ";\n"; }

}  // namespace bhs::mma
