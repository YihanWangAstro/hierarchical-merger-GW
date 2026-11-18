// Nuclear star cluster model and the timescales that drive the Monte Carlo.
#pragma once

#include <cmath>
#include <vector>

#include "constants.hpp"

namespace bhs {

/// The cluster as a binary sees it at one distance from the central black
/// hole: a power-law cusp of stars and black holes.
struct LocalCluster {
    double radius = 0.0;           // [AU]
    double starDensity = 0.0;      // [AU^-3]
    double bhDensity = 0.0;        // [AU^-3]
    double starDispersion = 0.0;   // [AU/yr]
    double bhDispersion = 0.0;     // [AU/yr]

    double radiusPc() const { return radius / kParsec; }

    /// Escape speed from the central black hole alone.
    double escapeSpeed() const { return std::sqrt(2.0 * kG * kSmbhMass / radius); }
};

inline LocalCluster clusterAt(double radius) {
    const double norm = kDensityNorm / (kStarMass * kParsec * kParsec * kParsec);
    auto density = [&](double slope) {
        return norm * std::pow(1.0 / kDensityNormRadius, -slope) * std::pow(radius / kParsec, -slope);
    };
    auto dispersion = [&](double slope) { return std::sqrt(kG * kSmbhMass / (radius * (1.0 + slope))); };
    return {radius, density(kCuspSlope), kBhNumberFraction * density(kBhCuspSlope), dispersion(kCuspSlope),
            dispersion(kBhCuspSlope)};
}

/// The log-spaced radii [AU] tabulated by the `profile` command.
inline std::vector<double> profileRadii() {
    const double logInner = std::log(kInnerRadius * kParsec), logOuter = std::log(kOuterRadius * kParsec);
    std::vector<double> radii(kNumRadii);
    for (int i = 0; i < kNumRadii; ++i) radii[i] = std::exp(logInner + (logOuter - logInner) * i / (kNumRadii - 1));
    return radii;
}

/// Mean mass of the field (1G) black holes, dN/dm ~ m^-2.
inline double meanFieldBhMass() {
    return std::log(kMaxFieldBhMass / kMinFieldBhMass) / (1.0 / kMinFieldBhMass - 1.0 / kMaxFieldBhMass);
}

/// Gravitational-wave merger time of the binary m1, m2 (Peters 1964,
/// high-eccentricity form).
inline double mergerTime(double m1, double m2, double a, double e) {
    const double total = m1 + m2;
    const double reduced = m1 * m2 / total;
    return 5.0 * std::pow(kSpeedOfLight, 5) * std::pow(a, 4) / (256.0 * kG * kG * kG * total * total * reduced) *
           std::pow(1.0 - e * e, 3.5);
}

/// Mean time between encounters of a binary (semi-major axis a) with members
/// of a population of given number density and velocity dispersion, for
/// pericentres within 2a (gravitational-focusing limit).  totalMass is the
/// mass of the binary plus the intruder.
inline double encounterTime(double totalMass, double density, double dispersion, double a) {
    return 1.0 / (kTwoPi * kG * totalMass * density / dispersion * 2.0 * a);
}

/// Distance from a black hole (mass bhMass) inside which its tide pulls a
/// star apart: r_t = R_star (M_bh / m_star)^(1/3).
inline double tidalRadius(double bhMass) {
    return kStarRadius * std::cbrt(bhMass / kStarMass);
}

/// True if a star on the orbit (a, e) about a black hole comes within
/// `factor` tidal radii of it.
inline bool tidallyDisrupted(double bhMass, double a, double e, double factor) {
    return a * (1.0 - e) <= factor * tidalRadius(bhMass);
}

}  // namespace bhs
