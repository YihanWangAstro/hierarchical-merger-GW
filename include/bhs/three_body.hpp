// The Newtonian three-body problem: state layout, equations of motion, and the
// small interface every integrator in this project implements.
#pragma once

#include <array>
#include <cmath>
#include <vector>

#include "constants.hpp"
#include "vec3.hpp"

namespace bhs {

using Masses = std::array<double, 3>;

/// Phase-space point: (r1, r2, r3, v1, v2, v3), three components each.
using State = std::array<double, 18>;

inline Vec3 position(const State& y, int body) { return {y[3 * body], y[3 * body + 1], y[3 * body + 2]}; }
inline Vec3 velocity(const State& y, int body) { return {y[9 + 3 * body], y[10 + 3 * body], y[11 + 3 * body]}; }

inline void setBody(State& y, int body, const Vec3& r, const Vec3& v) {
    y[3 * body] = r.x, y[3 * body + 1] = r.y, y[3 * body + 2] = r.z;
    y[9 + 3 * body] = v.x, y[10 + 3 * body] = v.y, y[11 + 3 * body] = v.z;
}

inline bool allFinite(const State& y) {
    for (double v : y)
        if (!std::isfinite(v)) return false;
    return true;
}

/// Time derivative of the state under mutual Newtonian gravity.
inline void derivative(const Masses& m, const State& y, State& dydt) {
    for (int i = 0; i < 9; ++i) dydt[i] = y[9 + i];
    for (int i = 9; i < 18; ++i) dydt[i] = 0.0;
    for (int i = 0; i < 3; ++i) {
        for (int j = i + 1; j < 3; ++j) {
            const double dx = y[3 * j] - y[3 * i];
            const double dy = y[3 * j + 1] - y[3 * i + 1];
            const double dz = y[3 * j + 2] - y[3 * i + 2];
            const double r2 = dx * dx + dy * dy + dz * dz;
            const double gOverR3 = kG / (r2 * std::sqrt(r2));
            dydt[9 + 3 * i] += m[j] * gOverR3 * dx;
            dydt[10 + 3 * i] += m[j] * gOverR3 * dy;
            dydt[11 + 3 * i] += m[j] * gOverR3 * dz;
            dydt[9 + 3 * j] -= m[i] * gOverR3 * dx;
            dydt[10 + 3 * j] -= m[i] * gOverR3 * dy;
            dydt[11 + 3 * j] -= m[i] * gOverR3 * dz;
        }
    }
}

/// Total energy (kinetic + potential).
inline double totalEnergy(const Masses& m, const State& y) {
    double energy = 0.0;
    for (int i = 0; i < 3; ++i) energy += 0.5 * m[i] * dot(velocity(y, i), velocity(y, i));
    for (int i = 0; i < 3; ++i)
        for (int j = i + 1; j < 3; ++j) energy -= kG * m[i] * m[j] / norm(position(y, j) - position(y, i));
    return energy;
}

// ----------------------------------------------------------------------------
// Integrator interface
// ----------------------------------------------------------------------------
//
// An integrator is any class providing
//
//     static const char* name();            short id, used in file names ("bs")
//     static const char* description();     one line for --help
//     static double defaultAccuracy();      default of the --acc option
//     static const char* accuracyHelp();    what --acc means for this method
//     explicit Integrator(double accuracy);
//     IntegrationStats integrate(const Masses&, State& y, double tEnd,
//                                const BinaryScales&, Sampler* sampler) const;
//
// integrate() advances y from t = 0 to exactly t = tEnd.  If `sampler` is not
// null it must also hand over the state at every sampler time (see Sampler).

/// Characteristic scales of the tight binary, used to choose steps and tolerances.
struct BinaryScales {
    double semiMajorAxis = 1.0;
    double eccentricity = 0.0;
    double period = 1.0;

    double pericentre() const { return semiMajorAxis * (1.0 - eccentricity); }
    double orbitalSpeed() const { return kTwoPi * semiMajorAxis / period; }
};

struct IntegrationStats {
    bool ok = true;    // false if the integrator gave up
    long steps = 0;
    int retries = 0;   // times the integrator had to start again with a looser tolerance
};

struct TrajectoryPoint {
    double t;
    State y;
};

/// Collects the state at the equally spaced times t_k = k * tEnd / n, k = 0..n.
/// Whoever owns the integration records the points in order.
class Sampler {
   public:
    Sampler(double tEnd, int intervals) : tEnd_(tEnd), intervals_(intervals) {}

    bool done() const { return static_cast<int>(points_.size()) > intervals_; }

    /// Time of the last point.
    double endTime() const { return tEnd_; }

    /// Time of the next point that has not been recorded yet.
    double nextTime() const { return tEnd_ * static_cast<double>(points_.size()) / intervals_; }

    /// Record the state at nextTime().
    void record(const State& y) { points_.push_back({nextTime(), y}); }

    const std::vector<TrajectoryPoint>& points() const { return points_; }

   private:
    double tEnd_;
    int intervals_;
    std::vector<TrajectoryPoint> points_;
};

}  // namespace bhs
