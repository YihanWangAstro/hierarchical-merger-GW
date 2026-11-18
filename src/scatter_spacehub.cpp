// BH/star binary scattering -- version 3: SpaceHub.
//
// Uses SpaceHub's default method: algorithmic-regularised chain with
// Bulirsch-Stoer extrapolation and compensated summation (AR_Chain_Plus).
// SpaceHub works in units G = 1, Msun, AU, so one year is 2 pi time units.

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <limits>
#include <vector>

// SpaceHub gives up on a step it cannot converge by printing a message and
// calling exit(), which would end a whole Monte Carlo run because of a single
// hard encounter.  While its headers are being read, exit() is therefore
// turned into an exception that the wrapper below can catch.
namespace {
struct SpaceHubAbort {};
[[noreturn]] inline void throwSpaceHubAbort(int) { throw SpaceHubAbort{}; }
}  // namespace
#define exit(code) throwSpaceHubAbort(code)
#include "spaceHub.hpp"   // from the SpaceHub repository (see the Makefile)
#undef exit

#include "bhs/driver.hpp"

namespace {

using namespace bhs;

class SpaceHubIntegrator {
   public:
    static const char* name() { return "spacehub"; }
    static const char* description() { return "SpaceHub AR-chain"; }
    static double defaultAccuracy() { return 1e-12; }
    static const char* accuracyHelp() { return "relative error tolerance per step                  "; }

    explicit SpaceHubIntegrator(double tolerance) : tolerance_(tolerance) {}

    IntegrationStats integrate(const Masses& m, State& y, double tEnd, const BinaryScales&, Sampler* sampler) const {
        IntegrationStats stats;
        if (!sampler) {
            stats.ok = advance(m, y, tEnd, stats);
            return stats;
        }
        double t = 0.0;
        while (stats.ok && !sampler->done()) {
            const double next = sampler->nextTime();
            stats.ok = advance(m, y, next - t, stats);
            if (stats.ok) sampler->record(y);
            t = next;
        }
        return stats;
    }

   private:
    using Solver = hub::methods::DefaultMethod<>;
    using Particle = Solver::Particle;

    static constexpr double kTimeUnitsPerYear = kTwoPi;
    static constexpr double kEndTimeTolerance = 1e-12;   // relative accuracy of the stop time
    static constexpr int kMaxRetries = 3;                // each retry loosens the tolerance tenfold
    static constexpr double kFirstStepShare = 0.1;       // first step, as a share of the shortest time-scale
    static constexpr int kMaxShorterFirstSteps = 6;      // each one a tenth of the one before

    /// Advance y by `duration`.  If SpaceHub cannot converge at the requested
    /// tolerance, start again from the same state with a looser one.
    ///
    /// The first step is chosen here.  SpaceHub would take a tenth of the
    /// shortest free-fall time of any pair, and its Bulirsch-Stoer iterator
    /// does not shorten a first step that fails to converge: it gives up.  A
    /// body that crosses the system in less than a free-fall time (a fast
    /// star passing a wide binary) therefore needs a shorter first step, and
    /// one that still fails is shortened, at the same tolerance.
    bool advance(const Masses& m, State& y, double duration, IntegrationStats& stats) const {
        double tolerance = tolerance_;
        double firstStep = kFirstStepShare * shortestTimescale(m, y);
        for (int loosened = 0, shortened = 0;;) {
            const long stepsBefore = stats.steps;
            try {
                Solver solver{0.0, particlesOf(m, y)};
                Solver::RunArgs args;
                args.rtol = tolerance;
                args.time_rtol = kEndTimeTolerance;
                args.step_size = firstStep * kTimeUnitsPerYear * solver.particles().step_scale();
                args.add_operation([&stats](auto&, auto) { ++stats.steps; });
                args.add_stop_condition(duration * kTimeUnitsPerYear);
                solver.run(args);
                y = stateOf(solver);
                return true;
            } catch (const SpaceHubAbort&) {
                const bool atFirstStep = stats.steps - stepsBefore <= 1;
                if (atFirstStep && shortened < kMaxShorterFirstSteps) {
                    ++shortened;
                    firstStep *= 0.1;
                } else if (loosened < kMaxRetries) {
                    ++loosened;
                    ++stats.retries;
                    tolerance *= 10.0;
                } else {
                    return false;
                }
            }
        }
    }

    /// The shortest time [yr] on which a pair of bodies changes its
    /// separation: free fall from rest, or crossing at the present speed.
    static double shortestTimescale(const Masses& m, const State& y) {
        double shortest = std::numeric_limits<double>::infinity();
        for (int i = 0; i < 3; ++i) {
            for (int j = i + 1; j < 3; ++j) {
                const double r = norm(position(y, j) - position(y, i));
                const double speed = norm(velocity(y, j) - velocity(y, i));
                const double freeFall = 0.5 * kPi * std::sqrt(r * r * r / (2.0 * kG * (m[i] + m[j])));
                shortest = std::min({shortest, freeFall, r / speed});
            }
        }
        return shortest;
    }

    static std::vector<Particle> particlesOf(const Masses& m, const State& y) {
        std::vector<Particle> particles;
        for (int i = 0; i < 3; ++i) {
            const Vec3 r = position(y, i), v = velocity(y, i) / kTimeUnitsPerYear;
            particles.emplace_back(m[i], r.x, r.y, r.z, v.x, v.y, v.z);
        }
        return particles;
    }

    static State stateOf(const Solver& solver) {
        const auto& particles = solver.particles();
        State y;
        for (int i = 0; i < 3; ++i) {
            const auto& r = particles.pos()[i];
            const auto& v = particles.vel()[i];
            setBody(y, i, {static_cast<double>(r.x), static_cast<double>(r.y), static_cast<double>(r.z)},
                    kTimeUnitsPerYear * Vec3{static_cast<double>(v.x), static_cast<double>(v.y), static_cast<double>(v.z)});
        }
        return y;
    }

    double tolerance_;
};

}  // namespace

int main(int argc, char** argv) { return bhs::runMain<SpaceHubIntegrator>(argc, argv); }
