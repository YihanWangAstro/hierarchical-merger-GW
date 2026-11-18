// BH/star binary scattering -- version 2: time-transformed leapfrog.
//
// A plain leapfrog with a fixed time step cannot follow these binaries
// (eccentricities up to 0.999), so the step is taken in a fictitious time s
// with ds = Omega dt, Omega = sum over pairs of 1 / r_ij (Mikkola & Aarseth
// 2002).  The physical step shrinks automatically at close approaches while
// the scheme stays a second-order, time-reversible drift-kick-drift leapfrog;
// an isolated binary is followed along its exact ellipse with only a phase
// error.

#include <algorithm>
#include <cmath>

#include "bhs/driver.hpp"

namespace {

using namespace bhs;

class TimeTransformedLeapfrog {
   public:
    static const char* name() { return "leapfrog"; }
    static const char* description() { return "time-transformed leapfrog (2nd order)"; }
    static double defaultAccuracy() { return 1000.0; }
    static const char* accuracyHelp() { return "steps per orbit of the initial binary              "; }

    explicit TimeTransformedLeapfrog(double stepsPerOrbit) : stepsPerOrbit_(stepsPerOrbit) {}

    IntegrationStats integrate(const Masses& m, State& y, double tEnd, const BinaryScales& binary,
                               Sampler* sampler) const {
        // Over one period of the binary the fictitious time advances by
        // integral(dt / r) = period / a.
        const double h = binary.period / binary.semiMajorAxis / stepsPerOrbit_;

        IntegrationStats stats;
        Phase now{y, 0.0, omega(y)};
        while (true) {
            Phase next = now;
            step(m, next, h);
            if (!(next.w > 0.0) || ++stats.steps > kMaxSteps) {
                stats.ok = false;
                break;
            }
            // Output points inside this step are side steps from `now`, so
            // that recording a trajectory does not alter the integration.
            while (sampler && !sampler->done() && sampler->nextTime() <= std::min(next.t, tEnd))
                sampler->record(advanceTo(m, now, sampler->nextTime()).y);
            if (next.t >= tEnd) {
                now = advanceTo(m, now, tEnd);
                break;
            }
            now = next;
        }
        y = now.y;
        return stats;
    }

   private:
    static constexpr long kMaxSteps = 2000000000;

    /// State extended by the physical time and by W, the integrated twin of Omega.
    struct Phase {
        State y;
        double t;
        double w;
    };

    static double omega(const State& y) {
        double sum = 0.0;
        for (int i = 0; i < 3; ++i)
            for (int j = i + 1; j < 3; ++j) sum += 1.0 / norm(position(y, j) - position(y, i));
        return sum;
    }

    static void drift(Phase& p, double h) {
        const double dt = h / p.w;
        for (int i = 0; i < 9; ++i) p.y[i] += dt * p.y[9 + i];
        p.t += dt;
    }

    static void kick(const Masses& m, Phase& p, double h) {
        double acceleration[9] = {}, gradOmega[9] = {};
        double omegaNow = 0.0;
        for (int i = 0; i < 3; ++i) {
            for (int j = i + 1; j < 3; ++j) {
                double d[3], r2 = 0.0;
                for (int c = 0; c < 3; ++c) {
                    d[c] = p.y[3 * j + c] - p.y[3 * i + c];
                    r2 += d[c] * d[c];
                }
                const double invR = 1.0 / std::sqrt(r2);
                const double invR3 = invR * invR * invR;
                omegaNow += invR;
                for (int c = 0; c < 3; ++c) {
                    acceleration[3 * i + c] += kG * m[j] * invR3 * d[c];
                    acceleration[3 * j + c] -= kG * m[i] * invR3 * d[c];
                    gradOmega[3 * i + c] += invR3 * d[c];
                    gradOmega[3 * j + c] -= invR3 * d[c];
                }
            }
        }
        // dv/ds = a / Omega,  dW/ds = (v . grad Omega) / Omega with v averaged over the kick.
        const double dt = h / omegaNow;
        double vDotGrad = 0.0;
        for (int i = 0; i < 9; ++i) {
            const double vOld = p.y[9 + i];
            const double vNew = vOld + dt * acceleration[i];
            vDotGrad += 0.5 * (vOld + vNew) * gradOmega[i];
            p.y[9 + i] = vNew;
        }
        p.w += dt * vDotGrad;
    }

    /// One drift-kick-drift step of fictitious length h.
    static void step(const Masses& m, Phase& p, double h) {
        drift(p, 0.5 * h);
        kick(m, p, h);
        drift(p, 0.5 * h);
    }

    /// A single step from `from` that lands exactly on physical time `target`
    /// (the fictitious step length is found by iteration).
    static Phase advanceTo(const Masses& m, const Phase& from, double target) {
        const double interval = target - from.t;
        if (interval <= 0.0) return from;
        double h = interval * from.w;
        Phase trial = from;
        for (int iter = 0; iter < 50; ++iter) {
            trial = from;
            step(m, trial, h);
            if (std::fabs(trial.t - target) <= 1e-10 * interval + 1e-15 * target) break;
            h *= interval / (trial.t - from.t);
        }
        return trial;
    }

    double stepsPerOrbit_;
};

}  // namespace

int main(int argc, char** argv) { return bhs::runMain<TimeTransformedLeapfrog>(argc, argv); }
