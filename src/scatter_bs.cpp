// BH/star binary scattering -- version 1: Bulirsch-Stoer integrator.
//
// Gragg's modified midpoint rule extrapolated to zero step size, with
// adaptive step and order.  General purpose and very accurate; plays the role
// of NDSolve's "ExplicitRungeKutta" with 12-digit goals in the notebook.

#include <algorithm>
#include <cmath>

#include "bhs/driver.hpp"

namespace {

using namespace bhs;

class BulirschStoer {
   public:
    static const char* name() { return "bs"; }
    static const char* description() { return "Bulirsch-Stoer extrapolation"; }
    static double defaultAccuracy() { return 1e-12; }
    static const char* accuracyHelp() { return "relative error tolerance per step                  "; }

    explicit BulirschStoer(double tolerance) : tolerance_(tolerance) {}

    IntegrationStats integrate(const Masses& m, State& y, double tEnd, const BinaryScales& binary,
                               Sampler* sampler) const {
        // Errors are measured relative to |y_i| plus a floor: the pericentre
        // distance for positions, the orbital speed for velocities.
        const double positionFloor = binary.pericentre();
        const double velocityFloor = binary.orbitalSpeed();

        IntegrationStats stats;
        double t = 0.0;
        double h = 1e-3 * binary.period;   // the controller adapts this within a few steps
        State dydt, errorScale, yNew;

        while (t < tEnd) {
            // Do not step over the next output time or the end.
            const double target = (sampler && !sampler->done()) ? std::min(sampler->nextTime(), tEnd) : tEnd;
            bool hitsTarget = t + h >= target;
            double hTry = hitsTarget ? target - t : h;

            derivative(m, y, dydt);
            for (int i = 0; i < 18; ++i)
                errorScale[i] = tolerance_ * (std::fabs(y[i]) + (i < 9 ? positionFloor : velocityFloor));

            Attempt attempt = extrapolate(m, y, dydt, hTry, errorScale, yNew);
            while (!attempt.converged()) {
                hTry *= shrinkFactor(attempt);
                hitsTarget = false;
                if (hTry < 1e-15 * tEnd) return {false, stats.steps};
                attempt = extrapolate(m, y, dydt, hTry, errorScale, yNew);
            }

            y = yNew;
            t = hitsTarget ? target : t + hTry;
            // A step shortened to hit the target says little about the natural step size.
            h = hitsTarget ? std::max(h, hTry * growthFactor(attempt)) : hTry * growthFactor(attempt);

            if (sampler && !sampler->done() && t == sampler->nextTime()) sampler->record(y);
            if (++stats.steps > kMaxSteps) return {false, stats.steps};
        }
        return stats;
    }

   private:
    static constexpr int kColumns = 8;   // extrapolation order up to 2 * kColumns
    static constexpr int kSubsteps[kColumns] = {2, 4, 6, 8, 10, 12, 14, 16};
    static constexpr int kFirstTestedColumn = 3;
    static constexpr long kMaxSteps = 20000000;

    /// Result of one trial step: the column that converged (-1 if none) and its error estimate.
    struct Attempt {
        int column = -1;
        double error = 0.0;   // in units of the tolerance
        bool converged() const { return column >= 0; }
    };

    /// Modified midpoint rule: advance y by H using n substeps.
    static void modifiedMidpoint(const Masses& m, const State& y, const State& dydt, double H, int n, State& out) {
        const double h = H / n;
        State previous = y, current, d;
        for (int i = 0; i < 18; ++i) current[i] = y[i] + h * dydt[i];
        for (int step = 1; step < n; ++step) {
            derivative(m, current, d);
            for (int i = 0; i < 18; ++i) {
                const double next = previous[i] + 2.0 * h * d[i];
                previous[i] = current[i];
                current[i] = next;
            }
        }
        derivative(m, current, d);
        for (int i = 0; i < 18; ++i) out[i] = 0.5 * (previous[i] + current[i] + h * d[i]);
    }

    /// Try a step of size H.  Builds the extrapolation tableau row by row and
    /// stops at the first column whose last correction is below the tolerance.
    Attempt extrapolate(const Masses& m, const State& y, const State& dydt, double H, const State& errorScale,
                        State& yNew) const {
        State previousRow[kColumns], row[kColumns];
        Attempt attempt;
        for (int k = 0; k < kColumns; ++k) {
            modifiedMidpoint(m, y, dydt, H, kSubsteps[k], row[0]);
            // Aitken-Neville extrapolation in h^2.
            for (int j = 1; j <= k; ++j) {
                const double ratio = static_cast<double>(kSubsteps[k]) / kSubsteps[k - j];
                const double weight = 1.0 / (ratio * ratio - 1.0);
                for (int i = 0; i < 18; ++i) row[j][i] = row[j - 1][i] + (row[j - 1][i] - previousRow[j - 1][i]) * weight;
            }
            if (k >= kFirstTestedColumn) {
                attempt.error = 0.0;
                for (int i = 0; i < 18; ++i)
                    attempt.error = std::max(attempt.error, std::fabs(row[k][i] - row[k - 1][i]) / errorScale[i]);
                if (attempt.error < 1.0) {
                    attempt.column = k;
                    yNew = row[k];
                    return attempt;
                }
            }
            for (int j = 0; j <= k; ++j) previousRow[j] = row[j];
        }
        return attempt;
    }

    /// Step-size factor after an accepted step.
    static double growthFactor(const Attempt& attempt) {
        double factor = 0.94 * std::pow(0.65 / std::max(attempt.error, 1e-12), 1.0 / (2 * attempt.column + 1));
        if (attempt.column <= 4) factor = std::max(factor, 1.5);   // converged early: the step was small
        return std::clamp(factor, 0.2, 4.0);
    }

    /// Step-size factor after a rejected step.
    static double shrinkFactor(const Attempt& attempt) {
        return std::clamp(0.94 * std::pow(0.65 / attempt.error, 1.0 / (2 * kColumns - 1)), 0.1, 0.7);
    }

    double tolerance_;
};

}  // namespace

int main(int argc, char** argv) { return bhs::runMain<BulirschStoer>(argc, argv); }
