// Code 0 -- a 2G black hole meets a binary of two stars.  The encounter may
// split the binary and leave one star bound to the black hole: a BH/star
// binary, the input of Code A.  (Notebook "2 Code 0".)
#pragma once

#include "pipeline.hpp"

namespace bhs {

/// A stellar binary and the black hole it is about to meet: the notebook's
/// Params0[[kk, jj]].
struct StellarBinary {
    Cell cell;
    double a = 0.0, e = 0.0;       // orbit of the two stars
    double breakupRadius = 0.0;    // where the black hole's tide splits the binary ("Rps") [AU]
    double formationTime = 0.0;    // analytic encounter timescale "TA" [yr]; handed on, not used here
};

/// One scattering experiment: a row of the notebook's Details0.
struct Capture {
    Cell cell;
    int ii = 0;                                  // experiment number within the cell
    std::string outcome = "UNRESOLVED / RESONANT";
    double a = 0.0, e = 0.0;                     // the binary left behind
    double vInfinity = 0.0, pericentre = 0.0;    // of the incoming stellar binary
    double duration = 0.0;                       // integrated time [yr]
    double energyError = 0.0;

    bool formedBhStar() const { return outcome == "m1/m2-m3" || outcome == "m1/m3-m2"; }
};

/// The notebook's rule for what an encounter left: an isolated binary needs a
/// bound pair whose third body is unbound, moving away, and farther than
/// kEscapeDistance times the larger of the pair's and the original binary's
/// semi-major axis.  Exactly one such pair gives a binary; otherwise the
/// result is "ALL_UNBOUND" if no pair is bound and unresolved if one is.
inline void classifyCapture(const Masses& m, const State& y, double aStars, Capture& out) {
    Vec3 r[3], v[3];
    for (int b = 0; b < 3; ++b) r[b] = position(y, b), v[b] = velocity(y, b);

    int isolatedBinaries = 0;
    bool anyBound = false;
    for (int p = 0; p < 3; ++p) {
        const Split& s = kSplits[p];
        const PairOrbit pair = pairOrbit(r[s.i], v[s.i], m[s.i], r[s.j], v[s.j], m[s.j]);
        if (pair.energy >= 0) continue;
        anyBound = true;
        const OuterOrbit third = outerOrbit(m, r, v, s);
        const double farAway = kEscapeDistance * std::max(pair.semiMajorAxis, aStars);
        if (pair.eccentricity < 1 && third.energy > 0 && dot(third.r, third.v) > 0 && norm(third.r) > farAway) {
            ++isolatedBinaries;
            out.outcome = outcomeName(static_cast<Outcome>(p));   // "m1/m2-m3", "m1/m3-m2" or "m2/m3-m1"
            out.a = pair.semiMajorAxis;
            out.e = pair.eccentricity;
        }
    }
    if (isolatedBinaries != 1) {
        out.outcome = anyBound ? "UNRESOLVED / RESONANT" : "ALL_UNBOUND";
        out.a = out.e = 0.0;
    }
}

/// One encounter of the black hole with the stellar binary.  The pericentre
/// of the binary's centre of mass is drawn within twice the breakup radius.
inline Capture scatterStellarBinary(const Scatter& scatter, const StellarBinary& in, int ii, double tidalFactor,
                                    Rng& rng) {
    constexpr double kLargestEnergyError = 1e-3;

    Capture out;
    out.cell = in.cell;
    out.ii = ii;
    const double sigma = clusterAt(in.cell.radius).starDispersion;
    const double maxPericentre = 2.0 * in.breakupRadius;
    if (!(in.a > 0 && in.e >= 0 && in.e < 1 && maxPericentre > 0 && sigma > 0)) {
        out.outcome = "INVALID_INPUT";
        return out;
    }

    const EncounterSetup setup =
        drawSingleBinary({blackHole(in.cell.m1), star(), star()}, in.a, in.e, maxPericentre, sigma, rng);
    out.vInfinity = setup.vInfinity;
    out.pericentre = setup.pericentre;

    const EncounterResult result = scatter(setup);
    if (!result.ok) {
        out.outcome = "INTEGRATION_FAILED";
        return out;
    }
    out.duration = result.time;
    // Measured against the binding energy of the stellar binary too: the total energy can be close to zero.
    const double energyScale = std::max({std::fabs(result.initialEnergy), kG * kStarMass * kStarMass / in.a, 1e-20});
    out.energyError = std::fabs(result.finalEnergy - result.initialEnergy) / energyScale;
    if (!(out.energyError <= kLargestEnergyError)) {
        out.outcome = "NUMERICAL_ERROR";
        return out;
    }

    classifyCapture(setup.masses, result.finalState, in.a, out);
    if (out.formedBhStar() && tidallyDisrupted(in.cell.m1, out.a, out.e, tidalFactor)) out.outcome = "TIDAL_DISRUPTION";
    return out;
}

/// The stage as the driver runs it.
class Code0 {
   public:
    using Record = Capture;
    static const char* name() { return "code0"; }
    static const char* previous() { return ""; }   // starts the pipeline: reads <data>/params0.csv
    static const char* jobName() { return "scatterings"; }

    /// Reads a table with one row for every mass pair at every radius, in the
    /// columns of the notebook's Params0: kk,jj,M1,MBH,R_AU,a_ss,e_ss,r_tide,t_a.
    /// The rows define the grid of the whole pipeline.
    Code0(const std::string& input, const Settings& settings) : settings_(settings) {
        CsvReader csv(input);
        const int kk = csv.column("kk"), jj = csv.column("jj"), m1 = csv.column("M1"), mbh = csv.column("MBH"),
                  radius = csv.column("R_AU"), a = csv.column("a_ss"), e = csv.column("e_ss"),
                  breakup = csv.column("r_tide"), formation = csv.column("t_a");
        while (csv.next()) {
            const Cell cell{csv.integer(kk), csv.integer(jj), csv.number(m1), csv.number(mbh), csv.number(radius)};
            if (cell.kk < 1 || cell.jj < 1) fail(input + ": kk and jj count from 1");
            if (cell.kk > grid_.numMass()) grid_.masses.resize(cell.kk);
            if (cell.jj > grid_.numRadii()) grid_.radii.resize(cell.jj);
            grid_.masses[cell.kk - 1] = {cell.m1, cell.mbh};
            grid_.radii[cell.jj - 1] = cell.radius;
            inputs_.push_back({cell, csv.number(a), csv.number(e), csv.number(breakup), csv.number(formation)});
        }
        std::sort(inputs_.begin(), inputs_.end(), [&](const StellarBinary& x, const StellarBinary& y) {
            return grid_.index(x.cell) < grid_.index(y.cell);
        });
        bool complete = static_cast<int>(inputs_.size()) == grid_.numCells();
        for (int c = 0; complete && c < grid_.numCells(); ++c) complete = grid_.index(inputs_[c].cell) == c;
        if (!complete) fail(input + ": needs exactly one row for every mass pair at every radius");
        grid_.keepMassPairs(settings.maxMassPairs);
        inputs_.resize(grid_.numCells());
    }

    const Grid& grid() const { return grid_; }
    std::size_t numJobs() const { return inputs_.size() * settings_.trials; }

    Record run(std::size_t job, const Scatter& scatter) const {
        const StellarBinary& in = inputs_[job / settings_.trials];
        const int ii = static_cast<int>(job % settings_.trials) + 1;
        Rng rng = seededRng(settings_.seed, {0, in.cell.kk, in.cell.jj, ii});
        return scatterStellarBinary(scatter, in, ii, settings_.tidalFactor, rng);
    }

    static std::string label(const Record& r) { return r.outcome; }

    static std::string csvHeader() { return format("kk,jj,ii,%s,outcome,a,e,v_inf,rp,t_end,dE\n", cellColumns()); }
    static std::string csvRow(const Record& r) {
        return csvLine({std::to_string(r.cell.kk), std::to_string(r.cell.jj), std::to_string(r.ii),
                        cellFields(r.cell), r.outcome, num(r.a), num(r.e), num(r.vInfinity), num(r.pericentre),
                        num(r.duration), num(r.energyError)});
    }

    /// The record of a row, from the fields csvRow() wrote.
    static Record parse(const std::vector<std::string>& f, const Grid& grid) {
        Record r;
        r.cell = grid.cell(toInt(f[0]), toInt(f[1]));
        r.ii = toInt(f[2]);
        r.outcome = f[6];
        r.a = toNumber(f[7]), r.e = toNumber(f[8]);
        r.vInfinity = toNumber(f[9]), r.pericentre = toNumber(f[10]);
        r.duration = toNumber(f[11]), r.energyError = toNumber(f[12]);
        return r;
    }

    /// For the plots: the fraction of scatterings that leave a BH/star binary
    /// (the notebook's "f1IN0"), and the number of scatterings.
    std::string mathematica(const std::vector<Record>& records) const {
        using namespace mma;
        const auto all = byCell(grid_, records);
        const auto formed = byCell(grid_, records, [](const Record& r) { return r.formedBhStar(); });
        auto formedFraction = [&](int c) { return real(fraction(formed[c].size(), all[c].size())); };
        return assign("Num0", integer(settings_.trials)) + assign("f1", perCell(grid_, formedFraction));
    }

   private:
    Settings settings_;
    Grid grid_;
    std::vector<StellarBinary> inputs_;   // one per cell, in grid order
};

}  // namespace bhs
