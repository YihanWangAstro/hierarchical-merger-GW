// Code B -- a BH/star binary that survived Code A meets a 1G black hole.  If
// the newcomer takes the star's place, a black-hole binary is born: the input
// of Code C.  (Notebook "2 Code B".)
#pragma once

#include "codeA.hpp"
#include "pipeline.hpp"

namespace bhs {

/// A BH/star binary at the moment its black-hole encounter comes due: the
/// n-th survivor of its cell, history ii of the Code 0 capture h.
struct BhStarAtBh {
    Cell cell;
    int n = 0, h = 0, ii = 0;
    double a = 0.0, e = 0.0;
    double time = 0.0;   // elapsed in Code A [yr]
};

/// One BH/star - BH scattering.
struct BhbFormation {
    Cell cell;
    int n = 0, h = 0, ii = 0;                  // as in BhStarAtBh
    int rep = 1;                               // which scattering of this binary
    double a0 = 0.0, e0 = 0.0, time = 0.0;     // the incoming binary
    std::string outcome = "UNRESOLVED / RESONANT";
    double a = 0.0, e = 0.0;                   // tightest bound pair after the encounter
    double aOuter = 0.0, eOuter = 0.0;         // third body's orbit (bound triples)
    double mutualInclination = 0.0;            // [deg] (bound triples)
    double recoilSpeed = 0.0;                  // of the pair, if the third body escaped [AU/yr]

    /// The two black holes are bound to each other, the star having escaped or orbiting them.
    bool formedBhb() const { return isolated() || outcome == "(m1/m3)/m2"; }
    bool isolated() const { return outcome == "m1/m3-m2"; }
};

/// One encounter of the binary (black hole cell.m1 + star) with the black hole cell.mbh.
inline BhbFormation scatterBhStar(const Scatter& scatter, const BhStarAtBh& in, double tidalFactor, Rng& rng) {
    const LocalCluster cluster = clusterAt(in.cell.radius);
    const Bodies bodies{blackHole(in.cell.m1), star(), blackHole(in.cell.mbh)};

    BhbFormation out;
    out.cell = in.cell, out.n = in.n, out.h = in.h, out.ii = in.ii;
    out.a0 = in.a, out.e0 = in.e, out.time = in.time;

    const EncounterSetup setup =
        drawBinarySingle(bodies, in.a, in.e, cluster.bhDispersion, ApproachAngles::Random, rng);
    const EncounterResult result = scatter(setup);
    if (!result.ok) {
        out.outcome = "NDSOLVE_FAILED";   // the notebook's name for an integration that did not finish
        return out;
    }
    const FinalConfiguration& final = result.final;
    out.outcome = starDisrupted(bodies, final, tidalFactor) ? "TIDAL_DISRUPTION" : outcomeName(final.outcome);
    out.a = final.aBinary, out.e = final.eBinary;
    out.aOuter = final.aOuter, out.eOuter = final.eOuter, out.mutualInclination = final.mutualInclination;
    out.recoilSpeed = final.recoilSpeed;
    return out;
}

/// The stage as the driver runs it.
class CodeB {
   public:
    using Record = BhbFormation;
    static const char* name() { return "codeB"; }
    static const char* previous() { return "codeA"; }
    static const char* jobName() { return "scatterings"; }

    /// Reads a Code A output: the histories that reached the black-hole
    /// encounter are scattered once each; the others only count.
    CodeB(const std::string& input, const Settings& settings) : settings_(settings) {
        StageReader table(input, settings.maxMassPairs);
        grid_ = table.grid();
        histories_.assign(grid_.numCells(), 0);
        const int h = table.column("h"), ii = table.column("ii"), fate = table.column("fate"), a = table.column("a"),
                  e = table.column("e"), time = table.column("tb");
        std::vector<int> found(grid_.numCells(), 0);   // per cell
        while (table.next()) {
            const Cell cell = table.cell();
            ++histories_[grid_.index(cell)];
            if (table.integer(fate) != static_cast<int>(Fate::BhEncounter)) continue;
            inputs_.push_back({cell, ++found[grid_.index(cell)], table.integer(h), table.integer(ii), table.number(a),
                               table.number(e), table.number(time)});
        }
    }

    const Grid& grid() const { return grid_; }
    std::size_t numJobs() const { return inputs_.size() * settings_.trialsB; }

    /// Every binary is scattered settings.trialsB times, each time by a
    /// black hole on another orbit.  The stage that forms the fewest systems
    /// is also the cheapest to repeat: one encounter, against the whole
    /// history that led to it.
    Record run(std::size_t job, const Scatter& scatter) const {
        const BhStarAtBh& in = inputs_[job / settings_.trialsB];
        const int rep = static_cast<int>(job % settings_.trialsB) + 1;
        // The first scattering draws what the only one drew before there were several.
        Rng rng = rep == 1 ? seededRng(settings_.seed, {2, in.cell.kk, in.cell.jj, in.n})
                           : seededRng(settings_.seed, {2, in.cell.kk, in.cell.jj, in.n, rep});
        Record out = scatterBhStar(scatter, in, settings_.tidalFactor, rng);
        out.rep = rep;
        return out;
    }

    static std::string label(const Record& r) { return r.outcome; }

    static std::string csvHeader() {
        return format("kk,jj,n,rep,h,ii,%s,a0,e0,tb,outcome,a,e,a_out,e_out,i_mut,v_recoil\n", cellColumns());
    }
    static std::string csvRow(const Record& r) {
        return csvLine({std::to_string(r.cell.kk), std::to_string(r.cell.jj), std::to_string(r.n),
                        std::to_string(r.rep), std::to_string(r.h), std::to_string(r.ii), cellFields(r.cell), num(r.a0),
                        num(r.e0),
                        num(r.time), r.outcome, num(r.a), num(r.e), num(r.aOuter), num(r.eOuter),
                        num(r.mutualInclination), num(r.recoilSpeed)});
    }

    /// The record of a row, from the fields csvRow() wrote.
    static Record parse(const std::vector<std::string>& f, const Grid& grid) {
        Record r;
        r.cell = grid.cell(toInt(f[0]), toInt(f[1]));
        r.n = toInt(f[2]), r.rep = toInt(f[3]), r.h = toInt(f[4]), r.ii = toInt(f[5]);
        r.a0 = toNumber(f[9]), r.e0 = toNumber(f[10]), r.time = toNumber(f[11]);
        r.outcome = f[12];
        r.a = toNumber(f[13]), r.e = toNumber(f[14]);
        r.aOuter = toNumber(f[15]), r.eOuter = toNumber(f[16]);
        r.mutualInclination = toNumber(f[17]), r.recoilSpeed = toNumber(f[18]);
        return r;
    }

    /// For the plots: the grid and the time each new black-hole binary took
    /// to form (the start of Code C), and the probability "f2" that a BH/star
    /// binary of Code A ends as one: it survives Code A and pairs up here.
    std::string mathematica(const std::vector<Record>& records) const {
        using namespace mma;
        const auto all = byCell(grid_, records);
        const auto formed = byCell(grid_, records, [](const Record& r) { return r.formedBhb(); });
        auto time = [](const Record* r) { return real(r->time); };
        auto survival = [&](int c) { return fraction(all[c].size() / settings_.trialsB, histories_[c]); };   // of Code A
        auto formation = [&](int c) { return fraction(formed[c].size(), all[c].size()); };                   // of Code B
        // The notebooks scatter every binary once and have no name for the number: written only if it is not 1.
        const std::string repeats = settings_.trialsB > 1 ? assign("NUMB", integer(settings_.trialsB)) : "";

        return repeats + assign("M12BHB", grid_.massList()) +
               assign("RBHB", table(grid_.numMass(), [&](int) { return grid_.radiusList(); })) +
               assign("TBHB", listPerCell(grid_, formed, time)) +
               assign("f2", perCell(grid_, [&](int c) { return real(survival(c) * formation(c)); }));
    }

   private:
    Settings settings_;
    Grid grid_;
    std::vector<BhStarAtBh> inputs_;
    std::vector<std::size_t> histories_;   // per cell: Code A histories, whatever their fate
};

}  // namespace bhs
