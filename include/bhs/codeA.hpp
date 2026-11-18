// Code A -- a BH/star binary meets field stars one after another, until it
// is destroyed or an encounter with another black hole comes due.
// (Notebook "Code A"; its survivors are the input of Code B.)
#pragma once

#include "pipeline.hpp"

namespace bhs {

/// A BH/star binary formed in Code 0: the h-th capture of its cell.
struct BhStar {
    Cell cell;
    int h = 0;
    double a = 0.0, e = 0.0;
};

/// How a history ends.  The numbers are the `fate` column of the table; 1 and
/// 2 stood for a merger of black hole and star by gravitational waves, which
/// is a tidal disruption: the star is pulled apart on its way in.
enum class Fate {
    BhEncounter = 0,       // a black-hole encounter comes due: a, e, t are recorded
    CompanionLost = 3,     // the black hole lost its companion (or is the outer body of a triple)
    Ionised = 4,           // all three bodies unbound
    Unresolved = 5,        // an encounter was not finished at the end of its integration
    TidalDisruption = 6    // the star is pulled apart: an encounter left its pericentre inside the tidal radius,
                           // or gravitational waves bring it in before the next star or the black hole arrives
};

inline const char* fateName(Fate fate) {
    switch (fate) {
        case Fate::BhEncounter: return "black-hole encounter reached";
        case Fate::CompanionLost: return "black hole lost its companion";
        case Fate::Ionised: return "ionised";
        case Fate::Unresolved: return "unresolved";
        case Fate::TidalDisruption: return "star tidally disrupted";
    }
    return "";
}

/// One history of one binary.
struct BhStarHistory {
    Cell cell;
    int h = 0, ii = 0;                         // which binary, which history
    double a0 = 0.0, e0 = 0.0;                 // initial orbit
    Fate fate = Fate::Unresolved;
    int numEncounters = 0;                     // stellar encounters integrated
    double a = 0.0, e = 0.0, time = 0.0;       // at the black-hole encounter; zero unless that is the fate
    double aLast = 0.0, eLast = 0.0;           // orbit when the history ended

    bool reachedBh() const { return fate == Fate::BhEncounter; }
};

/// The three clocks of a BH/star binary at a given place in the cluster.
struct Timescales {
    double inspiral;   // gravitational waves bring the star in, to be pulled apart
    double star;       // next encounter with a field star
    double blackHole;  // next encounter with another black hole
};

inline Timescales timescales(const Cell& cell, const LocalCluster& cluster, double a, double e) {
    return {mergerTime(cell.m1, kStarMass, a, e),
            encounterTime(cell.m1 + 2 * kStarMass, cluster.starDensity, cluster.starDispersion, a),
            encounterTime(cell.m1 + cell.mbh + kStarMass, cluster.bhDensity, cluster.bhDispersion, a)};
}

/// Follow the binary until it is destroyed or meets the black hole cell.mbh.
/// The star is taken to be disrupted if an encounter leaves its pericentre
/// within tidalFactor tidal radii of the black hole, and if gravitational
/// waves would bring it in before the next star or the black hole arrives.
inline BhStarHistory evolveBhStar(const Scatter& scatter, const BhStar& in, int ii, double tidalFactor, Rng& rng) {
    const LocalCluster cluster = clusterAt(in.cell.radius);
    const Bodies bodies{blackHole(in.cell.m1), star(), star()};   // black hole, companion, field star

    BhStarHistory out;
    out.cell = in.cell, out.h = in.h, out.ii = ii, out.a0 = in.a, out.e0 = in.e;
    double a = in.a, e = in.e;
    Timescales t = timescales(in.cell, cluster, a, e);

    // Each stellar encounter advances the "black-hole clock" by star/blackHole;
    // a black-hole encounter is due when the clock reaches 1.
    double bhClock = t.star / t.blackHole;
    double elapsed = t.star;

    while (true) {
        // As in the notebook, the incoming orbit's node and argument of pericentre are zero.
        const EncounterSetup setup = drawBinarySingle(bodies, a, e, cluster.starDispersion, ApproachAngles::Zero, rng);
        const FinalConfiguration final = scatter(setup).final;
        ++out.numEncounters;

        if (!final.boundPairIs(0) && !final.boundPairIs(1)) {   // the black hole is not in the tightest pair
            out.fate = final.outcome == Outcome::Ionised ? Fate::Ionised
                     : final.resolved() ? Fate::CompanionLost : Fate::Unresolved;
            break;
        }
        a = final.aBinary;
        e = final.eBinary;
        if (starDisrupted(bodies, final, tidalFactor)) {
            out.fate = Fate::TidalDisruption;
            break;
        }
        t = timescales(in.cell, cluster, a, e);

        const double nextBhClock = bhClock + t.star / t.blackHole;
        const double timeToBh = (1.0 - bhClock) * t.blackHole;

        if (nextBhClock >= 1.0) {   // the black hole arrives before the next star
            if (timeToBh <= t.inspiral) {
                out.fate = Fate::BhEncounter;
                out.a = a;
                out.e = e;
                out.time = elapsed + timeToBh;
            } else {
                out.fate = Fate::TidalDisruption;   // the star spirals in first
            }
            break;
        }
        if (t.inspiral <= t.star) {   // the star spirals in before the next star
            out.fate = Fate::TidalDisruption;
            break;
        }
        bhClock = nextBhClock;      // otherwise: another stellar encounter
        elapsed += t.star;
    }
    out.aLast = a;
    out.eLast = e;
    return out;
}

/// The BH/star binaries in a file: every row of a table with columns
/// kk,jj,a,e, or the captures among the rows of a Code 0 output.
inline std::vector<BhStar> readBhStars(const std::string& path, int maxMassPairs, Grid& grid) {
    StageReader table(path, maxMassPairs);
    grid = table.grid();
    const int a = table.column("a"), e = table.column("e");
    const int outcome = table.has("outcome") ? table.column("outcome") : -1;

    std::vector<BhStar> binaries;
    std::vector<int> found(grid.numCells(), 0);   // per cell
    while (table.next()) {
        if (outcome >= 0 && table.text(outcome) != "m1/m2-m3" && table.text(outcome) != "m1/m3-m2") continue;
        const Cell cell = table.cell();
        binaries.push_back({cell, ++found[grid.index(cell)], table.number(a), table.number(e)});
    }
    return binaries;
}

/// The stage as the driver runs it.
class CodeA {
   public:
    using Record = BhStarHistory;
    static const char* name() { return "codeA"; }
    static const char* previous() { return "code0"; }
    static const char* jobName() { return "histories"; }

    CodeA(const std::string& input, const Settings& settings)
        : settings_(settings), inputs_(readBhStars(input, settings.maxMassPairs, grid_)) {}

    const Grid& grid() const { return grid_; }
    std::size_t numJobs() const { return inputs_.size() * settings_.trials; }

    Record run(std::size_t job, const Scatter& scatter) const {
        const BhStar& in = inputs_[job / settings_.trials];
        const int ii = static_cast<int>(job % settings_.trials) + 1;
        Rng rng = seededRng(settings_.seed, {1, in.cell.kk, in.cell.jj, in.h, ii});
        return evolveBhStar(scatter, in, ii, settings_.tidalFactor, rng);
    }

    static std::string label(const Record& r) { return fateName(r.fate); }

    static std::string csvHeader() {
        return format("kk,jj,h,ii,%s,a0,e0,a,e,tb,fate,n_enc,a_last,e_last\n", cellColumns());
    }
    static std::string csvRow(const Record& r) {
        return csvLine({std::to_string(r.cell.kk), std::to_string(r.cell.jj), std::to_string(r.h),
                        std::to_string(r.ii), cellFields(r.cell), num(r.a0), num(r.e0), num(r.a), num(r.e),
                        num(r.time), std::to_string(static_cast<int>(r.fate)), std::to_string(r.numEncounters),
                        num(r.aLast), num(r.eLast)});
    }

    /// The record of a row, from the fields csvRow() wrote.
    static Record parse(const std::vector<std::string>& f, const Grid& grid) {
        Record r;
        r.cell = grid.cell(toInt(f[0]), toInt(f[1]));
        r.h = toInt(f[2]), r.ii = toInt(f[3]);
        r.a0 = toNumber(f[7]), r.e0 = toNumber(f[8]);
        r.a = toNumber(f[9]), r.e = toNumber(f[10]), r.time = toNumber(f[11]);
        r.fate = static_cast<Fate>(toInt(f[12]));
        r.numEncounters = toInt(f[13]);
        r.aLast = toNumber(f[14]), r.eLast = toNumber(f[15]);
        return r;
    }

    /// For the plots: the mean time a binary needs to reach its black-hole
    /// encounter ("X1"; zero where none did), and the number of histories.
    std::string mathematica(const std::vector<Record>& records) const {
        using namespace mma;
        const auto reached = byCell(grid_, records, [](const Record& r) { return r.reachedBh(); });
        auto meanTime = [&](int c) {
            double sum = 0.0;
            for (const Record* r : reached[c]) sum += r->time;
            return real(reached[c].empty() ? 0.0 : sum / reached[c].size());
        };
        return assign("NUMA", integer(settings_.trials)) + assign("X1", perCell(grid_, meanTime));
    }

   private:
    Settings settings_;
    Grid grid_;
    std::vector<BhStar> inputs_;
};

}  // namespace bhs
