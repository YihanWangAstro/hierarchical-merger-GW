// Code C -- a black-hole binary born in Code B is scattered by stars and by
// 1G black holes, in the order their encounters fall due, until it merges by
// gravitational waves or the 2G black hole loses its black-hole companion.
// (Notebook "2 Code C".)
#pragma once

#include <fstream>

#include "pipeline.hpp"

namespace bhs {

/// A black-hole binary from Code B: the bb-th of its cell.
struct Bhb {
    Cell cell;
    int bb = 0;
    double a = 0.0, e = 0.0;
};

/// One history of one black-hole binary.
struct BhbHistory {
    Cell cell;
    int bb = 0, ii = 0;                      // which binary, which history
    double a0 = 0.0, e0 = 0.0;               // initial orbit
    std::string outcome = "NOT_SET";
    double mergerTime = 0.0;                 // since the binary formed [yr]      } zero unless
    double companionMass = 0.0;              // of the black hole it merged with  } the outcome
    double remnantMass = 0.0;                // of the 3G black hole left         } is "MERGER"
    int numStar = 0, numBh = 0;              // encounters with stars and with black holes
    int numExchange = 0;                     // of them, black holes that took the companion's place

    bool merged() const { return outcome == "MERGER"; }
};

/// The masses of the 1G black holes a binary can meet.
class BhMassSample {
   public:
    explicit BhMassSample(const std::string& path) {
        std::ifstream in(path);
        if (!in) fail("cannot open " + path);
        for (double mass; in >> mass;) {
            masses_.push_back(mass);
            mean_ += mass;
            largest_ = std::max(largest_, mass);
        }
        if (masses_.empty()) fail(path + ": no masses found");
        mean_ /= masses_.size();
    }

    double mean() const { return mean_; }

    /// A mass drawn with probability proportional to (binaryMass + m):
    /// gravitational focusing favours the heavy ones.
    double drawFocused(double binaryMass, Rng& rng) const {
        while (true) {
            const double mass = masses_[uniformIndex(rng, masses_.size())];
            if (uniform(rng, 0.0, binaryMass + largest_) < binaryMass + mass) return mass;
        }
    }

   private:
    std::vector<double> masses_;
    double mean_ = 0.0, largest_ = 0.0;
};

/// Follow the binary of the 2G black hole cell.m1.  Three clocks run: the
/// merger time and the times to the next star and the next black hole.  The
/// clock of the population that did not come keeps the fraction it has run.
inline BhbHistory evolveBhb(const Scatter& scatter, const Bhb& in, int ii, const BhMassSample& fieldBhs,
                            double tidalFactor, Rng& rng) {
    constexpr int kMaxScatterings = 1000;

    const LocalCluster cluster = clusterAt(in.cell.radius);
    const double m1 = in.cell.m1;
    double m2 = in.cell.mbh;   // the companion; replaced when a black hole takes its place
    double a = in.a, e = in.e;

    BhbHistory out;
    out.cell = in.cell, out.bb = in.bb, out.ii = ii, out.a0 = a, out.e0 = e;
    if (!(a > 0 && e >= 0 && e < 1)) {
        out.outcome = "BAD_IC";
        return out;
    }

    double elapsed = 0.0;
    double starClock = 0.0, bhClock = 0.0;   // fraction of the waiting time already run
    while (true) {
        if (out.numStar + out.numBh >= kMaxScatterings) {
            out.outcome = "MAX_SCATTER";
            break;
        }
        const double tMerge = mergerTime(m1, m2, a, e);
        const double tStar = encounterTime(m1 + m2 + kStarMass, cluster.starDensity, cluster.starDispersion, a);
        const double tBh = encounterTime(m1 + m2 + fieldBhs.mean(), cluster.bhDensity, cluster.bhDispersion, a);
        const double toStar = std::max(0.0, (1.0 - starClock) * tStar);
        const double toBh = std::max(0.0, (1.0 - bhClock) * tBh);

        if (tMerge <= toStar && tMerge <= toBh) {
            out.outcome = "MERGER";
            out.mergerTime = elapsed + tMerge;
            out.companionMass = m2;
            out.remnantMass = kRemnantMassFraction * (m1 + m2);
            break;
        }

        // The next encounter: with a star or with a black hole, whichever is due first.
        const bool withStar = toStar <= toBh;
        const Body intruder = withStar ? star() : blackHole(fieldBhs.drawFocused(m1 + m2, rng));
        const double sigma = withStar ? cluster.starDispersion : cluster.bhDispersion;
        const double wait = withStar ? toStar : toBh;
        elapsed += wait;
        if (withStar) {
            starClock = 0.0;
            bhClock = std::min(1.0, bhClock + wait / tBh);
            ++out.numStar;
        } else {
            bhClock = 0.0;
            starClock = std::min(1.0, starClock + wait / tStar);
            ++out.numBh;
        }

        const Bodies bodies{blackHole(m1), blackHole(m2), intruder};
        const EncounterSetup setup = drawBinarySingle(bodies, a, e, sigma, ApproachAngles::Random, rng);
        const EncounterResult result = scatter(setup);
        if (!result.ok) {
            out.outcome = "NDSOLVE_FAILED";   // the notebook's name for an integration that did not finish
            break;
        }
        const FinalConfiguration& final = result.final;
        if (!final.resolved()) {
            out.outcome = withStar ? "STAR_UNRESOLVED_RESONANT" : "BH_UNRESOLVED_RESONANT";
            break;
        }
        if (starDisrupted(bodies, final, tidalFactor)) {   // a black hole took the star as its partner, and destroys it
            out.outcome = "TIDAL_DISRUPTION";
            break;
        }

        // Is the 2G black hole still in a black-hole binary?
        if (final.boundPairIs(0)) {                       // with the same companion
            a = final.aBinary, e = final.eBinary;
        } else if (!withStar && final.boundPairIs(1)) {   // with the newcomer
            m2 = intruder.mass;
            a = final.aBinary, e = final.eBinary;
            ++out.numExchange;
        } else {
            out.outcome = std::string(withStar ? "" : "TARGET_2G_NO_BHB_AFTER_BH:") + outcomeName(final.outcome);
            break;
        }
    }
    return out;
}

/// The stage as the driver runs it.
class CodeC {
   public:
    using Record = BhbHistory;
    static const char* name() { return "codeC"; }
    static const char* previous() { return "codeB"; }
    static const char* jobName() { return "histories"; }

    /// Reads a Code B output and follows every black-hole binary in it.
    CodeC(const std::string& input, const Settings& settings)
        : settings_(settings), fieldBhs_(settings.dataDir + "/M1G.txt") {
        StageReader table(input, settings.maxMassPairs);
        grid_ = table.grid();
        const int outcome = table.column("outcome"), a = table.column("a"), e = table.column("e");
        std::vector<int> found(grid_.numCells(), 0);   // per cell
        while (table.next()) {
            if (table.text(outcome) != "m1/m3-m2" && table.text(outcome) != "(m1/m3)/m2") continue;
            const Cell cell = table.cell();
            inputs_.push_back({cell, ++found[grid_.index(cell)], table.number(a), table.number(e)});
        }
    }

    const Grid& grid() const { return grid_; }
    std::size_t numJobs() const { return inputs_.size() * settings_.trials; }

    Record run(std::size_t job, const Scatter& scatter) const {
        const Bhb& in = inputs_[job / settings_.trials];
        const int ii = static_cast<int>(job % settings_.trials) + 1;
        Rng rng = seededRng(settings_.seed, {3, in.cell.kk, in.cell.jj, in.bb, ii});
        return evolveBhb(scatter, in, ii, fieldBhs_, settings_.tidalFactor, rng);
    }

    static std::string label(const Record& r) { return r.outcome; }

    static std::string csvHeader() {
        return format("kk,jj,bb,ii,%s,a0,e0,outcome,tm,m2_merge,m_3g,n_star,n_bh,n_exch\n", cellColumns());
    }
    static std::string csvRow(const Record& r) {
        return csvLine({std::to_string(r.cell.kk), std::to_string(r.cell.jj), std::to_string(r.bb),
                        std::to_string(r.ii), cellFields(r.cell), num(r.a0), num(r.e0), r.outcome, num(r.mergerTime),
                        num(r.companionMass), num(r.remnantMass), std::to_string(r.numStar), std::to_string(r.numBh),
                        std::to_string(r.numExchange)});
    }

    /// The record of a row, from the fields csvRow() wrote.
    static Record parse(const std::vector<std::string>& f, const Grid& grid) {
        Record r;
        r.cell = grid.cell(toInt(f[0]), toInt(f[1]));
        r.bb = toInt(f[2]), r.ii = toInt(f[3]);
        r.a0 = toNumber(f[7]), r.e0 = toNumber(f[8]);
        r.outcome = f[9];
        r.mergerTime = toNumber(f[10]), r.companionMass = toNumber(f[11]), r.remnantMass = toNumber(f[12]);
        r.numStar = toInt(f[13]), r.numBh = toInt(f[14]), r.numExchange = toInt(f[15]);
        return r;
    }

    /// For the plots.  The lists ending in BHBSS hold one entry per
    /// black-hole binary: [[kk, jj, bb]].
    ///   F3SS          per cell: the mean merger fraction of its binaries
    ///   F3BHBSS       per binary: the fraction of its histories that merge
    ///   TMBHBSS       per binary: the merger times                    } one per
    ///   M3GBHBSS      per binary: the masses of the remnants          } merger
    ///   M2MERGEBHBSS  per binary: the masses of the companions merged with
    ///   dataTM        per mass pair: {R [pc], mean merger time} at the radii with a merger
    std::string mathematica(const std::vector<Record>& records) const {
        using namespace mma;
        using Histories = std::vector<const Record*>;

        // binaries[cell][bb - 1]: the histories of one binary
        std::vector<std::vector<Histories>> binaries(grid_.numCells());
        for (const Record& r : records) {
            std::vector<Histories>& cell = binaries[grid_.index(r.cell)];
            if (static_cast<int>(cell.size()) < r.bb) cell.resize(r.bb);
            cell[r.bb - 1].push_back(&r);
        }
        auto mergers = [](const Histories& histories) {
            Histories merged;
            for (const Record* r : histories)
                if (r->merged()) merged.push_back(r);
            return merged;
        };
        auto mergerFraction = [&](const Histories& h) { return fraction(mergers(h).size(), h.size()); };
        auto perBinary = [&](auto f) { return listPerCell(grid_, binaries, f); };
        auto ofMergers = [&](auto value) {
            return perBinary([&](const Histories& h) { return list(mergers(h), value); });
        };

        auto meanMergerFraction = [&](int c) {
            double sum = 0.0;
            for (const Histories& h : binaries[c]) sum += mergerFraction(h);
            return real(binaries[c].empty() ? 0.0 : sum / binaries[c].size());
        };
        auto meanMergerTimes = [&](int k) {   // of one mass pair, against radius
            std::string points;
            for (int j = 0; j < grid_.numRadii(); ++j) {
                double sum = 0.0;
                std::size_t count = 0;
                for (const Histories& h : binaries[k * grid_.numRadii() + j])
                    for (const Record* r : mergers(h)) sum += r->mergerTime, ++count;
                if (count == 0) continue;
                points += (points.empty() ? "{" : ",{") + real(grid_.radii[j] / kParsec) + "," + real(sum / count) + "}";
            }
            return "{" + points + "}";
        };

        return assign("F3SS", perCell(grid_, meanMergerFraction)) +
               assign("F3BHBSS", perBinary([&](const Histories& h) { return real(mergerFraction(h)); })) +
               assign("TMBHBSS", ofMergers([](const Record* r) { return real(r->mergerTime); })) +
               assign("M3GBHBSS", ofMergers([](const Record* r) { return real(r->remnantMass); })) +
               assign("M2MERGEBHBSS", ofMergers([](const Record* r) { return real(r->companionMass); })) +
               assign("NUMC", integer(settings_.trials)) + assign("dataTM", table(grid_.numMass(), meanMergerTimes));
    }

   private:
    Settings settings_;
    BhMassSample fieldBhs_;
    Grid grid_;
    std::vector<Bhb> inputs_;
};

}  // namespace bhs
