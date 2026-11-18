// Tidal capture -- the other way a 2G black hole comes by a star.  A single
// star that passes close enough loses its orbital energy to the tides the
// black hole raises on it, and stays.  Nothing is integrated: the orbit of
// the new BH/star binary follows from the mass of the black hole and the
// velocity dispersion alone.  This stage stands where Code 0 stands, and
// Code A follows the binaries of either.
// (Notebooks "Code A" and "TC Code A": the lists A1sS and E1sS.)
#pragma once

#include <algorithm>

#include "pipeline.hpp"

namespace bhs {

// The notebooks hold the orbits as numbers, 1000 systems at 10 radii.  All
// of them follow the two rules below, to 13 digits:
//
//     a   = G M1 / sigma^2                          the black hole holds the star as tightly
//                                                    as the cluster moves: hard/soft boundary
//     r_p = C (M1 / Msun)^(1/3) sigma^(-2/10.65)    the capture pericentre: weakly larger where
//                                                    the stars are slower
//
// with sigma in AU/yr.  The exponent is that of a tide whose energy falls
// off as r_p^-10.65; C is fitted to the lists.
inline constexpr double kTidalCapturePericentre = 0.01015607651953;   // C [AU]
inline constexpr double kTidalCaptureSlope = -2.0 / 10.65;

struct CaptureOrbit {
    double a, e;
};

inline CaptureOrbit tidalCaptureOrbit(double m1, double dispersion) {
    const double a = kG * m1 / (dispersion * dispersion);
    const double pericentre = kTidalCapturePericentre * std::cbrt(m1) * std::pow(dispersion, kTidalCaptureSlope);
    return {a, 1.0 - pericentre / a};
}

/// The BH/star binary tidal capture leaves in one cell.
struct TidalBinary {
    Cell cell;
    double a = 0.0, e = 0.0;
    double pericentreInTidalRadii = 0.0;   // a (1 - e) over the radius inside which the star is pulled apart
};

/// The stage as the driver runs it: one binary for every mass pair at every
/// radius of its input, a table with the columns kk, jj, M1, MBH and R_AU
/// (the input of Code 0 will do: its other columns are not read).
class TidalCapture {
   public:
    using Record = TidalBinary;
    static const char* name() { return "tc"; }
    static const char* previous() { return ""; }   // starts a pipeline, as Code 0 does
    static const char* jobName() { return "binaries"; }

    TidalCapture(const std::string& input, const Settings& settings) {
        CsvReader csv(input);
        const int kk = csv.column("kk"), jj = csv.column("jj"), m1 = csv.column("M1"), mbh = csv.column("MBH"),
                  radius = csv.column("R_AU");
        while (csv.next()) {
            const Cell cell{csv.integer(kk), csv.integer(jj), csv.number(m1), csv.number(mbh), csv.number(radius)};
            if (cell.kk < 1 || cell.jj < 1) fail(input + ": kk and jj count from 1");
            if (cell.kk > grid_.numMass()) grid_.masses.resize(cell.kk);
            if (cell.jj > grid_.numRadii()) grid_.radii.resize(cell.jj);
            grid_.masses[cell.kk - 1] = {cell.m1, cell.mbh};
            grid_.radii[cell.jj - 1] = cell.radius;
            cells_.push_back(cell);
        }
        std::sort(cells_.begin(), cells_.end(),
                  [&](const Cell& x, const Cell& y) { return grid_.index(x) < grid_.index(y); });
        bool complete = static_cast<int>(cells_.size()) == grid_.numCells();
        for (int c = 0; complete && c < grid_.numCells(); ++c) complete = grid_.index(cells_[c]) == c;
        if (!complete) fail(input + ": needs exactly one row for every mass pair at every radius");
        grid_.keepMassPairs(settings.maxMassPairs);
        cells_.resize(grid_.numCells());
    }

    const Grid& grid() const { return grid_; }
    std::size_t numJobs() const { return cells_.size(); }

    Record run(std::size_t job, const Scatter&) const {
        const Cell& cell = cells_[job];
        const CaptureOrbit orbit = tidalCaptureOrbit(cell.m1, clusterAt(cell.radius).starDispersion);
        if (!(orbit.e >= 0.0 && orbit.e < 1.0))
            fail(format("tidal capture: no orbit for M1 = %g at %g pc: this close to the central black hole the "
                        "capture pericentre lies outside the orbit", cell.m1, cell.radius / kParsec));
        return {cell, orbit.a, orbit.e, orbit.a * (1.0 - orbit.e) / tidalRadius(cell.m1)};
    }

    /// Code A tests the star only after an encounter, as the notebook does: a
    /// binary that starts inside the tidal radius is disrupted there, unless
    /// its first encounter lifts it out.
    static std::string label(const Record& r) {
        return r.pericentreInTidalRadii <= 1.0 ? "pericentre inside the tidal radius" : "pericentre outside the tidal radius";
    }

    static std::string csvHeader() { return format("kk,jj,%s,a,e,rp_rt\n", cellColumns()); }
    static std::string csvRow(const Record& r) {
        return csvLine({std::to_string(r.cell.kk), std::to_string(r.cell.jj), cellFields(r.cell), num(r.a), num(r.e),
                        num(r.pericentreInTidalRadii)});
    }
    static Record parse(const std::vector<std::string>& f, const Grid& grid) {
        return {grid.cell(toInt(f[0]), toInt(f[1])), toNumber(f[5]), toNumber(f[6]), toNumber(f[7])};
    }

    /// The notebook's lists of the orbits: [[kk, jj]].
    std::string mathematica(const std::vector<Record>& records) const {
        using namespace mma;
        return assign("A1sS", perCell(grid_, [&](int c) { return real(records[c].a); })) +
               assign("E1sS", perCell(grid_, [&](int c) { return real(records[c].e); }));
    }

   private:
    Grid grid_;
    std::vector<Cell> cells_;   // in grid order
};

}  // namespace bhs
