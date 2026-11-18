// What the four stages of the pipeline share: the (mass pair, radius) grid
// they work on, the files they hand to each other, the run settings, and
// reproducible random numbers for every job.
//
// A stage is a class providing
//
//     using Record = ...;                        what one job produces: a row of its table
//     static const char* name();                 "code0": its command and the stem of its files
//     static const char* previous();             the stage whose table it reads ("" for the first)
//     static const char* jobName();              "scatterings", for the progress line
//     Stage(const std::string& input, const Settings&);                reads the input table
//     const Grid& grid() const;
//     std::size_t numJobs() const;
//     Record run(std::size_t job, const Scatter&) const;               any thread, any order
//     static std::string label(const Record&);                         how it ended, for the summary
//     static std::string csvHeader();
//     static std::string csvRow(const Record&);
//     static Record parse(const std::vector<std::string>& fields, const Grid&);   the row read back
//     std::string mathematica(const std::vector<Record>&) const;       the notebook's lists
//
// and runStage() in driver.hpp does the rest.
#pragma once

#include <array>
#include <initializer_list>
#include <string>
#include <vector>

#include "encounter.hpp"
#include "io.hpp"
#include "mathematica.hpp"

namespace bhs {

// ----------------------------------------------------------------------------
// The grid
// ----------------------------------------------------------------------------

/// One point of the grid, numbered from 1 as in the notebooks.
struct Cell {
    int kk = 0;            // mass pair
    int jj = 0;            // radius
    double m1 = 0.0;       // the 2G black hole followed through the pipeline
    double mbh = 0.0;      // the 1G black hole it pairs with in Code B
    double radius = 0.0;   // distance from the central black hole [AU]
};

/// Every mass pair at every radius.  Each stage file carries the grid in its
/// first line, so that cells in which nothing happened are not forgotten.
struct Grid {
    std::vector<std::array<double, 2>> masses;   // per mass pair: {m1, mbh}
    std::vector<double> radii;                   // [AU]

    int numMass() const { return static_cast<int>(masses.size()); }
    int numRadii() const { return static_cast<int>(radii.size()); }
    int numCells() const { return numMass() * numRadii(); }
    bool contains(int kk, int jj) const { return kk >= 1 && kk <= numMass() && jj >= 1 && jj <= numRadii(); }

    Cell cell(int kk, int jj) const { return {kk, jj, masses[kk - 1][0], masses[kk - 1][1], radii[jj - 1]}; }
    int index(const Cell& c) const { return (c.kk - 1) * numRadii() + c.jj - 1; }   // cells in notebook order

    /// Keep only the first `count` mass pairs (0: all of them).
    void keepMassPairs(int count) {
        if (count > 0 && count < numMass()) masses.resize(count);
    }

    // The grid as Mathematica lists: for the notebooks, and (with plain numbers) for the first line of a stage file.
    std::string massList(std::string (*number)(double) = mma::real) const {
        return mma::list(masses, [&](const auto& m) { return "{" + number(m[0]) + "," + number(m[1]) + "}"; });
    }
    std::string radiusList(std::string (*number)(double) = mma::real) const { return mma::list(radii, number); }
    std::string meta() const { return "M12S=" + massList(num) + " Routs=" + radiusList(num); }
};

/// The numbers in a Mathematica list such as "{{1.5,2},{3,4}}".
inline std::vector<double> numbersIn(std::string list) {
    for (char& c : list)
        if (c == '{' || c == '}' || c == ',') c = ' ';
    std::vector<double> numbers;
    for (const char* at = list.c_str();;) {
        char* end;
        const double value = std::strtod(at, &end);
        if (end == at) return numbers;
        numbers.push_back(value);
        at = end;
    }
}

// ----------------------------------------------------------------------------
// Stage files
// ----------------------------------------------------------------------------

/// Reads the table of a stage: its grid, then the rows that lie on it.
class StageReader : public CsvReader {
   public:
    /// maxMassPairs > 0 cuts the grid down to its first mass pairs.
    StageReader(const std::string& path, int maxMassPairs) : CsvReader(path), kk_(column("kk")), jj_(column("jj")) {
        const std::vector<double> masses = numbersIn(meta("M12S"));
        for (std::size_t i = 0; i + 1 < masses.size(); i += 2) grid_.masses.push_back({masses[i], masses[i + 1]});
        grid_.radii = numbersIn(meta("Routs"));
        if (grid_.numCells() == 0) fail(path + ": the first line does not give the grid (M12S=... Routs=...)");
        grid_.keepMassPairs(maxMassPairs);
    }

    const Grid& grid() const { return grid_; }

    /// Move on to the next row of a cell of the grid; false at the end.
    bool next() {
        while (CsvReader::next())
            if (grid_.contains(integer(kk_), integer(jj_))) return true;
        return false;
    }

    /// The cell of the current row.
    Cell cell() const { return grid_.cell(integer(kk_), integer(jj_)); }

   private:
    Grid grid_;
    int kk_, jj_;
};

/// The cell columns every stage table has after its index columns.
inline const char* cellColumns() { return "R_pc,M1,MBH"; }
inline std::string cellFields(const Cell& c) { return num(c.radius / kParsec) + "," + num(c.m1) + "," + num(c.mbh); }

// ----------------------------------------------------------------------------
// From a flat list of records to the notebooks' nested lists
// ----------------------------------------------------------------------------

/// The records of each cell, in the order given: result[grid.index(cell)].
/// With `keep`, only the records it accepts.
template <class Record, class Keep>
std::vector<std::vector<const Record*>> byCell(const Grid& grid, const std::vector<Record>& records, Keep keep) {
    std::vector<std::vector<const Record*>> cells(grid.numCells());
    for (const Record& record : records)
        if (keep(record)) cells[grid.index(record.cell)].push_back(&record);
    return cells;
}
template <class Record>
std::vector<std::vector<const Record*>> byCell(const Grid& grid, const std::vector<Record>& records) {
    return byCell(grid, records, [](const Record&) { return true; });
}

/// {{f(1,1),f(1,2),...},{f(2,1),...},...}: one entry per cell, nested by mass
/// pair as in the notebooks.  f takes the index of the cell.
template <class F>
std::string perCell(const Grid& grid, F f) {
    return mma::table(grid.numMass(), [&](int k) {
        return mma::table(grid.numRadii(), [&](int j) { return f(k * grid.numRadii() + j); });
    });
}

/// The same with a list per cell: value(record) of every record in cells[index].
template <class Cells, class Value>
std::string listPerCell(const Grid& grid, const Cells& cells, Value value) {
    return perCell(grid, [&](int c) { return mma::list(cells[c], value); });
}

/// part / whole, and zero if there is nothing to take a fraction of.
inline double fraction(std::size_t part, std::size_t whole) {
    return whole == 0 ? 0.0 : static_cast<double>(part) / static_cast<double>(whole);
}

// ----------------------------------------------------------------------------
// Settings and random numbers
// ----------------------------------------------------------------------------

/// What the user can set for a run of the pipeline.
struct Settings {
    unsigned long long seed = 12345;
    int trials = 100;           // experiments per input system (Code 0, A and C)
    int trialsB = 1;            // Code B scatterings per BH/star binary
    double tidalFactor = 1.0;   // a star is disrupted within this many tidal radii of a black hole
    int maxMassPairs = 0;       // use only the first so many mass pairs; 0: all
    std::string dataDir = "data";
};

/// A generator that depends only on the seed and on the key of the job, so
/// that results do not depend on the number of threads or the order of work.
inline Rng seededRng(unsigned long long seed, std::initializer_list<int> key) {
    std::vector<unsigned long long> words{seed};
    words.insert(words.end(), key.begin(), key.end());
    std::seed_seq sequence(words.begin(), words.end());
    return Rng(sequence);
}

}  // namespace bhs
