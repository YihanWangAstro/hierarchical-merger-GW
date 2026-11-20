// Command-line driver shared by all integrator versions.  Each executable is
//
//     int main(int argc, char** argv) { return bhs::runMain<MyIntegrator>(argc, argv); }
//
// and offers the four stages of the pipeline (code0, codeA, codeB, codeC, or
// `all` for one after the other), the tidal-capture binaries that can stand
// in for Code 0 (tc, all-tc), `lists` to gather what they found for the
// notebooks, and three tools to look at an integrator (encounters, traj,
// profile).
#pragma once

#include <chrono>
#include <fstream>
#include <map>
#include <string>
#include <vector>

#include "code0.hpp"
#include "codeA.hpp"
#include "codeB.hpp"
#include "codeC.hpp"
#include "tc.hpp"

namespace bhs {

struct Options {
    std::string command;
    std::string dir = ".";                // where the stage files are
    std::string inFile;                   // empty: the default input of the command
    std::string outFile;                  // empty: the default output of the command, in `dir`
    double accuracy = 0.0;                // integrator-specific, see Integrator::accuracyHelp()
    unsigned threads = availableThreads();
    Settings settings;                    // what the stages themselves need to know
    int trials0 = 0, trialsA = 0, trialsC = 0;   // settings.trials of one stage only (0: as the others)
    int part = 0, parts = 0;              // a stage: run only part `part` of `parts` (0: all of it)
    int merge = 0;                        // a stage: put so many parts together instead of running
    bool check = false;                   // a stage: only say whether its table is the one it would write

    int numEncounters = 2000;             // encounters: how many
    int encounter = 0;                    // traj: which encounter of the `encounters` list
    int samplesPerOrbit = 500;            // traj: output points per initial binary period
};

template <class Integrator>
void printUsage(const char* program) {
    std::printf(
        "Black holes and stars in a nuclear star cluster -- %s\n"
        "\n"
        "usage: %s <command> [options]\n"
        "\n"
        "the pipeline: each stage reads the file of the one before and writes a\n"
        "table of its results (.csv) and the lists the notebooks plot (.m)\n"
        "  code0        a 2G black hole meets a stellar binary and may capture a star\n"
        "  codeA        the black hole/star binary meets field stars until the star is\n"
        "               disrupted or lost, or an encounter with a black hole comes due\n"
        "  codeB        it meets that 1G black hole: a black-hole binary may form\n"
        "  codeC        the black-hole binary meets stars and black holes until it\n"
        "               merges or the 2G black hole loses its black-hole companion\n"
        "  all          the four stages one after the other, then `lists`\n"
        "  lists        the lists of the four stages in one file, lists_%s.m\n"
        "\n"
        "the same with stars captured by tides instead of from binaries\n"
        "  tc           the black hole/star binary tidal capture leaves for every mass\n"
        "               pair of the input at ten radii, 0.01 to 2 pc: no encounter, an\n"
        "               orbit from two rules (the radii of the input are not read)\n"
        "  all-tc       tc, codeA, codeB, codeC, then `lists`\n"
        "\n"
        "tools\n"
        "  encounters   integrate N independent Code A encounters, one row each\n"
        "  traj         integrate one of those encounters and write its trajectory\n"
        "  profile      write the cluster density / velocity-dispersion table\n"
        "\n"
        "options\n"
        "  --dir DIR        where the stage files code0_%s.csv ... are     [.]\n"
        "  --in FILE        input instead of the previous stage's file in DIR\n"
        "                   (code0, tc: DATA/params0.csv; tools: DATA/bhstar_tidal_capture.csv)\n"
        "  --out FILE       output table instead of <command>_%s.csv in DIR\n"
        "  --data DIR       directory with params0.csv and M1G.txt              [data]\n"
        "  --seed S         random seed                                        [12345]\n"
        "  --acc X          %s [%g]\n"
        "  --threads N      worker threads                   [the cores of the batch job,\n"
        "                                                     or all cores of the machine]\n"
        "  --part I/N       a stage: run only the I-th of N equal parts of its jobs and\n"
        "                   write them to <stage>_%s.partIofN.csv (for job arrays)\n"
        "  --merge N        a stage: put its N parts together into the table and\n"
        "                   write the lists, instead of running anything\n"
        "  --check          a stage: run nothing; say whether the table that is there\n"
        "                   is the one this command would write, from the same input\n"
        "                   and with the same numbers (exit status 0 if it is)\n"
        "  --nreal N        code0: scatterings per stellar binary; codeA, codeC:\n"
        "                   histories per binary                               [100]\n"
        "  --nreal0 N, --nrealA N, --nrealC N    the same for one stage only\n"
        "  --nrealB N       codeB: scatterings per black hole/star binary      [1]\n"
        "  --tidal X        after an encounter a star left bound to a black hole is\n"
        "                   disrupted if its pericentre is within X tidal radii [1]\n"
        "  --nsys N         use only the first N mass pairs of the input       [all]\n"
        "  --n N            encounters: number of encounters                   [2000]\n"
        "  --encounter I    traj: index of the encounter (column n)            [0]\n"
        "  --samples N      traj: output points per initial binary period      [500]\n",
        Integrator::description(), program, Integrator::name(), Integrator::name(), Integrator::name(),
        Integrator::accuracyHelp(), Integrator::defaultAccuracy(), Integrator::name());
}

template <class Integrator>
Options parseOptions(int argc, char** argv) {
    Options opt;
    opt.accuracy = Integrator::defaultAccuracy();
    if (argc < 2 || std::string(argv[1]) == "--help" || std::string(argv[1]) == "-h") {
        printUsage<Integrator>(argv[0]);
        std::exit(argc < 2 ? 1 : 0);
    }
    opt.command = argv[1];
    for (int i = 2; i < argc; ++i) {
        const std::string flag = argv[i];
        auto value = [&]() -> std::string {
            if (i + 1 >= argc) fail("missing value for " + flag);
            return argv[++i];
        };
        if (flag == "--dir") opt.dir = value();
        else if (flag == "--in") opt.inFile = value();
        else if (flag == "--out") opt.outFile = value();
        else if (flag == "--data") opt.settings.dataDir = value();
        else if (flag == "--seed") opt.settings.seed = std::stoull(value());
        else if (flag == "--acc") opt.accuracy = std::stod(value());
        else if (flag == "--threads") opt.threads = static_cast<unsigned>(std::stoul(value()));
        else if (flag == "--part") {
            const std::vector<std::string> ofParts = split(value(), '/');
            if (ofParts.size() != 2) fail("--part takes I/N, as in --part 3/16");
            opt.part = toInt(ofParts[0]), opt.parts = toInt(ofParts[1]);
            if (opt.part < 1 || opt.part > opt.parts) fail("--part I/N needs 1 <= I <= N");
        }
        else if (flag == "--merge") {
            opt.merge = std::stoi(value());
            if (opt.merge < 1) fail("--merge takes the number of parts");
        }
        else if (flag == "--check") opt.check = true;
        else if (flag == "--nreal") opt.settings.trials = std::stoi(value());
        else if (flag == "--nreal0") opt.trials0 = std::stoi(value());
        else if (flag == "--nrealA") opt.trialsA = std::stoi(value());
        else if (flag == "--nrealC") opt.trialsC = std::stoi(value());
        else if (flag == "--nrealB") opt.settings.trialsB = std::stoi(value());
        else if (flag == "--tidal") opt.settings.tidalFactor = std::stod(value());
        else if (flag == "--nsys") opt.settings.maxMassPairs = std::stoi(value());
        else if (flag == "--n") opt.numEncounters = std::stoi(value());
        else if (flag == "--encounter") opt.encounter = std::stoi(value());
        else if (flag == "--samples") opt.samplesPerOrbit = std::stoi(value());
        else if (flag == "--help" || flag == "-h") { printUsage<Integrator>(argv[0]); std::exit(0); }
        else fail("unknown option " + flag + " (try --help)");
    }
    if (opt.threads == 0) opt.threads = 1;
    if (opt.settings.trials < 1) fail("--nreal must be at least 1");
    if (opt.settings.trialsB < 1) fail("--nrealB must be at least 1");
    if (opt.trials0 < 0 || opt.trialsA < 0 || opt.trialsC < 0) fail("--nreal0, --nrealA and --nrealC must be at least 1");
    if (!(opt.settings.tidalFactor > 0.0)) fail("--tidal must be positive: the tidal check is always applied");
    if (opt.part > 0 && opt.merge > 0) fail("--part runs a part, --merge puts the parts together: not both");
    if (opt.check && (opt.part > 0 || opt.merge > 0)) fail("--check looks at the whole table: not with --part or --merge");
    return opt;
}

/// The file a command writes or a later stage reads: <dir>/<command>_<integrator>.csv.
template <class Integrator>
std::string defaultFile(const Options& opt, const std::string& command) {
    return opt.dir + "/" + command + "_" + Integrator::name() + ".csv";
}

/// The path with its extension (if it has one) replaced: results.csv -> results.m.
inline std::string withExtension(const std::string& path, const std::string& extension) {
    const std::size_t dot = path.rfind('.');
    const bool hasExtension = dot != std::string::npos && path.find('/', dot) == std::string::npos;
    return (hasExtension ? path.substr(0, dot) : path) + extension;
}

/// How a file was produced: the start of its first line.
template <class Integrator>
std::string provenance(const Options& opt) {
    return format("integrator=%s accuracy=%g seed=%llu", Integrator::name(), opt.accuracy, opt.settings.seed);
}

// ----------------------------------------------------------------------------
// The pipeline: one function runs any of its stages
// ----------------------------------------------------------------------------

// A stage can be run in parts, by as many jobs of a batch system as there are
// parts (a job array), each on a node of its own: `--part I/N` runs the I-th
// of N equal shares of the jobs and writes them to a file of its own, and
// `--merge N` puts the N files together.  Every job draws its random numbers
// from the seed and its own indices, so the table is the one a single run
// would have written.

/// The file of one part of a table: codeA_spacehub.csv -> codeA_spacehub.part03of16.csv.
inline std::string partFile(const std::string& table, int part, int parts) {
    const int width = static_cast<int>(std::to_string(parts).size());
    return withExtension(table, format(".part%0*dof%d.csv", width, part, parts));
}

/// The last line of a part that ran to its end, before the number of its rows.
inline constexpr const char* kPartEnd = "# end of part, rows: ";

/// True if the file is a part written to its end under this head (its first two lines).
inline bool partComplete(const std::string& path, const std::string& head) {
    std::ifstream in(path);
    std::string start(head.size(), '\0'), line, last;
    if (!in || !in.read(start.data(), static_cast<std::streamsize>(head.size())) || start != head) return false;
    in.seekg(0, std::ios::end);
    const std::streamoff size = in.tellg();
    in.seekg(std::max<std::streamoff>(static_cast<std::streamoff>(head.size()), size - 200));
    while (std::getline(in, line)) last = line;
    return last.rfind(kPartEnd, 0) == 0;
}

/// Take from the first line of a part how the parts were run, so that the
/// lists written when they are put together say the same.
inline void adoptRun(Options& opt, const std::string& part) {
    const CsvReader made(part);
    auto has = [&](const char* key) { return !made.meta(key).empty(); };
    if (has("nreal")) opt.settings.trials = toInt(made.meta("nreal"));
    if (has("nrealB")) opt.settings.trialsB = toInt(made.meta("nrealB"));
    if (has("tidal")) opt.settings.tidalFactor = toNumber(made.meta("tidal"));
    if (has("seed")) opt.settings.seed = std::stoull(made.meta("seed"));
    if (has("accuracy")) opt.accuracy = toNumber(made.meta("accuracy"));
}

/// Put the parts of a stage together into its table, and return the records
/// of all rows.  Every part must have been written to its end, under the
/// same head, and every row must read back as it was written.
template <class Stage>
std::vector<typename Stage::Record> mergeParts(const Stage& stage, const std::string& tablePath,
                                               const std::string& head, int parts) {
    std::vector<typename Stage::Record> records;
    records.reserve(stage.numJobs());
    // Written under another name first: a table that is there already stays
    // as it is if the parts turn out not to be in order.
    const std::string unfinished = tablePath + ".merging";
    {
    OutputFile table(unfinished);
    table.write(head);
    std::string block;   // rows on their way into the table
    for (int part = 1; part <= parts; ++part) {
        const std::string path = partFile(tablePath, part, parts);
        std::ifstream in(path);
        if (!in) fail("cannot open " + path + ": run that part first");
        std::string line, partHead;
        for (int i = 0; i < 2 && std::getline(in, line); ++i) partHead += line + "\n";
        if (partHead != head)
            fail(path + ": not made as part 1 was, or from another input (its first two lines differ)");

        bool complete = false;
        std::size_t rows = 0;
        while (std::getline(in, line)) {
            if (line.rfind(kPartEnd, 0) == 0) {
                complete = std::stoull(line.substr(std::string(kPartEnd).size())) == rows;
                break;
            }
            records.push_back(Stage::parse(split(line, ','), stage.grid()));
            if (Stage::csvRow(records.back()) != line + "\n") fail(path + ": a row does not read back as written");
            block += line + "\n";
            ++rows;
            if (block.size() > (1u << 20)) table.write(block), block.clear();
        }
        if (!complete) fail(path + " was not written to its end: run that part again");
    }
    table.write(block);
    if (records.size() != stage.numJobs())
        fail(format("%s: the parts hold %zu rows, the stage has %zu jobs", tablePath.c_str(), records.size(),
                    stage.numJobs()));
    }
    if (std::rename(unfinished.c_str(), tablePath.c_str()) != 0) fail("cannot write " + tablePath);
    return records;
}

/// The grid a stage works on, without its jobs.  A stage that reads the table
/// of another works on the grid in the first line of that table, however
/// many rows follow; one that starts a pipeline makes the grid of its input.
template <class Stage>
Grid gridOf(const std::string& input, const Settings& settings) {
    if (std::string(Stage::previous()).empty()) return Stage(input, settings).grid();
    return StageReader(input, settings.maxMassPairs).grid();
}

/// --check: is the table at `path` the one whose first two lines are `head`?
/// Says so, and how it was made if it was made otherwise; 0 if it is.
inline int checkTable(const char* stage, const std::string& path, const std::string& head) {
    std::ifstream in(path);
    std::string start(head.size(), '\0');
    if (in && in.read(start.data(), static_cast<std::streamsize>(head.size())) && start == head) {
        std::printf("%s: %s is the table this command would write\n", stage, path.c_str());
        return 0;
    }
    in.clear();
    in.seekg(0);
    std::string first;
    if (!in || !std::getline(in, first)) {
        std::printf("%s: there is no table %s\n", stage, path.c_str());
        return 1;
    }
    // A first line is "# how it was made M12S=... Routs=...": say which of the three differ.
    const std::string asked = head.substr(0, head.find('\n'));
    auto how = [](const std::string& line) { return line.substr(0, line.find(" M12S=")); };
    auto value = [](const std::string& line, const std::string& key) {
        const std::size_t at = line.find(" " + key + "=");
        if (at == std::string::npos) return std::string();
        const std::size_t from = at + key.size() + 2, end = line.find(' ', from);
        return line.substr(from, end == std::string::npos ? end : end - from);
    };
    std::printf("%s: %s is not the table this command would write\n", stage, path.c_str());
    if (how(first) != how(asked))
        std::printf("  it was made as  %s\n  asked for is    %s\n", how(first).c_str(), how(asked).c_str());
    if (value(first, "M12S") != value(asked, "M12S")) std::printf("  its mass pairs are not those of the input\n");
    if (value(first, "Routs") != value(asked, "Routs"))
        std::printf("  its radii are   %s\n  asked for are   %s\n", value(first, "Routs").c_str(),
                    value(asked, "Routs").c_str());
    return 1;
}

/// Run every job of a stage on the worker threads, or one part of them, or
/// put the parts together.  The table is written row by row while the run
/// lasts, the notebook's lists when the whole table is there.
template <class Stage, class Integrator>
int runStage(Options opt) {
    using Record = typename Stage::Record;
    const std::string tablePath = opt.outFile.empty() ? defaultFile<Integrator>(opt, Stage::name()) : opt.outFile;
    const std::string name = Stage::name();
    const int own = name == "code0" ? opt.trials0 : name == "codeA" ? opt.trialsA : name == "codeC" ? opt.trialsC : 0;
    if (own > 0) opt.settings.trials = own;
    if (opt.merge > 0) adoptRun(opt, partFile(tablePath, 1, opt.merge));   // the parts know best how they were run

    const std::string previous = Stage::previous();
    const std::string input = !opt.inFile.empty() ? opt.inFile
                            : previous.empty() ? opt.settings.dataDir + "/params0.csv"
                            : defaultFile<Integrator>(opt, previous);
    const std::string made = format("code=%s %s tidal=%g rstar=%g nreal=%d nrealB=%d", Stage::name(),
                                    provenance<Integrator>(opt).c_str(), opt.settings.tidalFactor, kStarRadius,
                                    opt.settings.trials, opt.settings.trialsB);
    auto headOn = [&](const Grid& grid) { return "# " + made + " " + grid.meta() + "\n" + Stage::csvHeader(); };
    if (opt.check) return checkTable(Stage::name(), tablePath, headOn(gridOf<Stage>(input, opt.settings)));
    const Stage stage(input, opt.settings);
    const std::string head = headOn(stage.grid());

    const std::size_t all = stage.numJobs();
    const std::size_t first = opt.part > 0 ? all * (opt.part - 1) / opt.parts : 0;
    const std::size_t last = opt.part > 0 ? all * opt.part / opt.parts : all;
    const std::string path = opt.part > 0 ? partFile(tablePath, opt.part, opt.parts) : tablePath;
    const std::string share = opt.part > 0 ? format(" (part %d of %d)", opt.part, opt.parts) : std::string();
    if (opt.part > 0 && partComplete(path, head)) {
        std::printf("%s%s: %s is there already, written to its end\n", Stage::name(), share.c_str(), path.c_str());
        return 0;
    }

    std::vector<Record> records;
    Progress progress(opt.merge > 0 ? 0 : last - first, Stage::jobName());
    if (opt.merge > 0) {
        records = mergeParts(stage, tablePath, head, opt.merge);
    } else {
        const Integrator integrator(opt.accuracy);
        const Scatter scatter = [&](const EncounterSetup& setup) { return runEncounter(integrator, setup); };
        OutputFile table(path);
        table.write(head);
        records.resize(last - first);
        OrderedWriter rows(table, records.size());
        parallelFor(records.size(), opt.threads, [&](std::size_t i) {
            records[i] = stage.run(first + i, scatter);
            rows.put(i, Stage::csvRow(records[i]));
            progress.tick();
        });
        if (opt.part > 0) table.write(format("%s%zu\n", kPartEnd, records.size()));
    }

    std::string wrote = path;
    if (opt.part == 0) {   // the whole table is there: the lists can be written
        OutputFile lists(withExtension(tablePath, ".m"));
        lists.write("(* " + made + " *)\n" + stage.mathematica(records));
        wrote += " and " + lists.path();
    }

    std::map<std::string, std::size_t> counts;
    for (const Record& record : records) ++counts[Stage::label(record)];
    std::printf("%s%s: %zu %s from %s\n", Stage::name(), share.c_str(), records.size(), Stage::jobName(),
                input.c_str());
    for (const auto& [label, count] : counts) std::printf("  %8zu  %s\n", count, label.c_str());
    std::printf("%s %s (%.0f s)\n", opt.merge > 0 ? "put together" : "wrote", wrote.c_str(), progress.seconds());
    std::fflush(stdout);   // so that a log file shows each stage of `all` as it ends
    return 0;
}

// ----------------------------------------------------------------------------
// lists: what the four stages found, in one file for the notebooks
// ----------------------------------------------------------------------------

/// The lists in the order of the notebook cell they are pasted into, each
/// group under a line that says where it comes from.
struct ListGroup {
    const char* comment;
    std::vector<const char*> names;
};
inline const std::vector<ListGroup>& listGroups() {
    static const std::vector<ListGroup> groups = {
        {"from Code B: where Code C starts", {"M12BHB", "RBHB", "TBHB"}},
        {"from Code C", {"F3SS", "F3BHBSS", "TMBHBSS", "M3GBHBSS", "M2MERGEBHBSS"}},
        {"sample sizes: Code 0 scatterings per mass pair and radius, Code A histories per black hole/star "
         "binary, Code C histories per black-hole binary; NUMB, if there, Code B scatterings per black "
         "hole/star binary", {"Num0", "NUMA", "NUMC", "NUMB"}},
        {"from Code 0: f1IN0", {"f1"}},
        {"from Code B", {"f2"}},
        {"from Code A, for the plots: the mean of TB", {"X1"}},
        {"from Code C, for the plots: {R [pc], mean merger time}", {"dataTM"}},
    };
    return groups;
}

/// Gather the .m files of the four stages in `dir` into one.  Every list is
/// a line "name=value;" there, and is copied as it is.
template <class Integrator>
int runLists(const Options& opt) {
    std::map<std::string, std::string> lists;   // by name
    std::string made;                           // how the stages were run: their first lines
    // Without the lists of Code 0 the stars were captured by tides (tc): then
    // there is no capture experiment, and so no Num0 and no f1.
    const std::string code0Lists = withExtension(defaultFile<Integrator>(opt, Code0::name()), ".m");
    const bool tidalCapture = !std::ifstream(code0Lists) &&
                              std::ifstream(withExtension(defaultFile<Integrator>(opt, TidalCapture::name()), ".m"));
    for (const char* stage : {tidalCapture ? TidalCapture::name() : Code0::name(), CodeA::name(), CodeB::name(),
                              CodeC::name()}) {
        const std::string path = withExtension(defaultFile<Integrator>(opt, stage), ".m");
        std::ifstream in(path);
        if (!in) fail("cannot open " + path + ": run " + stage + " first");
        for (std::string line; std::getline(in, line);) {
            const std::size_t equals = line.find('=');
            if (line.rfind("(*", 0) == 0) made += line + "\n";
            else if (equals != std::string::npos) lists[line.substr(0, equals)] = std::move(line);
        }
    }

    OutputFile out(opt.outFile.empty() ? opt.dir + "/lists_" + Integrator::name() + ".m" : opt.outFile);
    out.write(made);
    if (tidalCapture)
        out.write("(* the stars were captured by tides, not from binaries: there is no Code 0, no Num0 and no f1 *)\n");
    for (const ListGroup& group : listGroups()) {
        std::string found;
        for (const char* name : group.names) {
            const std::string key = name;
            const bool optional = key == "NUMB" || (tidalCapture && (key == "Num0" || key == "f1"));
            if (lists.count(key)) found += lists[key] + "\n";
            else if (!optional) fail("the stage files do not hold the list " + key);
        }
        if (!found.empty()) out.write(std::string("\n(* ") + group.comment + " *)\n" + found);
    }
    std::printf("wrote %s\n", out.path().c_str());
    return 0;
}

/// The whole pipeline.  The BH/star binaries Code A follows come from Code 0
/// (a stellar binary broken up) or, with tidalCapture, from the tc stage.
template <class Integrator>
int runPipeline(Options opt, bool tidalCapture) {
    if (!opt.outFile.empty()) fail("`all` writes its files into one directory: choose it with --dir, not --out");
    if (opt.part > 0 || opt.merge > 0 || opt.check)
        fail("--part, --merge and --check are for one stage at a time, not for `all`");
    if (tidalCapture) {
        runStage<TidalCapture, Integrator>(opt);
        opt.inFile = defaultFile<Integrator>(opt, TidalCapture::name());   // Code A looks for Code 0's file otherwise
    } else {
        runStage<Code0, Integrator>(opt);
        opt.inFile.clear();
    }
    runStage<CodeA, Integrator>(opt);
    opt.inFile.clear();   // from here on every stage reads the file of the one before
    runStage<CodeB, Integrator>(opt);
    runStage<CodeC, Integrator>(opt);
    return runLists<Integrator>(opt);
}

// ----------------------------------------------------------------------------
// profile
// ----------------------------------------------------------------------------
inline int runProfile(const Options& opt) {
    OutputFile out(opt.outFile.empty() ? opt.dir + "/profile.csv" : opt.outFile);
    out.write("# cluster model; densities in AU^-3, dispersions in AU/yr\n");
    out.write("R_pc,n_star,n_bh,sigma_star,sigma_bh\n");
    for (double radius : kGridRadii) {
        const LocalCluster c = clusterAt(radius);
        out.write(format("%.16g,%.16g,%.16g,%.16g,%.16g\n", c.radiusPc(), c.starDensity, c.bhDensity, c.starDispersion,
                         c.bhDispersion));
    }
    std::printf("wrote %s\n", out.path().c_str());
    return 0;
}

// ----------------------------------------------------------------------------
// encounters / traj: a reproducible list of Code A encounters
// ----------------------------------------------------------------------------

/// The BH/star binaries the list cycles through.
inline std::vector<BhStar> listedBinaries(const Options& opt) {
    Grid grid;
    const std::string input = opt.inFile.empty() ? opt.settings.dataDir + "/bhstar_tidal_capture.csv" : opt.inFile;
    const std::vector<BhStar> binaries = readBhStars(input, opt.settings.maxMassPairs, grid);
    if (binaries.empty()) fail(input + ": no black hole/star binaries");
    return binaries;
}

/// The n-th encounter of the list: the random draws depend only on (seed, n),
/// so every integrator is handed exactly the same initial conditions.
struct ListedEncounter {
    BhStar binary;
    EncounterSetup setup;
};

inline ListedEncounter listedEncounter(int n, unsigned long long seed, const std::vector<BhStar>& binaries) {
    const BhStar& binary = binaries[n % binaries.size()];
    const double sigma = clusterAt(binary.cell.radius).starDispersion;
    Rng rng = seededRng(seed, {0xE1C0, n});
    return {binary, drawBinarySingle({blackHole(binary.cell.m1), star(), star()}, binary.a, binary.e, sigma,
                                     ApproachAngles::Zero, rng)};
}

template <class Integrator>
int runEncounters(const Options& opt) {
    const std::vector<BhStar> binaries = listedBinaries(opt);
    const Integrator integrator(opt.accuracy);

    std::vector<std::string> rows(opt.numEncounters);
    Progress progress(opt.numEncounters, "encounters");
    parallelFor(opt.numEncounters, opt.threads, [&](std::size_t n) {
        const ListedEncounter listed = listedEncounter(static_cast<int>(n), opt.settings.seed, binaries);
        const EncounterSetup& setup = listed.setup;
        const Cell& cell = listed.binary.cell;
        const auto start = std::chrono::steady_clock::now();
        const EncounterResult r = runEncounter(integrator, setup);
        const double ms = std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - start).count();
        rows[n] = format("%zu,%d,%d,%.10g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,"
                         "%d,%s,%.17g,%.17g,%.6g,%ld,%d,%d,%.6g\n",
                         n, cell.kk, cell.jj, cell.radius / kParsec, setup.masses[0], setup.binary.semiMajorAxis,
                         setup.binary.eccentricity, setup.vInfinity, setup.pericentre / setup.binary.semiMajorAxis,
                         setup.inclination * 180.0 / kPi, r.time,
                         static_cast<int>(r.final.outcome), outcomeName(r.final.outcome), r.final.aBinary,
                         r.final.eBinary, r.energyError(), r.steps, r.ok ? 1 : 0, r.retries, ms);
        progress.tick();
    });

    OutputFile out(opt.outFile.empty() ? defaultFile<Integrator>(opt, "encounters") : opt.outFile);
    out.write("# " + provenance<Integrator>(opt) + "\n");
    out.write("n,kk,jj,R_pc,M1,a0,e0,v_inf,rp_over_a,inc_deg,t_end,"
              "outcome,outcome_name,a_bin,e_bin,dE,steps,ok,retries,wall_ms\n");
    for (const std::string& row : rows) out.write(row);
    std::printf("wrote %s (%.1f s)\n", out.path().c_str(), progress.seconds());
    return 0;
}

template <class Integrator>
int runTrajectory(const Options& opt) {
    const Integrator integrator(opt.accuracy);
    const ListedEncounter listed = listedEncounter(opt.encounter, opt.settings.seed, listedBinaries(opt));
    const EncounterSetup& setup = listed.setup;
    const double radiusPc = listed.binary.cell.radius / kParsec;
    // The trajectory covers what `encounters` integrates: up to the look that finds the encounter settled.
    const double duration = runEncounter(integrator, setup).time;
    const double orbits = duration / setup.binary.period;
    const int intervals = static_cast<int>(std::min(2.0e5, std::ceil(orbits * opt.samplesPerOrbit)));
    Sampler sampler(duration, intervals);
    const EncounterResult r = runEncounter(integrator, setup, &sampler);

    OutputFile out(opt.outFile.empty() ? defaultFile<Integrator>(opt, "traj") : opt.outFile);
    out.write(format("# %s encounter=%d M1=%.10g a0=%.10g e0=%.10g R_pc=%.6g outcome=%s a_bin=%.10g e_bin=%.10g\n",
                     provenance<Integrator>(opt).c_str(), opt.encounter, setup.masses[0], setup.binary.semiMajorAxis,
                     setup.binary.eccentricity, radiusPc, outcomeName(r.final.outcome), r.final.aBinary,
                     r.final.eBinary));
    out.write("t,x1,y1,z1,x2,y2,z2,x3,y3,z3,vx1,vy1,vz1,vx2,vy2,vz2,vx3,vy3,vz3\n");
    std::string block;
    for (const TrajectoryPoint& p : sampler.points()) {
        block += format("%.17g", p.t);
        for (double value : p.y) block += format(",%.17g", value);
        block += '\n';
    }
    out.write(block);

    std::printf("encounter %d: system %d at %.3g pc, M1 = %.3f, a0 = %.5g AU, e0 = %.5f\n", opt.encounter,
                listed.binary.cell.kk, radiusPc, setup.masses[0], setup.binary.semiMajorAxis,
                setup.binary.eccentricity);
    std::printf("  v_inf = %.4g AU/yr, rp = %.3g a0, inclination = %.1f deg, duration = %.4g yr (%.1f orbits)\n",
                setup.vInfinity, setup.pericentre / setup.binary.semiMajorAxis, setup.inclination * 180.0 / kPi,
                duration, orbits);
    std::printf("  outcome %s, a_bin = %.10g AU, e_bin = %.10g, steps = %ld, |dE/E| = %.2e\n",
                outcomeName(r.final.outcome), r.final.aBinary, r.final.eBinary, r.steps, r.energyError());
    std::printf("wrote %s (%zu points)\n", out.path().c_str(), sampler.points().size());
    return 0;
}

// ----------------------------------------------------------------------------
template <class Integrator>
int runMain(int argc, char** argv) {
    const Options opt = parseOptions<Integrator>(argc, argv);
    if (opt.command == "code0") return runStage<Code0, Integrator>(opt);
    if (opt.command == "codeA") return runStage<CodeA, Integrator>(opt);
    if (opt.command == "codeB") return runStage<CodeB, Integrator>(opt);
    if (opt.command == "codeC") return runStage<CodeC, Integrator>(opt);
    if (opt.command == "all") return runPipeline<Integrator>(opt, false);
    if (opt.command == "tc") return runStage<TidalCapture, Integrator>(opt);
    if (opt.command == "all-tc") return runPipeline<Integrator>(opt, true);
    if (opt.command == "lists") return runLists<Integrator>(opt);
    if (opt.command == "encounters") return runEncounters<Integrator>(opt);
    if (opt.command == "traj") return runTrajectory<Integrator>(opt);
    if (opt.command == "profile") return runProfile(opt);
    fail("unknown command '" + opt.command + "' (try --help)");
}

}  // namespace bhs
