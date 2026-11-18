// Reading and writing CSV files, and running jobs on several threads.
#pragma once

#include <atomic>
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <initializer_list>
#include <iostream>
#include <map>
#include <mutex>
#include <string>
#include <thread>
#include <vector>

#include <unistd.h>   // isatty
#ifdef __linux__
#include <sched.h>    // sched_getaffinity
#endif

namespace bhs {

[[noreturn]] inline void fail(const std::string& message) {
    std::cerr << "error: " << message << "\n";
    std::exit(1);
}

/// printf into a std::string.
template <class... Args>
std::string format(const char* fmt, Args... args) {
    const int size = std::snprintf(nullptr, 0, fmt, args...);
    std::string out(size, '\0');
    std::snprintf(out.data(), size + 1, fmt, args...);
    return out;
}

/// The shortest decimal text that reads back as exactly x.
inline std::string num(double x) {
    for (int digits = 15; digits < 17; ++digits) {
        const std::string text = format("%.*g", digits, x);
        if (std::strtod(text.c_str(), nullptr) == x) return text;
    }
    return format("%.17g", x);
}

/// One line of a CSV file.
inline std::string csvLine(std::initializer_list<std::string> fields) {
    std::string line;
    for (const std::string& field : fields) line += (line.empty() ? "" : ",") + field;
    return line + "\n";
}

inline double toNumber(const std::string& text) { return std::strtod(text.c_str(), nullptr); }
inline int toInt(const std::string& text) { return static_cast<int>(std::strtol(text.c_str(), nullptr, 10)); }

inline std::vector<std::string> split(const std::string& text, char separator) {
    std::vector<std::string> parts;
    std::size_t start = 0;
    for (std::size_t end; (end = text.find(separator, start)) != std::string::npos; start = end + 1)
        parts.push_back(text.substr(start, end - start));
    parts.push_back(text.substr(start));
    return parts;
}

/// Reads, row by row, a CSV file in the layout this code writes: an optional
/// first line "# key=value key=value ...", a line of column names, the rows.
class CsvReader {
   public:
    explicit CsvReader(const std::string& path) : path_(path), in_(path) {
        if (!in_) fail("cannot open " + path);
        std::string line;
        while (readLine(line) && (line.empty() || line[0] == '#'))
            for (const std::string& word : split(line, ' '))
                if (const std::size_t eq = word.find('='); eq != std::string::npos)
                    meta_[word.substr(0, eq)] = word.substr(eq + 1);
        columns_ = split(line, ',');
    }

    const std::string& path() const { return path_; }

    /// Value of a key of the first line; empty if it is not there.
    std::string meta(const std::string& key) const {
        const auto found = meta_.find(key);
        return found == meta_.end() ? std::string() : found->second;
    }

    bool has(const std::string& column) const { return find(column) >= 0; }

    /// Position of a column, to be handed to text() and number().
    int column(const std::string& name) const {
        const int index = find(name);
        if (index < 0) fail(path_ + ": no column '" + name + "'");
        return index;
    }

    /// Move on to the next row; false at the end of the file.
    bool next() {
        std::string line;
        do {
            if (!readLine(line)) return false;
        } while (line.empty() || line[0] == '#');
        fields_ = split(line, ',');
        if (fields_.size() != columns_.size())
            fail(format("%s: a row has %zu fields, the header %zu", path_.c_str(), fields_.size(), columns_.size()));
        return true;
    }

    const std::string& text(int column) const { return fields_[column]; }
    double number(int column) const { return std::strtod(fields_[column].c_str(), nullptr); }
    int integer(int column) const { return static_cast<int>(std::strtol(fields_[column].c_str(), nullptr, 10)); }

   private:
    bool readLine(std::string& line) {
        if (!std::getline(in_, line)) return false;
        if (!line.empty() && line.back() == '\r') line.pop_back();
        return true;
    }
    int find(const std::string& name) const {
        for (std::size_t i = 0; i < columns_.size(); ++i)
            if (columns_[i] == name) return static_cast<int>(i);
        return -1;
    }

    std::string path_;
    std::ifstream in_;
    std::map<std::string, std::string> meta_;
    std::vector<std::string> columns_, fields_;
};

/// Output file that several threads may append whole blocks to.
class OutputFile {
   public:
    explicit OutputFile(const std::string& path) : path_(path), file_(std::fopen(path.c_str(), "w")) {
        if (!file_) fail("cannot write " + path);
    }
    ~OutputFile() { std::fclose(file_); }
    OutputFile(const OutputFile&) = delete;
    OutputFile& operator=(const OutputFile&) = delete;

    void write(const std::string& text) {
        std::lock_guard<std::mutex> lock(mutex_);
        std::fputs(text.c_str(), file_);
        std::fflush(file_);
    }
    const std::string& path() const { return path_; }

   private:
    std::string path_;
    std::FILE* file_;
    std::mutex mutex_;
};

/// Merges blocks that worker threads finish in any order into one file, in
/// index order: block i is written as soon as blocks 0..i-1 have arrived.
class OrderedWriter {
   public:
    OrderedWriter(OutputFile& out, std::size_t count) : out_(out), blocks_(count), arrived_(count, false) {}

    void put(std::size_t index, std::string block) {
        std::lock_guard<std::mutex> lock(mutex_);
        blocks_[index] = std::move(block);
        arrived_[index] = true;
        std::string ready;
        for (; next_ < blocks_.size() && arrived_[next_]; ++next_) {
            ready += blocks_[next_];
            std::string().swap(blocks_[next_]);   // release the memory
        }
        if (!ready.empty()) out_.write(ready);
    }

   private:
    OutputFile& out_;
    std::vector<std::string> blocks_;
    std::vector<bool> arrived_;
    std::size_t next_ = 0;
    std::mutex mutex_;
};

/// Progress on stderr; tick() may be called from any thread.  On a terminal
/// one line is rewritten twice a second; in the log file of a batch job a
/// new line is written once a minute.
class Progress {
   public:
    Progress(std::size_t total, std::string what)
        : total_(total), what_(std::move(what)), onTerminal_(isatty(fileno(stderr)) != 0),
          interval_(onTerminal_ ? 0.5 : 60.0), lastPrint_(-interval_) {}

    void tick() {
        std::lock_guard<std::mutex> lock(mutex_);
        ++done_;
        const double elapsed = seconds();
        if (done_ < total_ && elapsed < lastPrint_ + interval_) return;
        lastPrint_ = elapsed;
        std::fprintf(stderr, onTerminal_ ? "\r%zu / %zu %s, %.0f s elapsed, ~%.0f s left   "
                                         : "%zu / %zu %s, %.0f s elapsed, ~%.0f s left\n",
                     done_, total_, what_.c_str(), elapsed, elapsed / done_ * (total_ - done_));
        if (onTerminal_ && done_ == total_) std::fprintf(stderr, "\n");
    }
    double seconds() const { return std::chrono::duration<double>(Clock::now() - start_).count(); }

   private:
    using Clock = std::chrono::steady_clock;
    std::size_t total_, done_ = 0;
    std::string what_;
    bool onTerminal_;
    double interval_;     // seconds between two progress lines
    double lastPrint_;
    std::mutex mutex_;
    Clock::time_point start_ = Clock::now();
};

/// The number of cores this process may use: what a batch system granted
/// the job (SLURM_CPUS_PER_TASK, or the cores the process is bound to), and
/// otherwise all cores of the machine.  A compute node has more cores than
/// a job is given, and using them all would slow down everyone on it.
inline unsigned availableThreads() {
    if (const char* granted = std::getenv("SLURM_CPUS_PER_TASK")) {
        const long n = std::strtol(granted, nullptr, 10);
        if (n > 0) return static_cast<unsigned>(n);
    }
#ifdef __linux__
    cpu_set_t bound;
    if (sched_getaffinity(0, sizeof(bound), &bound) == 0 && CPU_COUNT(&bound) > 0)
        return static_cast<unsigned>(CPU_COUNT(&bound));
#endif
    const unsigned n = std::thread::hardware_concurrency();
    return n > 0 ? n : 1;
}

/// Call task(i) for every i in [0, count) using `threads` worker threads.
template <class Task>
void parallelFor(std::size_t count, unsigned threads, Task task) {
    std::atomic<std::size_t> next{0};
    auto worker = [&] {
        for (std::size_t i; (i = next.fetch_add(1)) < count;) task(i);
    };
    std::vector<std::thread> pool;
    for (unsigned t = 1; t < threads; ++t) pool.emplace_back(worker);
    worker();
    for (std::thread& thread : pool) thread.join();
}

}  // namespace bhs
