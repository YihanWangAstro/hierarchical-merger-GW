# Black holes and stars in a nuclear star cluster

> This is the public copy of the code. The Mathematica notebooks the stages
> were written from, their text copies and the documents under `docs/` are
> not part of it, so the links to them below lead nowhere.

C++ version of four Mathematica notebooks that follow a second-generation
(2G) black hole through a stellar cusp around a 4×10⁶ M☉ black hole, from
the capture of a star to a merger with another black hole. Every encounter is
integrated as a Newtonian three-body problem.

| stage | notebook | what happens | it ends with |
|---|---|---|---|
| `code0` | `2 Code 0.nb` | the 2G black hole meets a binary of two stars | a BH/star binary, if a star stays bound |
| `codeA` | `Code A.nb` | the BH/star binary meets field stars, one after another | the binary as it is when an encounter with a black hole comes due |
| `codeB` | `2 Code B.nb` | it meets that 1G black hole | a black-hole binary, if the newcomer takes the star's place |
| `codeC` | `2 Code C.nb` | the black-hole binary meets stars and 1G black holes | a gravitational-wave merger |

Each stage reads what the stage before it wrote. Units: AU, yr, M☉.

The star of `codeA` can also come from tidal capture instead of from a
stellar binary: the `tc` stage (`TC Code A.nb`) then stands where `code0`
stands. See [Two ways to capture a star](#two-ways-to-capture-a-star).

[docs/experiment.md](docs/experiment.md) (also as PDF) explains the Code A
experiment step by step and compares three integrators on it. Next to each
notebook lies a copy that reads without Mathematica (`2 Code 0.txt`, …), made
with `python nb2txt.py notebook.nb notebook.txt`. The PDF is built in `docs/`
with

    pandoc experiment.md -o experiment.pdf --pdf-engine=xelatex -V geometry:margin=1in -V fontsize=11pt

## Build

    make                                   # into bin/
    make SPACEHUB=/path/to/SpaceHub/src    # if SpaceHub is not in ~/Repositories/SpaceHub

Needs a C++17 compiler and the SpaceHub headers at commit `1a0ca3fd` or
later: before it, about one encounter in ten million sends SpaceHub's
Bulirsch–Stoer iterator into a step that never ends.

## Run

    bin/scatter_spacehub all --dir out/pipeline      # the four stages, one after the other

    bin/scatter_spacehub code0 --dir out/pipeline    # or one at a time:
    bin/scatter_spacehub codeA --dir out/pipeline    #   each reads the file of the
    bin/scatter_spacehub codeB --dir out/pipeline    #   stage before it in --dir
    bin/scatter_spacehub codeC --dir out/pipeline
    bin/scatter_spacehub lists --dir out/pipeline    #   and their lists in one file

`--in FILE` hands a stage another input than the previous stage's file, and
`--help` lists every option. The ones that change the result:

| option | meaning | default |
|---|---|---|
| `--tidal X` | a star left bound to a black hole is disrupted if its pericentre is within X tidal radii | 1 |
| `--nreal N` | `code0`: scatterings per stellar binary; `codeA`, `codeC`: histories per binary | 100 |
| `--nreal0 N`, `--nrealA N`, `--nrealC N` | the same for one stage only, so that one command can give every stage its own number | `--nreal` |
| `--nrealB N` | `codeB`: scatterings per BH/star binary, each by a black hole on another orbit | 1 |
| `--seed S` | random seed | 12345 |
| `--acc X` | relative tolerance of the integrator | 1e-12 |
| `--nsys N` | use only the first N mass pairs of the input | all |

Work is shared out among all cores (`--threads N` to use fewer); in a batch
job, among the cores the job was given. The result does not depend on the
number of threads: every scattering and every history draws its random
numbers from the seed and its own indices.

On 8 cores the two mass pairs of `2 Code 0.nb` take about two minutes, nearly
all of it in `codeA` and `codeC`.

SpaceHub prints `Reach max iteration loop number!` when a step does not
converge. If it was the first step, the integration is started again with a
ten times shorter first step; otherwise with a ten times looser tolerance.
Either way the run goes on.

### Where an encounter starts and ends

The incoming body starts ten tidal radii of the binary away,
`10 a (m_single / m_binary)^(1/3)`, and no closer than `10 a`. From infinity
to that point it moves on the hyperbola given by its speed at infinity and
its pericentre, so nothing is integrated before.

The outcome is looked at for the first time when a body leaving as it came
would be twice as far out as it started, and then each time the time since
pericentre has doubled. The integration ends at the first look that finds
the encounter settled: three single bodies, or a binary and an unbound body
that recedes and is as far from the binary as an incoming body starts. An
encounter not settled by then (a bound triple, a resonance) ends when a body
leaving as it came would be `100 a` away, and is classified as it is. The
times of the looks follow from the initial conditions alone, so all three
integrators stop at the same time.

The notebooks start every encounter `100 a` away and integrate for a fixed
time, in to pericentre and back out to `100 a`; so did this program until
17 November 2026. On the same inputs the two ways give the same statistics
(outcome counts: χ² = 4.5, 4.2, 3.4 and 6.7 for 5, 5, 7 and 11 degrees of
freedom in `code0`, `codeA`, `codeB` and `codeC`), and the new one needs a
fifth of the integration steps in `codeA`.

### The tidal check

After every encounter, in all four stages: if the tightest bound pair is a
black hole and a star, and the pericentre `a(1-e)` of their orbit lies inside
`X` tidal radii `R★ (M_BH/M☉)^(1/3)`, the star is disrupted. The radius of
a star is `R★ = 0.005 AU`, the `RS` of `2 Code 0 - Yihan.nb` (1.075 solar
radii). Runs made before 18 November 2026 took the solar radius, 0.00465 AU,
the run behind `docs/results.md` among them: with the larger one, 1.2% more
of the stars captured in Code 0 are disrupted at once (5.6% at the innermost
radius). Every table says in its first line which radius it was made with
(`rstar`). So it is in
`codeA` when its orbit would shrink by gravitational waves before anything
else happens; the notebooks call that a merger.

| stage | where a star can end up bound to a black hole | what disruption means |
|---|---|---|
| `code0` | the captured star | no BH/star binary forms: outcome `TIDAL_DISRUPTION` |
| `codeA` | the companion after each stellar encounter | the history ends (`fate` 6) |
| `codeA` | the companion, if gravitational waves bring it in before the next star or the black hole arrives | the same: a star does not merge with a black hole, it is pulled apart on the way in |
| `codeB` | the BH/star binary survives, or the star pairs with the 1G black hole | outcome `TIDAL_DISRUPTION` |
| `codeC` | a star takes the place of one of the two black holes | the history ends: outcome `TIDAL_DISRUPTION` |

`2 Code 0 - Yihan.nb` makes the same check on the captured star and calls
the outcome `TDE_BHS`; on its 1200 scatterings the six outcomes agree with
this program's (χ² = 3.3 for 5 degrees of freedom). The other notebooks have
no such check: this is where the C++ results differ from theirs by more than
chance. A star left on a wide orbit around a new
black-hole binary (`(m1/m3)/m2` in `codeB`) is not checked.

## Two ways to capture a star

Code A follows a black hole and the star bound to it. That binary can be
made in two ways, and the rest of the pipeline is the same for both:

| | `code0` | `tc` |
|---|---|---|
| the star comes from | a stellar binary that the black hole breaks up | a single star that loses its orbital energy to the tide of the black hole |
| how the orbit is found | a three-body integration for every scattering | two rules, no integration |
| binaries per mass pair and radius | 30% of `--nreal0` scatterings | one |
| pericentre of the star | 5.6 tidal radii in the median | 0.97 to 1.28 tidal radii at 0.11 to 2 pc |

    bin/scatter_spacehub all    --in data/params0_1000.csv --dir out/binary ...   # code0, codeA, codeB, codeC, lists
    bin/scatter_spacehub all-tc --in data/params0_1000.csv --dir out/tidal  ...   # tc,    codeA, codeB, codeC, lists

or stage by stage, `tc` first and then `codeA --in out/tidal/tc_spacehub.csv`.
`tc` reads the mass pairs and radii of the same input table as `code0` (its
other columns are not read), so the two runs are on the same grid. With
`slurm/submit.sh`, set `CAPTURE=tidal`.

The notebooks hold the tidal-capture orbits as numbers: the lists `A1sS` and
`E1sS`, 1000 systems at ten radii from 0.01 to 2 pc (also in
`data/bhstar_tidal_capture.csv`). They follow two rules, which `tc` uses
(`include/bhs/tc.hpp`):

    a   = G M1 / sigma^2                                  orbital speed of the star = velocity dispersion
    r_p = 0.010156 AU (M1 / Msun)^(1/3) (sigma yr/AU)^(-2/10.65)
    e   = 1 - r_p / a

These reproduce all 10,000 orbits of the notebooks to 13 digits. The
constant and the exponent are read off the numbers: the formula itself is
not in the notebooks.

Checked on 18 November 2026 against `TC Code A.nb`, which runs the first two
systems at its ten radii with 1000 histories each: the fraction of histories
in which the star is disrupted agrees in the 20 cells (χ² = 21.8), and so
does the mean time to the black-hole encounter.

What follows from the orbits:

- A star captured by tides sits at about one tidal radius, where the tidal
  check disrupts it. On the notebook's grid half of the binaries start inside
  the tidal radius (all of them at 0.01 to 0.11 pc), 88% of the histories end
  with the star disrupted, 62% of those at the first encounter, and 9% reach
  the black-hole encounter. As in the notebook, the star is tested after an
  encounter and not before the first.
- There is no capture experiment, hence no `Num0` and no `f1` in the lists
  file: how often a black hole captures a star by tides has to come from the
  rate of that process, which is not in the notebooks.
- The notebook's ends of a history are this program's `fate`: `REACH_BH` 0,
  `TDE_BHS` and `BHS_MERGER` 6, `BHS_LOST` 3, 4 and 5.
- `analysis.py` expects a Code 0 table and does not yet read a run of
  `all-tc`.

## Output

Every stage writes two files, named after the stage and the integrator:

- `code0_spacehub.csv`: a table with one row per scattering or history. This is
  what the next stage reads, and what Python reads most easily. Everything a
  stage found is here.
- `code0_spacehub.m`: the lists the notebooks plot, as Mathematica input.

`lists` (and `all`, when its stages are done) gathers the four `.m` files
into `lists_spacehub.m`: the lists in the order of the notebook cell they
are pasted into, one per line as `name=value;`. In Mathematica,
`Get["lists_spacehub.m"]` defines them all.

| list | from | what it holds |
|---|---|---|
| `M12BHB` | `codeB` | the mass pairs `{M1, MBH}` |
| `RBHB` | `codeB` | the radii [AU], once per mass pair |
| `TBHB` | `codeB` | for every black-hole binary formed: the time its BH/star binary spent in Code A [yr] |
| `F3SS` | `codeC` | the mean merger fraction of the black-hole binaries of a cell |
| `F3BHBSS` | `codeC` | per black-hole binary: the fraction of its histories that merge |
| `TMBHBSS` | `codeC` | per black-hole binary: the merger times [yr] |
| `M3GBHBSS` | `codeC` | per black-hole binary: the masses of the remnants |
| `M2MERGEBHBSS` | `codeC` | per black-hole binary: the masses of the companions merged with |
| `Num0`, `NUMA`, `NUMC` | `code0`, `codeA`, `codeC` | the `--nreal` of the three stages |
| `NUMB` | `codeB` | its `--nrealB`; written only if that is not 1, which the notebooks have no name for |
| `f1` | `code0` | the fraction of scatterings that leave a BH/star binary (the notebook's `f1IN0`) |
| `f2` | `codeB` | the fraction of Code A histories that end as a black-hole binary: `PBHS PBHB` |
| `X1` | `codeA` | the mean time to the black-hole encounter [yr]; 0 where no history reached it |
| `dataTM` | `codeC` | per mass pair: `{R [pc], mean merger time [yr]}` at the radii with a merger |

Lists are nested as in the notebooks: `[[kk]]` mass pair, `[[jj]]` radius,
then, where there is more, `[[bb]]` black-hole binary and its entries. The
other lists of the notebooks are no longer written: the tables hold what
they held.

Outcomes carry the notebooks' labels (`m1/m3-m2`, `UNRESOLVED / RESONANT`,
`MERGER`, `TARGET_2G_NO_BHB_AFTER_BH:m1-m2-m3`, …). Two labels have another
cause here than in the notebooks: `NDSOLVE_FAILED` and `INTEGRATION_FAILED`
mean that SpaceHub gave up on an encounter. `EARLY_STOP` does not occur, and
`TIDAL_DISRUPTION` is new.

### Table columns

Every table starts with `kk`, `jj` and the indices of the row within its cell,
then `R_pc`, `M1` (the 2G black hole) and `MBH` (the 1G black hole of Code B).
Its first line says how it was made and holds the grid of masses and radii.

| table | further columns |
|---|---|
| `code0` | `ii` scattering · `outcome` · `a`, `e` of the binary left behind · `v_inf`, `rp` of the incoming stellar binary · `t_end` integrated time · `dE` energy error |
| `codeA` | `h` which Code 0 capture, `ii` history · `a0`, `e0` initial orbit · `a`, `e`, `tb` orbit and elapsed time at the black-hole encounter, zero unless `fate` is 0 · `fate` · `n_enc` stellar encounters · `a_last`, `e_last` orbit when the history ended |
| `codeB` | `n` which BH/star binary, `rep` which scattering of it, `h`, `ii` the Code A history it came from · `a0`, `e0`, `tb` the incoming binary · `outcome` · `a`, `e` tightest bound pair · `a_out`, `e_out`, `i_mut` third body of a bound triple · `v_recoil` of the pair if the third body escaped |
| `codeC` | `bb` which black-hole binary, `ii` history · `a0`, `e0` · `outcome` · `tm` merger time, `m2_merge` companion at merger, `m_3g` remnant mass, zero unless `MERGER` · `n_star`, `n_bh` encounters · `n_exch` black holes that took the companion's place |

`fate` in `codeA`: 0 black-hole encounter reached · 3 black hole lost its
companion · 4 ionised · 5 unresolved · 6 star tidally disrupted. Tables
written before 17 November 2026 also hold 1 and 2, a merger of black hole
and star by gravitational waves before the black hole or the next star
arrives: these are now 6.

## Input

`data/params0.csv` starts the pipeline: one row for every mass pair at every
radius, with the columns of the notebook's `Params0`.

| column | meaning |
|---|---|
| `kk`, `jj` | mass pair and radius, counted from 1 |
| `M1`, `MBH` | 2G black hole; 1G black hole it meets in Code B [M☉] |
| `R_AU` | distance from the central black hole |
| `a_ss`, `e_ss` | orbit of the stellar binary |
| `r_tide` | radius inside which the black hole splits the binary (`Rps`); the pericentre is drawn within twice this |
| `t_a` | analytic encounter time `TA` [yr]; handed on to `TA0`, not used |

From the notebook:

    Export["params0.csv", Prepend[
      Flatten[Table[Join[{kk, jj}, Params0[[kk, jj, 1 ;; 7]]], {kk, nMass}, {jj, nR}], 1],
      {"kk", "jj", "M1", "MBH", "R_AU", "a_ss", "e_ss", "r_tide", "t_a"}]]

Or, without Mathematica, from the text copy of a notebook that holds the
lists `M12S`, `Routs`, `ABIs`, `EBIs`, `Rps` and `TA`:

    python nb2txt.py "2 Code 0 - Yihan.nb" "2 Code 0 - Yihan.txt"
    python make_params0.py --notebook "2 Code 0 - Yihan.txt"      # -> data/params0_1000.csv

`data/params0_1000.csv` is that table: the 1000 mass pairs of
`2 Code 0 - Yihan.nb` at its six radii, with its numbers as they stand.
`data/params0.csv` holds the first two of them, the mass pairs of
`2 Code 0.nb`. The notebook's columns follow closed formulas (the stellar
binary scales with the radius, its breakup radius and encounter time with the
black-hole mass), so the table can also be written for any list of masses:

    python make_params0.py masses.txt -o params0_mine.csv

For the 1000 mass pairs of `data/M12S.txt`, which are those of the notebook,
the formulas give the notebook's table to 15 digits.

`data/M1G.txt` holds the masses of the 1G black holes `codeC` draws from
(`M1G` in the notebook), one per line.

`data/bhstar_tidal_capture.csv` is the input of the original `Code A.nb`: 1000
BH/star binaries from tidal capture at ten radii. To run Code A on it instead
of on the output of Code 0:

    bin/scatter_spacehub codeA --in data/bhstar_tidal_capture.csv --nsys 100
    python export_aet.py codeA_spacehub.csv     # -> codeA_spacehub.txt

The second line writes the result as the three lists of the original notebook,
`AINSS`, `EINSS` and `TBsS`, indexed `[[system]][[radius]][[realisation]]` with
0 where the black-hole encounter was not reached.

Any table with columns `kk,jj,a,e` and the grid in its first line works the
same way.

## The full run and its figures

The 1000 mass pairs of `data/params0_1000.csv` at six radii:

    B=bin/scatter_spacehub
    D=out/run_200_5_8_20; mkdir -p $D
    caffeinate -i $B code0 --in data/params0_1000.csv --nreal 200 --dir $D
    caffeinate -i $B codeA --nreal 5  --dir $D
    caffeinate -i $B codeB --nrealB 8 --dir $D
    caffeinate -i $B codeC --nreal 20 --dir $D
    $B lists --dir $D                       # -> lists_spacehub.m

    python analysis.py $D                   # seven figures into docs/figures/, numbers to the terminal

This is the run behind [docs/results.md](docs/results.md) (also as PDF), which
shows the figures and says what they mean. `analysis.py` works on any run of
the four stages; it needs numpy, pandas and matplotlib. Its errors take the
Code 0 captures as the independent draws: what descends from one capture
shares its history.

The run took four and a half hours on 8 cores (Apple M2) and 1.8 GB of
tables, and gave 45,496 black-hole binaries and 284,258 mergers: the
probability that a Code 0 scattering ends in a merger to 3.1% at the
innermost radius and 1.6% at the outermost, or 4.4% and 2.4% with the choice
of the mass pairs counted in. A stage that is interrupted starts again from its
beginning, which is why the stages are run one at a time; `caffeinate -i`
keeps a Mac from going to sleep meanwhile. An earlier run with 20, 8 and 40
at the three stages and one Code B scattering, made when encounters still
started `100 a` away, is in `out/full/`; the two agree within their errors.

### Choosing the numbers

What limits the result is the number of black-hole binaries, of which Code B
forms one in 120 scatterings, and what a binary is worth depends on how the
repeats are spread. Measured on 17 November 2026:

- Black-hole binaries differ from each other more than histories of one
  binary do (merger probability 0.32 with a scatter of 0.24), so a few
  histories of many binaries beat many histories of few.
- Captured stars differ too, and Code 0 costs next to nothing: more
  scatterings with fewer Code A histories each.
- Code B is the cheapest stage to repeat. With `--nrealB` every BH/star
  binary that reaches it is scattered several times.
- The yield of black-hole binaries differs between mass pairs by a factor of
  its own size (relative scatter 0.5 to 1.1). An average over the pairs
  cannot be better than that scatter over the root of their number: 5 to 11%
  for 100 pairs, 1.5 to 3.4% for 1000.

These led to the numbers of the full run above: 200 Code 0 scatterings, 5
Code A histories, 8 Code B scatterings and 20 Code C histories. The
notebooks' 100 at every stage with one Code B scattering (`--nreal 100`, the
default) gives the same accuracy in 34 hours instead of four and a half.

## On a cluster (SLURM)

### Build and submit

    slurm/build.sh                 # with the g++ that is there
    slurm/build.sh gcc/12.2.0      # or with a module: `module avail gcc` lists what there is

    slurm/submit_runs.sh           # the two runs: out/binary and out/tidal
    slurm/submit_runs.sh tidal     # or one of them

`slurm/build.sh` clones SpaceHub into `../SpaceHub` if it is not there (or
takes `$SPACEHUB`, its `src` directory), builds the three programs and runs a
test of twelve scatterings. A module it loaded is written to
`slurm/modules.sh`, and every job loads it again.

`slurm/submit_runs.sh` submits the run in which the black hole takes its star
from a stellar binary and the run in which it captures the star by tides
(see [Two ways to capture a star](#two-ways-to-capture-a-star)), both with the
notebooks' 100 at every stage and one Code B scattering. Its first lines hold
what can be set, and the environment sets it too:

    SBATCH_PARTITION=WorkA CPUS=48 TIME=06:00:00 slurm/submit_runs.sh
    NUMA=1000 slurm/submit_runs.sh tidal

| setting | meaning | default |
|---|---|---|
| `NUM0`, `NUMA`, `NUMB`, `NUMC` | `--nreal0`, `--nrealA`, `--nrealB`, `--nrealC` | 100, 100, 1, 100 |
| `PARAMS`, `OUT` | input table of mass pairs and radii; where `binary/` and `tidal/` go | `data/params0_1000.csv`, `out` |
| `SBATCH_PARTITION` | the partition the jobs go to | `WorkC` |
| `SBATCH_ACCOUNT` | the account they are charged to | the one SLURM ties you to on that partition |
| `PARTS_BINARY`, `PARTS_TIDAL` | tasks of the job arrays, a node each | 16, 2 |
| `CPUS` | cores of every task | what every node of the partition has |
| `TIME`, `MEM` | time limit of every job; memory of the jobs that hold a whole table | 12:00:00, 16G |

The first run is about 240 core-hours on Apple M2 cores: 18 million Code A
histories and 5.6 million Code C histories. The second is a few core-hours:
600,000 Code A histories, of which 13% reach the black hole, and some 190
black-hole binaries. If a job fails, run the script again: what is finished
is left alone.

By hand, the build is

    git clone https://github.com/YihanWangAstro/SpaceHub.git ../SpaceHub     # commit 1a0ca3fd or later
    make CXX=g++ SPACEHUB=$PWD/../SpaceHub/src

A C++17 compiler is all it needs. Build once, on the kind of node the jobs
run on, and use that build for a whole run. With a compiler from a module the
programs need that module when they run as well, which is what
`slurm/modules.sh` is for; `LDFLAGS="-static-libstdc++ -static-libgcc"` links
the compiler's libraries in instead, where the static library is installed
(`ld: cannot find -lstdc++` says it is not).

### One node

The program takes its number of threads from the job
(`SLURM_CPUS_PER_TASK`, or the cores it is bound to), and in a log file its
progress is a line a minute. `slurm/run.sbatch` runs whatever command it is
given:

    sbatch --cpus-per-task=64 --time=04:00:00 slurm/run.sbatch \
        bin/scatter_spacehub all --in data/params0_1000.csv --dir out/run \
        --nreal0 200 --nrealA 5 --nrealB 8 --nrealC 20

That run is about 35 core-hours on an Apple M2. The cores of a compute node
are usually slower one by one: expect one to two hours on 64 of them.

### Many nodes

A stage can be run in parts, each by a job of its own:

    bin/scatter_spacehub codeA --nreal 5 --dir out/run --part 3/16   # -> codeA_spacehub.part03of16.csv
    bin/scatter_spacehub codeA --dir out/run --merge 16              # the 16 parts -> the table and the lists

`--part I/N` runs the I-th of N equal shares of the jobs of the stage.
`--merge N` puts the N files together; it takes `--nreal` and the like from
the parts, refuses parts that were not written to their end or not made
alike, and leaves a table that is there untouched if it refuses. The table
is, byte for byte, the one a single run writes. A part that is there and
complete is not run again.

`slurm/submit.sh` does this with job arrays for the two long stages, every
job waiting for the one before it:

    export SBATCH_PARTITION=...  SBATCH_ACCOUNT=...       # what your cluster asks for
    PARTS=16 CPUS=32 TIME=12:00:00 slurm/submit.sh out/run data/params0_1000.csv \
        --nreal0 200 --nrealA 5 --nrealB 8 --nrealC 20

With `CAPTURE=tidal` the `tc` stage takes the place of Code 0;
`slurm/submit_runs.sh` calls this script once for each of the two.

It submits Code 0, an array of `PARTS` tasks for Code A and a job that puts
them together, Code B, the same pair for Code C, and `lists`. Logs go to
`out/run/log/`. If a task fails (out of time, a node down), the jobs behind
it never start: submit again with the same arguments. Stages whose lists are
there are left alone, and so are the parts that are complete. The part files
can be deleted once their stage is put together.

### What carries over from one machine to another

A seed means the same initial conditions under every compiler: the random
distributions are written out in `encounter.hpp`, because the standard
library's differ (GCC's gave other speeds than clang's for the same seed).
What follows from them does not carry over exactly. Compilers round
differently in the last digit, and the three-body problem amplifies that:
between GCC 15 and Apple clang 21 on the same Mac, all of 120 Code 0
scatterings had the same outcome and orbits equal to 11 digits or more,
while the histories of Code A, each a chain of dozens of encounters, went
different ways. Tables
from two builds agree statistically, not row by row; the parts of one run
should come from one build.

Checked on 17 November 2026: the build with GCC 15 and Apple clang 21 (on
macOS), parts and merging against a single run, and `slurm/submit.sh`
against a stand-in for `sbatch`, with a failing task and a second
submission. Not yet run on Linux or under a real SLURM.

## Agreement with the notebooks

The random numbers differ, so results agree statistically. The comparison
below was made on 16 November 2026 with the tidal check switched off, as the
notebooks have none, and with encounters starting `100 a` away, as in the
notebooks; the program no longer offers either. `code0`, `codeB` and
`codeC` were run on the lists the notebooks themselves contain:

| stage | input | notebook | C++ |
|---|---|---|---|
| `code0` | `Params0`, 1200 scatterings | 401 BH/star binaries, 595 stellar binaries survive, 193 all unbound, 11 unresolved | 387, 609, 190, 14 |
| `codeA` | the BH/star binaries of its own Code 0, 100 histories each | 36.5% reach the black-hole encounter; `PBHS` from 0.16 at the innermost radius to 0.56 at the outermost | 36.7%; 0.15 to 0.59 |
| `codeB` | 14,654 BH/star binaries of the notebook | 91 black-hole binaries | 101 |
| `codeC` | 91 black-hole binaries of the notebook, 9100 histories | 2676 mergers, 2496 lost the companion to a black hole, 1386 at the limit of 1000 scatterings | 2540, 2500, 1522 |

For `codeC` the merger fractions of the 91 binaries, compared one by one,
differ by no more than their sampling noise (χ² = 88). The deterministic
parts (cluster model, merger time, and the diagnostics of Code B, which the
program still wrote then) reproduce the numbers printed in the notebooks to
all digits.

The notebook of Code A in its pipeline form was not at hand: what `codeA`
reads and writes was worked out from the lists `2 Code B.nb` takes from it.

`out/pipeline/` holds a run of the whole pipeline.

## Other integrators, and looking at single encounters

The stages do not know which integrator runs. Two more are built, for
comparison:

| program | integrator | `--acc` means |
|---|---|---|
| `bin/scatter_spacehub` | SpaceHub regularised chain (`methods::DefaultMethod`) | relative tolerance, default 1e-12 |
| `bin/scatter_bs` | Bulirsch–Stoer extrapolation, Cartesian coordinates | relative tolerance, default 1e-12 |
| `bin/scatter_leapfrog` | time-transformed leapfrog, 2nd order | steps per orbit of the initial binary, default 1000 |

All three offer the same commands, and three tools besides the stages:

    bin/scatter_bs encounters --n 20000 --nsys 100   # 20000 independent Code A encounters, one row each
    bin/scatter_bs traj --encounter 1950 --nsys 100  # trajectory of one of them
    bin/scatter_bs profile                           # cluster density / velocity-dispersion table

With the same `--seed` every integrator is given identical initial conditions,
so `encounters` and `traj` can be compared row by row, and so can the
scatterings of `code0` and `codeB`. Histories (`codeA`, `codeC`) start
identically but drift apart, because the three-body problem is chaotic;
compare those statistically.

    python plot_results.py profile
    python plot_results.py results codeA_bs.csv
    python plot_results.py traj traj_bs.csv

    python compare.py traj        traj_bs.csv traj_leapfrog.csv traj_spacehub.csv
    python compare.py encounters  encounters_bs.csv encounters_leapfrog.csv encounters_spacehub.csv
    python compare.py mc          codeA_bs.csv codeA_leapfrog.csv codeA_spacehub.csv

Needs numpy, matplotlib and scipy.

## Layout

    include/bhs/constants.hpp     units and model parameters
    include/bhs/vec3.hpp          3-vector
    include/bhs/three_body.hpp    state, equations of motion, integrator interface
    include/bhs/cluster.hpp       cluster model, merger and encounter timescales, tidal radius
    include/bhs/encounter.hpp     initial conditions of an encounter, outcome classification
    include/bhs/pipeline.hpp      what the stages share: grid, stage files, settings, random numbers
    include/bhs/code0.hpp         stage 0: black hole + stellar binary
    include/bhs/tc.hpp            in its place: the binary tidal capture leaves
    include/bhs/codeA.hpp         stage A: BH/star binary among field stars
    include/bhs/codeB.hpp         stage B: BH/star binary + black hole
    include/bhs/codeC.hpp         stage C: black-hole binary until merger
    include/bhs/driver.hpp        command line; runs any stage on all cores, or in parts
    include/bhs/io.hpp            CSV files, thread pool
    include/bhs/mathematica.hpp   writing Mathematica lists
    src/scatter_*.cpp             one integrator each
    slurm/build.sh                builds on a cluster: module, SpaceHub, make, a test
    slurm/submit_runs.sh          submits the two runs, binary and tidal
    slurm/submit.sh               submits the stages of one run to SLURM, the long ones as job arrays
    slurm/run.sbatch              batch script: runs a command, as a task of a job array one part of a stage
    data/                         input tables
    make_params0.py               the input table of Code 0 for a list of mass pairs
    analysis.py                   figures and numbers of a run of the whole pipeline
    plot_results.py, compare.py   figures of one Code A run; comparison of the integrators
    export_aet.py                 the lists of the original Code A.nb from a codeA table
    nb2txt.py                     a notebook as readable text
