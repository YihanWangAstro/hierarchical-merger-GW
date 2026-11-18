#!/usr/bin/env python3
"""Compare the integrator versions of the BH/star binary scattering code.

    python compare.py traj        traj_bs.csv traj_leapfrog.csv traj_spacehub.csv
    python compare.py encounters  out/encounters_*.csv
    python compare.py mc          codeA_bs.csv codeA_leapfrog.csv codeA_spacehub.csv

Give one output file per integrator, produced with the same command and seed.
Differences are measured against the reference integrator (--ref, default
spacehub).  Each command prints a summary table and writes a PNG.
"""
import argparse

import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import ks_2samp

from plot_results import (BLUE, FATES, G, GREEN, GREY, ORANGE, PURPLE, SKY, VERMILLION, finish, load_csv,
                          panel_label, set_style)

# One colour and line style per integrator, the same in every figure.
STYLES = {
    "bs": dict(color=BLUE, ls="-", name="Bulirsch-Stoer", unit="tol"),
    "leapfrog": dict(color=VERMILLION, ls="--", name="leapfrog", unit="steps/orbit"),
    "spacehub": dict(color=GREEN, ls="-.", name="SpaceHub", unit="tol"),
}
SPARE_COLOURS = [PURPLE, ORANGE, SKY, GREY]


class Run:
    """One output file together with how to draw it."""

    def __init__(self, path, spare):
        self.data, self.meta = load_csv(path)
        self.key = self.meta.get("integrator", path)
        style = STYLES.get(self.key, dict(color=GREY, ls="-", name=self.key, unit="acc"))
        accuracy = self.meta.get("accuracy", "?")
        self.name = style["name"]
        self.label = (f"{self.name} ({accuracy} {style['unit']})" if style["unit"] != "tol"
                      else f"{self.name} (tol {accuracy})")
        self.line = dict(color=style["color"], ls=style["ls"])
        if spare is not None:   # second file of the same integrator
            self.line = dict(color=spare, ls=":")


def load_runs(paths, ref_key):
    runs, seen = [], set()
    spare = iter(SPARE_COLOURS)
    for path in paths:
        key = load_csv(path)[1].get("integrator", path)
        runs.append(Run(path, next(spare) if key in seen else None))
        seen.add(key)
    ref = next((r for r in runs if r.key == ref_key), runs[0])
    return runs, ref


def cdf(ax, values, floor=None, **line):
    """Cumulative distribution as a step curve; values below `floor` are drawn at the floor."""
    x = np.sort(np.asarray(values, dtype=float))
    if floor is not None:
        x = np.maximum(x, floor)
    ax.step(x, np.arange(1, len(x) + 1) / len(x), where="post", **line)


# ----------------------------------------------------------------------------
# One encounter
# ----------------------------------------------------------------------------
def body(d, i):
    return (np.column_stack([d[f"{c}{i}"] for c in "xyz"]),
            np.column_stack([d[f"v{c}{i}"] for c in "xyz"]))


def energy(d, masses):
    r, v = zip(*(body(d, i) for i in (1, 2, 3)))
    e = sum(0.5 * m * np.sum(vi**2, axis=1) for m, vi in zip(masses, v))
    for i, j in [(0, 1), (0, 2), (1, 2)]:
        e -= G * masses[i] * masses[j] / np.linalg.norm(r[i] - r[j], axis=1)
    return e


def cmd_traj(args):
    runs, ref = load_runs(args.files, args.ref)
    n = min(len(r.data) for r in runs)
    a0 = float(ref.meta["a0"])
    masses = (float(ref.meta["M1"]), 1.0, 1.0)
    t = ref.data["t"][:n]

    print(f"encounter {ref.meta.get('encounter')}: M1 = {masses[0]:.4g}, a0 = {a0:.5g} AU, e0 = {ref.meta.get('e0')}")
    for r in runs:
        print(f"  {r.label:36s} outcome {r.meta['outcome']:12s} a_bin = {float(r.meta['a_bin']):.10g}"
              f"  e_bin = {float(r.meta['e_bin']):.10g}")

    fig, axes = plt.subplots(3, 1, figsize=(7.5, 8.5), sharex=True)

    # A: the encounter itself (reference integrator)
    ax = axes[0]
    pos = [body(ref.data, i)[0][:n] for i in (1, 2, 3)]
    for (i, j), col, lab in zip([(0, 1), (0, 2), (1, 2)], [VERMILLION, "black", BLUE], ["m1-m2", "m1-m3", "m2-m3"]):
        ax.semilogy(t, np.linalg.norm(pos[i] - pos[j], axis=1) / a0, color=col, lw=0.7, label=lab)
    ax.set_ylabel(r"separation / $a_0$")
    ax.legend(ncol=3, loc="upper center")
    ax.set_title(f"encounter {ref.meta.get('encounter')}: outcome {ref.meta['outcome']} ({ref.label})", fontsize=10)

    # B: difference of the binary's separation vector from the reference
    ax = axes[1]
    sep_ref = pos[1] - pos[0]
    for r in runs:
        if r is ref:
            continue
        sep = body(r.data, 2)[0][:n] - body(r.data, 1)[0][:n]
        ax.semilogy(t, np.maximum(np.linalg.norm(sep - sep_ref, axis=1) / a0, 1e-16), lw=0.8, label=r.label, **r.line)
    ax.set_ylabel(r"$|\Delta\mathbf{r}_{12}|\,/\,a_0$" + f"  vs {ref.name}")
    ax.legend(loc="lower right")

    # C: energy error of each integrator
    ax = axes[2]
    for r in runs:
        e = energy(r.data, masses)[:n]
        ax.semilogy(t, np.maximum(np.abs(e / e[0] - 1), 1e-16), lw=0.8, label=r.label, **r.line)
    ax.set_ylabel(r"$|\Delta E / E_0|$")
    ax.set_xlabel(r"$t$ (yr)")
    ax.legend(loc="lower right")

    for ax, letter in zip(axes, "ABC"):
        panel_label(ax, letter)
    finish(fig, args.out or "compare_traj.png", args.show)


# ----------------------------------------------------------------------------
# Many independent encounters with identical initial conditions
# ----------------------------------------------------------------------------
def cmd_encounters(args):
    runs, ref = load_runs(args.files, args.ref)
    n = min(len(r.data) for r in runs)
    for r in runs:
        r.data = r.data[:n]
        if not np.array_equal(r.data["v_inf"], ref.data["v_inf"]):
            raise SystemExit("the files do not contain the same encounters (different seed or options?)")

    outcomes, first = np.unique(ref.data["outcome"], return_index=True)
    names = {o: ref.data["outcome_name"][i] for o, i in zip(outcomes, first)}
    all_codes = sorted(set(np.concatenate([np.unique(r.data["outcome"]) for r in runs])))
    for r in runs:
        for o, i in zip(*np.unique(r.data["outcome"], return_index=True)):
            names.setdefault(o, r.data["outcome_name"][i])

    print(f"{n} encounters; differences relative to {ref.label}")
    print(f"{'integrator':36s} {'same outcome':>12s} {'median|da/a|':>13s} {'99%|da/a|':>10s} "
          f"{'median|dE/E|':>13s} {'max|dE/E|':>10s} {'failed':>7s} {'retries':>8s} {'ms/enc':>7s}")
    diffs = {}
    for r in runs:
        same = r.data["outcome"] == ref.data["outcome"]
        bound = same & (ref.data["a_bin"] > 0)
        diffs[r] = np.abs(r.data["a_bin"][bound] / ref.data["a_bin"][bound] - 1)
        da = diffs[r] if r is not ref else np.zeros(1)
        print(f"{r.label:36s} {same.mean():12.4%} {np.median(da):13.1e} {np.percentile(da, 99):10.1e} "
              f"{np.median(r.data['dE']):13.1e} {r.data['dE'].max():10.1e} {np.sum(r.data['ok'] == 0):7d} "
              f"{np.sum(r.data['retries']):8d} {r.data['wall_ms'].mean():7.2f}")

    fig, axes = plt.subplots(2, 2, figsize=(10, 7.5))

    ax = axes[0, 0]   # A: energy conservation
    for r in runs:
        cdf(ax, r.data["dE"][r.data["ok"] == 1], floor=1e-16, label=r.label, **r.line)
    ax.set_xscale("log")
    ax.set_xlabel(r"energy error $|\Delta E/E|$ per encounter")
    ax.set_ylabel("cumulative fraction of encounters")
    ax.legend(loc="lower right")

    ax = axes[0, 1]   # B: final semi-major axis against the reference
    for r in runs:
        if r is not ref:
            cdf(ax, diffs[r], floor=1e-16, label=r.label, **r.line)
    ax.set_xscale("log")
    ax.set_xlabel(r"$|a_{\rm bin}/a_{\rm bin,ref}-1|$ (same outcome), ref = " + ref.name)
    ax.set_ylabel("cumulative fraction of encounters")
    ax.legend(loc="upper left")

    ax = axes[1, 0]   # C: outcome statistics
    width = 0.8 / len(runs)
    for k, r in enumerate(runs):
        counts = [np.sum(r.data["outcome"] == o) for o in all_codes]
        ax.bar(np.arange(len(all_codes)) + (k - (len(runs) - 1) / 2) * width, counts, width, color=r.line["color"],
               edgecolor="white", linewidth=0.4, label=r.label)
    ax.set_yscale("log")
    ax.set_xticks(np.arange(len(all_codes)))
    ax.set_xticklabels([names[o] for o in all_codes], rotation=30, ha="right")
    ax.tick_params(axis="x", which="both", top=False)
    ax.set_ylabel("number of encounters")

    ax = axes[1, 1]   # D: cost
    for r in runs:
        cdf(ax, r.data["wall_ms"], label=f"{r.name}: mean {r.data['wall_ms'].mean():.2f} ms", **r.line)
    ax.set_xscale("log")
    ax.set_xlabel("wall time per encounter (ms)")
    ax.set_ylabel("cumulative fraction of encounters")
    ax.legend(loc="lower right")

    for ax, letter in zip(axes.flat, "ABCD"):
        panel_label(ax, letter)
    finish(fig, args.out or "compare_encounters.png", args.show)


# ----------------------------------------------------------------------------
# Monte Carlo results
# ----------------------------------------------------------------------------
def cmd_mc(args):
    runs, ref = load_runs(args.files, args.ref)
    radii = np.unique(ref.data["R_pc"])

    print(f"{'fate':28s}" + "".join(f"{r.label:>34s}" for r in runs))
    for i, (name, _) in enumerate(FATES):
        print(f"{name:28s}" + "".join(f"{np.mean(r.data['fate'] == i):34.3%}" for r in runs))
    print(f"{'realisations':28s}" + "".join(f"{len(r.data):34d}" for r in runs))
    print(f"{'stellar encounters, mean':28s}" + "".join(f"{r.data['n_enc'].mean():34.2f}" for r in runs))

    reached = {r: r.data["fate"] == 0 for r in runs}
    quantities = [("a / a0", lambda r: r.data["a"][reached[r]] / r.data["a0"][reached[r]]),
                  ("1 - e", lambda r: 1 - r.data["e"][reached[r]]),
                  ("time to BH", lambda r: r.data["tb"][reached[r]])]
    print(f"\nKolmogorov-Smirnov test against {ref.label} (binaries that reach the BH encounter):")
    for r in runs:
        if r is ref:
            continue
        line = f"  {r.label:36s}"
        for name, f in quantities:
            statistic, pvalue = ks_2samp(f(r), f(ref))[:2]
            line += f"  {name}: D = {statistic:.4f}, p = {pvalue:.2f}"
        print(line)

    fig, axes = plt.subplots(2, 3, figsize=(13.5, 7.8))

    ax = axes[0, 0]   # A: fates
    width = 0.8 / len(runs)
    for k, r in enumerate(runs):
        frac = np.array([np.mean(r.data["fate"] == i) for i in range(len(FATES))])
        err = np.sqrt(frac * (1 - frac) / len(r.data))
        ax.bar(np.arange(len(FATES)) + (k - (len(runs) - 1) / 2) * width, frac, width, yerr=err, color=r.line["color"],
               edgecolor="white", linewidth=0.4, label=r.label, error_kw=dict(lw=0.8, capsize=1.5))
    ax.set_yscale("log")
    ax.set_xticks(np.arange(len(FATES)))
    ax.set_xticklabels([name for name, _ in FATES], rotation=30, ha="right")
    ax.tick_params(axis="x", which="both", top=False)
    ax.set_ylabel("fraction of realisations")
    fig.legend(*ax.get_legend_handles_labels(), loc="outside upper center", ncol=len(runs))

    ax = axes[0, 1]   # B: fraction reaching the BH encounter, by radius
    for k, r in enumerate(runs):
        frac, err = [], []
        for radius in radii:
            sel = r.data["R_pc"] == radius
            p = np.mean(reached[r][sel])
            frac.append(p)
            err.append(np.sqrt(p * (1 - p) / sel.sum()))
        ax.errorbar(radii * (1 + 0.03 * (k - 1)), frac, yerr=err, marker="o", ms=3, lw=1.1, capsize=2, **r.line)
    ax.set_xscale("log")
    ax.set_xlabel(r"$R$ (pc)")
    ax.set_ylabel("fraction reaching the BH encounter")

    ax = axes[0, 2]   # C: number of stellar encounters
    bins = np.arange(0.5, max(r.data["n_enc"].max() for r in runs) + 1.5)
    for r in runs:
        ax.hist(r.data["n_enc"], bins=bins, histtype="step", lw=1.2, **r.line)
    ax.set_yscale("log")
    ax.set_xlabel("stellar encounters per realisation")
    ax.set_ylabel("number of realisations")

    labels = [r"$a/a_0$ at the BH encounter", r"$1-e$ at the BH encounter", "time to the BH encounter (yr)"]
    for ax, (_, f), label in zip(axes[1], quantities, labels):   # D-F: distributions
        for r in runs:
            cdf(ax, f(r), **r.line)
        ax.set_xscale("log")
        ax.set_xlabel(label)
        ax.set_ylabel("cumulative fraction")

    for ax, letter in zip(axes.flat, "ABCDEF"):
        panel_label(ax, letter)
    finish(fig, args.out or "compare_mc.png", args.show)


# ----------------------------------------------------------------------------
def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = parser.add_subparsers(dest="cmd", required=True)
    common = argparse.ArgumentParser(add_help=False)
    common.add_argument("files", nargs="+", help="one output file per integrator")
    common.add_argument("--ref", default="spacehub", help="reference integrator [spacehub]")
    common.add_argument("-o", "--out", help="output image file")
    common.add_argument("--show", action="store_true", help="also open a window")
    for name, func in [("traj", cmd_traj), ("encounters", cmd_encounters), ("mc", cmd_mc)]:
        sub.add_parser(name, parents=[common]).set_defaults(func=func)
    args = parser.parse_args()
    set_style()
    args.func(args)


if __name__ == "__main__":
    main()
