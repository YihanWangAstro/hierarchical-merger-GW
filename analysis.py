#!/usr/bin/env python3
"""Figures and summary numbers for a run of the whole pipeline.

    python analysis.py out/full                  # figures into docs/figures/, numbers to the terminal
    python analysis.py out/full -o my_figures

Reads the four tables code0_*.csv ... codeC_*.csv of the run and writes seven
figures, each as PDF (vector, for a manuscript) and PNG:

    pipeline_outcomes   how every scattering and history ends, stage by stage
    pipeline_funnel     the probability of passing each stage, and of the whole chain
    pipeline_orbits     the binaries on their way through the stages
    pipeline_times      how long the stages take
    pipeline_mergers    which black-hole binaries merge
    pipeline_masses     dependence on the masses, and the masses that come out
    pipeline_spectrum   mass spectrum of the merging black holes, one generation per panel

Figures are built at the column and text widths of a two-column AASTeX
manuscript, with its fonts, so that they need no rescaling there.
Needs numpy, pandas and matplotlib.
"""
import argparse
import os

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D
from matplotlib.ticker import NullFormatter, ScalarFormatter

import make_params0 as params

# ---- units and constants (AU, yr, Msun), as in include/bhs/constants.hpp ----
PC = 206264.8
G = 4 * np.pi**2
SPEED_OF_LIGHT = 6.32397263e4
STAR_RADIUS_BEFORE = 0.00465047   # AU: of tables that do not say which radius the stars had (the Sun's)
SMBH_MASS, CUSP_SLOPE = 4.0e6, 1.75
HUBBLE_TIME = 1.38e10

# ---- figure geometry: two-column AASTeX, in inches --------------------------
PT = 1.0 / 72.27
COLUMN, PAGE = 242.26653 * PT, 513.11743 * PT

# ---- colours: Okabe-Ito, safe for colour-blind readers ----------------------
BLUE, VERMILLION, GREEN, ORANGE, PURPLE, SKY, YELLOW = (
    "#0072B2", "#D55E00", "#009E73", "#E69F00", "#CC79A7", "#56B4E9", "#F0E442")
GREY, LIGHT = "#7f7f7f", "#d0d0d0"
MOST_BINARIES = 3000   # black-hole binaries drawn one by one in a figure: a random share, if there are more

# How a scattering or history can end, from the bottom of the bar upwards:
# (legend label, outcomes of the table, colour).  The first entry of a stage
# is the one that goes on to the next stage.
STAGES = {
    "0": ("Code 0: black hole + stellar binary", "scatterings", [
        ("star captured", {"m1/m2-m3", "m1/m3-m2"}, BLUE),
        ("captured star disrupted", {"TIDAL_DISRUPTION"}, VERMILLION),
        ("stellar binary survives", {"m2/m3-m1"}, SKY),
        ("all three unbound", {"ALL_UNBOUND"}, ORANGE),
        ("unresolved", {"UNRESOLVED / RESONANT", "NUMERICAL_ERROR", "INTEGRATION_FAILED"}, LIGHT)]),
    "A": ("Code A: BH/star binary among stars", "histories", [
        ("meets a black hole", {0}, BLUE),
        ("star disrupted", {6, 1, 2}, VERMILLION),   # 1, 2: spiralling in by gravitational waves, in older tables
        ("star lost", {3, 4}, ORANGE),
        ("unresolved", {5}, LIGHT)]),
    "B": ("Code B: BH/star binary + black hole", "scatterings", [
        ("black-hole binary forms", {"m1/m3-m2", "(m1/m3)/m2"}, BLUE),
        ("star disrupted", {"TIDAL_DISRUPTION"}, VERMILLION),
        ("BH/star binary survives", {"m1/m2-m3", "(m1/m2)/m3"}, SKY),
        ("star goes to the newcomer", {"m2/m3-m1", "(m2/m3)/m1"}, GREEN),
        ("all three unbound", {"m1-m2-m3"}, ORANGE),
        ("unresolved", {"UNRESOLVED / RESONANT", "NDSOLVE_FAILED"}, LIGHT)]),
    "C": ("Code C: black-hole binary among stars and black holes", "histories", [
        ("merges", {"MERGER"}, BLUE),
        ("broken up by a black hole", {"TARGET_2G_NO_BHB_AFTER_BH:m1-m2-m3", "TARGET_2G_NO_BHB_AFTER_BH:m2/m3-m1",
                                       "TARGET_2G_NO_BHB_AFTER_BH:(m2/m3)/m1"}, ORANGE),
        ("broken up by a star", {"m1-m2-m3", "m1/m3-m2", "(m1/m3)/m2", "m2/m3-m1", "(m2/m3)/m1",
                                 "TIDAL_DISRUPTION"}, GREEN),
        ("still scattering after 1000", {"MAX_SCATTER"}, GREY),
        ("unresolved", {"STAR_UNRESOLVED_RESONANT", "BH_UNRESOLVED_RESONANT", "NDSOLVE_FAILED"}, LIGHT)]),
}


def set_style():
    plt.rcParams.update({
        "figure.constrained_layout.use": True, "savefig.bbox": None,
        "font.family": "serif", "font.serif": ["cmr10", "Computer Modern Roman", "STIXGeneral"],
        "mathtext.fontset": "cm", "axes.unicode_minus": False, "axes.formatter.use_mathtext": True,
        "font.size": 9, "axes.labelsize": 9, "axes.titlesize": 9, "legend.fontsize": 8,
        "xtick.labelsize": 8, "ytick.labelsize": 8,
        "axes.linewidth": 0.6, "lines.linewidth": 1.2, "lines.markersize": 4,
        "xtick.major.width": 0.6, "ytick.major.width": 0.6, "xtick.minor.width": 0.45, "ytick.minor.width": 0.45,
        "xtick.direction": "in", "ytick.direction": "in", "xtick.top": True, "ytick.right": True,
        "legend.frameon": False, "legend.handlelength": 1.6, "legend.borderaxespad": 0.3,
    })


# ----------------------------------------------------------------------------
# Data
# ----------------------------------------------------------------------------

def merger_time(m1, m2, a, e):
    """Gravitational-wave merger time (Peters), as mergerTime() in cluster.hpp."""
    return (5 * SPEED_OF_LIGHT**5 * a**4 / (256 * G**3 * (m1 + m2)**2 * (m1 * m2 / (m1 + m2)))
            * (1 - e * e)**3.5)


def dispersion(radius_pc):
    return np.sqrt(G * SMBH_MASS / (radius_pc * PC * (1 + CUSP_SLOPE)))


def load(run, integrator):
    """The four tables of a run, with the columns the figures need added."""
    t = {s: pd.read_csv(os.path.join(run, f"code{s}_{integrator}.csv"), comment="#") for s in "0ABC"}
    for s, table in t.items():
        column = "fate" if s == "A" else "outcome"
        for rank, (_, outcomes, _) in enumerate(STAGES[s][2]):
            table.loc[table[column].isin(outcomes), "end"] = rank
        if table["end"].isna().any():
            raise SystemExit(f"code{s}: outcomes without a category: {set(table[column][table['end'].isna()])}")
        table["passed"] = table["end"] == 0

    c0, cA, cB, cC = t["0"], t["A"], t["B"], t["C"]
    with open(os.path.join(run, f"code0_{integrator}.csv")) as handle:   # the first line says how the run was made
        made = dict(word.split("=", 1) for word in handle.readline().split() if "=" in word)
    star_radius = float(made.get("rstar", STAR_RADIUS_BEFORE))
    c0["rp_rt"] = c0.a * (1 - c0.e) / (star_radius * np.cbrt(c0.M1))      # pericentre in tidal radii
    cA["rp_rt"] = cA.a * (1 - cA.e) / (star_radius * np.cbrt(cA.M1))
    # The Code 0 capture a row descends from.  Whatever comes of one capture shares its history: its
    # Code A histories, their Code B scatterings and the black-hole binaries those form are not
    # independent draws, the captures are.
    for table in (cA, cB):
        table["seed"] = (table.kk * 100 + table.jj) * 100000 + table.h
    # Code B: speed of the star on its orbit over the typical speed of the incoming black hole
    cB["speed_ratio"] = np.sqrt(G * (cB.M1 + 1.0) / cB.a0) / dispersion(cB.R_pc)
    cC["binary"] = (cC.kk * 100 + cC.jj) * 100000 + cC.bb

    # The black-hole binaries: Code B's successes, numbered within their cell as Code C numbers them.
    bhb = cB[cB.passed].copy()
    bhb["bb"] = bhb.groupby(["kk", "jj"]).cumcount() + 1
    bhb["binary"] = (bhb.kk * 100 + bhb.jj) * 100000 + bhb.bb
    reduced = (bhb.M1 + bhb.MBH) / (bhb.M1 + bhb.MBH + 1.0)
    bhb["hardness"] = G * bhb.M1 * bhb.MBH / (bhb.a * reduced * dispersion(bhb.R_pc)**2)   # "H12s" of the notebook
    bhb["t_gw"] = merger_time(bhb.M1, bhb.MBH, bhb.a, bhb.e)
    per_binary = cC.groupby("binary").agg(f3=("passed", "mean"), histories=("passed", "size"), a_c=("a0", "first"))
    bhb = bhb.join(per_binary, on="binary")
    if len(bhb) and not np.allclose(bhb.a, bhb.a_c):
        raise SystemExit("Code C does not follow the black-hole binaries of Code B: tables of different runs?")
    t["bhb"] = bhb
    cC["seed"] = cC.binary.map(bhb.set_index("binary").seed)
    t["mergers"] = cC[cC.passed].join(bhb.set_index("binary")[["t_gw", "tb", "hardness"]], on="binary")
    return t


def proportion(passed, groups=None):
    """Fraction of True values and its standard error.  Rows of one group
    (all that descends from one captured star) are not independent draws."""
    passed = np.asarray(passed, float)
    if len(passed) == 0:
        return np.nan, np.nan
    p = passed.mean()
    if groups is None:
        return p, np.sqrt(p * (1 - p) / len(passed))
    sums = pd.DataFrame({"s": passed, "g": np.asarray(groups)}).groupby("g").s.agg(["sum", "size"])
    if len(sums) < 2:
        return p, np.nan
    residual = sums["sum"] - p * sums["size"]
    return p, np.sqrt(len(sums) / (len(sums) - 1) * np.sum(residual**2)) / len(passed)


def capture_time(m1, radius_pc, captured_fraction):
    """Mean time for a black hole to capture a star: the time between its
    encounters with stellar binaries (t_a of the Code 0 input) over the
    fraction of them that leave a star bound."""
    return np.array([params.row(0, 0, m, 0, r * PC)[-1] for m, r in zip(m1, radius_pc)]) / captured_fraction


def chain_complete(t, radius, times):
    """Probability that a single 2G black hole at this radius has captured a
    star, met a black hole, formed a binary with it and merged, all within
    each of `times`.  The capture is a Poisson process; the later stages take
    the times the simulation recorded.  One attempt per black hole."""
    p = stage_probabilities(t, lambda table: table.R_pc == radius)
    mergers = t["mergers"][t["mergers"].R_pc == radius]
    if len(mergers) == 0 or not p["0"][0] > 0:
        return np.zeros(len(times))
    wait = capture_time(mergers.M1, mergers.R_pc, p["0"][0])
    later = (mergers.tb + mergers.tm).to_numpy()
    left = np.clip(np.asarray(times)[:, None] - later[None, :], 0, None)       # time left for the capture
    return p["A"][0] * p["B"][0] * p["C"][0] * np.mean(1 - np.exp(-left / wait[None, :]), axis=1)


def stage_probabilities(t, select=lambda table: slice(None)):
    """Probability of passing each stage, and of passing all four, with
    standard errors: {"0": (p, se), ..., "all": (p, se)}.

    The four are not independent (a hard binary both forms a black-hole
    binary more easily and merges more easily), so the error of their product
    is not taken from theirs: the product is the number of mergers per Code 0
    scattering, and the scatterings are independent."""
    out = {}
    for s, groups in (("0", None), ("A", "seed"), ("B", "seed"), ("C", "seed")):
        table = t[s][select(t[s])]
        out[s] = proportion(table.passed, table[groups] if groups else None)
    p = np.prod([out[s][0] for s in "0ABC"])
    mergers = t["mergers"][select(t["mergers"])]
    per_capture = mergers.groupby("seed").size().to_numpy(float)
    scatterings = len(t["0"][select(t["0"])])
    if p > 0 and len(per_capture) > 1:
        relative = np.sqrt(max(np.sum(per_capture**2) / per_capture.sum()**2 - 1.0 / scatterings, 0.0))
    else:
        relative = np.nan
    out["all"] = (p, p * relative)
    return out


# ----------------------------------------------------------------------------
# Figure helpers
# ----------------------------------------------------------------------------

def panel_label(ax, letter, x=0.04, y=0.95, **kwargs):
    ax.text(x, y, f"({letter})", transform=ax.transAxes, va="top", ha="left", **kwargs)


def radius_colours(radii):
    return dict(zip(radii, plt.cm.viridis(np.linspace(0.0, 0.9, len(radii)))))


def radius_label(radius):
    return f"{radius:.2g}" if radius < 1 else f"{radius:.1f}"


def save(fig, out, name):
    os.makedirs(out, exist_ok=True)
    fig.savefig(os.path.join(out, name + ".pdf"), dpi=600)   # dpi: for the rasterised point clouds only
    fig.savefig(os.path.join(out, name + ".png"), dpi=300)
    plt.close(fig)
    print(f"wrote {os.path.join(out, name)}.pdf and .png")


def binned(x, y, edges, groups=None):
    """Mean of y and its standard error in bins of x: centres, means, errors."""
    centres, means, errors = [], [], []
    for lo, hi in zip(edges[:-1], edges[1:]):
        inside = (x >= lo) & (x < hi)
        if inside.sum() < 5:
            continue
        mean, error = proportion(y[inside], None if groups is None else groups[inside])
        centres.append(np.sqrt(lo * hi))
        means.append(mean)
        errors.append(error)
    return np.array(centres), np.array(means), np.array(errors)


# ----------------------------------------------------------------------------
# Figures
# ----------------------------------------------------------------------------

def fig_outcomes(t, radii, out):
    """How every scattering and history ends: one stacked bar per radius and stage."""
    fig, axes = plt.subplots(2, 2, figsize=(PAGE, 0.62 * PAGE), sharey=True)
    x = np.arange(len(radii))
    for ax, s, letter in zip(axes.flat, "0ABC", "abcd"):
        title, unit, categories = STAGES[s]
        table = t[s]
        bottom = np.zeros(len(radii))
        for rank, (label, _, colour) in enumerate(categories):
            share = np.array([np.mean(table.end[table.R_pc == r] == rank) if np.any(table.R_pc == r) else 0.0
                              for r in radii])
            ax.bar(x, share, 0.78, bottom=bottom, color=colour, edgecolor="white", linewidth=0.4, label=label)
            bottom += share
        ax.set_title(f"({letter}) {title}", loc="left")
        ax.set_xticks(x, [radius_label(r) for r in radii])
        ax.set_xlim(-0.6, len(radii) - 0.4)
        ax.set_ylim(0, 1)
        ax.tick_params(axis="x", which="both", top=False, length=0)
        handles, labels = ax.get_legend_handles_labels()
        ax.legend(handles[::-1], labels[::-1], loc="center left", bbox_to_anchor=(1.01, 0.5), handlelength=1.0,
                  title=f"{len(table):,} {unit}", title_fontsize=8, alignment="left")
    for ax in axes[1]:
        ax.set_xlabel("Distance from the central black hole (pc)")
    for ax in axes[:, 0]:
        ax.set_ylabel("Fraction")
    save(fig, out, "pipeline_outcomes")


def fig_funnel(t, radii, out):
    """The probability of passing each stage, why Code B is the narrow one, and the whole chain against time."""
    colours = radius_colours(radii)
    fig, axes = plt.subplots(1, 3, figsize=(PAGE, 0.36 * PAGE))

    # (a) Probability of passing each stage, and all four, against radius.
    ax = axes[0]
    series = [("0", "star captured (0)", SKY, "o", "-"), ("A", "meets a black hole (A)", GREEN, "s", "-"),
              ("B", "black-hole binary forms (B)", ORANGE, "^", "-"), ("C", "binary merges (C)", PURPLE, "D", "-"),
              ("all", "all four", "black", "o", "--")]
    by_radius = [stage_probabilities(t, lambda table, r=r: table.R_pc == r) for r in radii]
    for key, label, colour, marker, style in series:
        p = np.array([b[key][0] for b in by_radius])
        se = np.array([b[key][1] for b in by_radius])
        ax.errorbar(radii, p, yerr=se, color=colour, marker=marker, ls=style, mec="black", mew=0.4, capsize=1.5,
                    elinewidth=0.7, label=label)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_ylim(1e-8, 1.5)
    ax.set_xlabel("Distance from the central black hole (pc)")
    ax.set_ylabel("Probability")
    ax.set_title("(a) Passing each stage", loc="left")
    ax.legend(loc="lower left", handlelength=2.4)

    # (b) Code B: what the encounter with the black hole does, against how
    # fast the star moves on its orbit compared with the incoming black hole.
    ax = axes[1]
    cB = t["B"]
    edges = np.unique(cB.speed_ratio.quantile(np.linspace(0, 1, 9)).to_numpy())
    edges[-1] *= 1.0001
    for rank, colour, marker in ((2, SKY, "s"), (4, ORANGE, "D"), (0, BLUE, "o")):
        centre, mean, error = binned(cB.speed_ratio.to_numpy(), (cB.end == rank).to_numpy(), edges)
        ax.errorbar(centre, mean, yerr=error, color=colour, marker=marker, mec="black", mew=0.4, capsize=1.5,
                    elinewidth=0.7, label=STAGES["B"][2][rank][0])
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_ylim(3e-5, 2.0)
    ax.set_xlabel(r"Orbital speed of the star / $\sigma$")
    ax.set_ylabel("Fraction of scatterings")
    ax.set_title("(b) The encounter of Code B", loc="left")
    ax.legend(loc="lower right")

    # (c) Probability that one 2G black hole has gone through all four stages, against time.
    ax = axes[2]
    times = np.logspace(8.5, 11.5, 200)
    for r in radii:
        ax.plot(times, chain_complete(t, r, times), color=colours[r], label=radius_label(r))
    ax.axvline(HUBBLE_TIME, color="black", ls="--", lw=0.8)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(times[0], times[-1])
    ax.set_ylim(1e-8, 3e-3)
    ax.set_xlabel("Time (yr)")
    ax.set_ylabel("Probability of a 3G black hole")
    ax.set_title("(c) The whole chain", loc="left")
    ax.legend(title="radius (pc)", title_fontsize=8, loc="lower right", ncol=2, handlelength=1.2, columnspacing=0.8)
    save(fig, out, "pipeline_funnel")


def fig_orbits(t, radii, out):
    """The binaries on their way: captured, at the black-hole encounter, and as black-hole binaries."""
    colours = radius_colours(radii)
    rng = np.random.default_rng(1)
    fig, axes = plt.subplots(1, 3, figsize=(PAGE, 0.36 * PAGE))

    captured = t["0"][t["0"].end <= 1]            # with the stars disrupted at once
    reached = t["A"][t["A"].passed]
    for ax, table, letter, title in ((axes[0], captured, "a", "Stars captured (Code 0)"),
                                     (axes[1], reached, "b", "Meeting a black hole (Code A)")):
        table = table.iloc[rng.permutation(len(table))[:30000]]   # in random order: no radius hides another
        ax.scatter(table.a, table.rp_rt, s=1.5, c=[colours[r] for r in table.R_pc], alpha=0.5, linewidths=0,
                   rasterized=True)
        ax.axhline(1.0, color=VERMILLION, lw=1.0)
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlim(3e-3, 3e2)
        ax.set_ylim(0.05, 3e4)
        ax.set_xlabel("Semi-major axis (AU)")
        ax.set_title(f"({letter}) {title}", loc="left")
    axes[0].set_ylabel("Pericentre / tidal radius")
    axes[0].text(0.97, 0.115, "star disrupted", color=VERMILLION, ha="right", va="top", fontsize=8,
                 transform=axes[0].transAxes, bbox=dict(facecolor="white", edgecolor="none", alpha=0.85, pad=1.5))
    axes[1].tick_params(labelleft=False)
    handles = [Line2D([], [], marker="o", ls="", color=colours[r], ms=3.5, label=radius_label(r)) for r in radii]
    axes[0].legend(handles=handles, title="radius (pc)", title_fontsize=8, loc="upper left", ncol=2,
                   columnspacing=0.8, handletextpad=0.1)

    # Black-hole binaries, coloured by the fraction of their histories that merge.
    bhb = t["bhb"].dropna(subset=["f3"])
    shown = bhb.iloc[rng.permutation(len(bhb))[:MOST_BINARIES]].sort_values("f3")
    ax = axes[2]
    points = ax.scatter(shown.a, 1 - shown.e, s=9, c=shown.f3, cmap="viridis", vmin=0, vmax=1, edgecolors="black",
                        linewidths=0.2)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(bhb.a.min() / 2, bhb.a.max() * 2)
    ax.set_ylim(min(1e-3, (1 - bhb.e).min() / 2), 1.5)
    ax.set_xlabel("Semi-major axis (AU)")
    ax.set_ylabel(r"$1-e$")
    ax.set_title("(c) Black-hole binaries (Code B)", loc="left")
    fig.colorbar(points, ax=ax, label="Fraction of histories that merge", pad=0.02)

    # The orbits that merge by gravitational waves alone in 13.8 Gyr, for the median masses:
    # (1 - e^2)^(7/2) = time / merger_time(e = 0).  The time is written along the line.
    a = np.logspace(-2, 3.5, 400)
    one_minus_e2 = (HUBBLE_TIME / merger_time(bhb.M1.median(), bhb.MBH.median(), a, 0.0))**(2 / 7)
    one_minus_e = 1 - np.sqrt(np.clip(1 - one_minus_e2, 0, 1))
    ax.plot(a[one_minus_e2 < 1], one_minus_e[one_minus_e2 < 1], color="black", lw=0.8)
    fig.canvas.draw()   # fixes the layout, so that the slope of the line on the page is known
    i = np.searchsorted(a, 300.0)   # where no binaries are
    (x0, y0), (x1, y1) = ax.transData.transform([(a[i], one_minus_e[i]), (a[i + 5], one_minus_e[i + 5])])
    ax.text(a[i], one_minus_e[i] * 0.55, "13.8 Gyr", fontsize=8, ha="center", va="top", rotation_mode="anchor",
            rotation=np.degrees(np.arctan2(y1 - y0, x1 - x0)))
    save(fig, out, "pipeline_orbits")


def fig_times(t, radii, out):
    """How long the stages take, and how the mergers beat the gravitational-wave clock."""
    colours = radius_colours(radii)
    fig, axes = plt.subplots(1, 3, figsize=(PAGE, 0.36 * PAGE))
    mergers, reached = t["mergers"], t["A"][t["A"].passed]

    # (a) Time spent in each stage against radius: median and the range holding 68% of the cases.
    ax = axes[0]
    cells = t["0"].drop_duplicates(["kk", "jj"])
    wait = pd.Series([params.row(0, 0, m1, 0, r * PC)[-1] for m1, r in zip(cells.M1, cells.R_pc)], index=cells.R_pc)
    for values, label, colour, marker in ((wait, "to meet a stellar binary", GREY, "o"),
                                          (reached.set_index("R_pc").tb, "BH/star binary (Code A)", GREEN, "s"),
                                          (mergers.set_index("R_pc").tm, "black-hole binary (Code C)", PURPLE, "D")):
        q = np.array([values[values.index == r].quantile([0.16, 0.5, 0.84]).to_numpy() for r in radii])
        ax.fill_between(radii, q[:, 0], q[:, 2], color=colour, alpha=0.2, lw=0)
        ax.plot(radii, q[:, 1], color=colour, marker=marker, mec="black", mew=0.4, label=label)
    # The BH/star binaries that go on to form a black-hole binary are the tight ones, which wait longest.
    formers = t["bhb"].groupby("R_pc").tb.median().reindex(radii)
    ax.plot(radii, formers, color=GREEN, marker="s", mfc="white", mec=GREEN, mew=0.8, ls="--",
            label="... of those that form a binary")
    ax.axhline(HUBBLE_TIME, color="black", ls="--", lw=0.8)
    ax.text(radii[0], HUBBLE_TIME * 1.2, "age of the Universe", fontsize=8, va="bottom")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_ylim(4e6, 6e10)
    ax.set_xlabel("Distance from the central black hole (pc)")
    ax.set_ylabel("Time (yr)")
    ax.set_title("(a) Time in each stage", loc="left")
    ax.legend(loc="lower left")

    # (b) When the black-hole binaries merge, counted from their formation.
    ax = axes[1]
    for r in radii:
        times = np.sort(mergers.tm[mergers.R_pc == r])
        if len(times) >= 10:
            ax.step(times, np.arange(1, len(times) + 1) / len(times), where="post", color=colours[r],
                    label=radius_label(r))
    ax.axvline(HUBBLE_TIME, color="black", ls="--", lw=0.8)
    ax.set_xscale("log")
    ax.set_xlim(1e5, 3e11)
    ax.set_ylim(0, 1)
    ax.set_xlabel("Time from formation to merger (yr)")
    ax.set_ylabel("Fraction of mergers")
    ax.set_title("(b) Merger times", loc="left")
    ax.legend(title="radius (pc)", title_fontsize=8, loc="upper left", handlelength=1.2)

    # (c) Merger time against the time the same binary would have needed alone.
    ax = axes[2]
    shuffled = mergers.iloc[np.random.default_rng(2).permutation(len(mergers))[:30000]]
    ax.scatter(shuffled.t_gw, shuffled.tm, s=2, c=[colours[r] for r in shuffled.R_pc], alpha=0.5, linewidths=0,
               rasterized=True)
    span = [1e6, 1e21]
    ax.plot(span, span, color="black", lw=0.8, label="equal")
    ax.axvline(HUBBLE_TIME, color="black", ls="--", lw=0.8, label="13.8 Gyr")
    ax.axhline(HUBBLE_TIME, color="black", ls="--", lw=0.8)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(*span)
    ax.set_ylim(1e5, 3e11)
    ax.set_xticks([1e6, 1e9, 1e12, 1e15, 1e18, 1e21])
    ax.set_xlabel("Merger time if left alone (yr)")
    ax.set_ylabel("Merger time in the cluster (yr)")
    ax.set_title("(c) Effect of the encounters", loc="left")
    ax.legend(loc="lower right", handlelength=1.6)
    save(fig, out, "pipeline_times")


def fig_mergers(t, radii, out):
    """Which black-hole binaries merge: hardness decides, eccentricity at birth does not."""
    bhb = t["bhb"].dropna(subset=["f3"])
    histories = t["C"].join(bhb.set_index("binary")[["hardness", "e", "a"]], on="binary", rsuffix="_b")
    fig, axes = plt.subplots(1, 3, figsize=(PAGE, 0.36 * PAGE))

    # (a) Merger fraction against hardness: every binary, and the mean in bins.
    ax = axes[0]
    shown = bhb.iloc[np.random.default_rng(4).permutation(len(bhb))[:MOST_BINARIES]]
    ax.scatter(shown.hardness, shown.f3, s=4, color=LIGHT, linewidths=0, rasterized=True, label="one binary")
    edges = np.logspace(np.log10(bhb.hardness.quantile(0.01)), np.log10(bhb.hardness.quantile(0.99)), 9)
    centre, mean, error = binned(histories.hardness.to_numpy(), histories.passed.to_numpy(), edges,
                                 histories.seed.to_numpy())
    ax.errorbar(centre, mean, yerr=error, color=BLUE, marker="o", mec="black", mew=0.4, capsize=1.5, label="mean")
    ax.axvline(1.0, color="black", ls=":", lw=0.8)
    ax.set_xscale("log")
    ax.set_xlim(edges[0] / 2, edges[-1] * 2)
    ax.set_ylim(-0.03, 1.03)
    ax.set_xlabel("Hardness $H$")
    ax.set_ylabel("Fraction of histories that merge")
    ax.set_title("(a) Merger fraction", loc="left")
    ax.legend(loc="upper left", frameon=True, framealpha=0.9, edgecolor="none", borderpad=0.2)

    # (b) How the histories end, in the same bins of hardness.
    ax = axes[1]
    categories = STAGES["C"][2]
    inside = [(histories.hardness >= lo) & (histories.hardness < hi) for lo, hi in zip(edges[:-1], edges[1:])]
    bottom = np.zeros(len(inside))
    for rank, (label, _, colour) in enumerate(categories):
        share = np.array([np.mean(histories.end[sel] == rank) if sel.any() else 0.0 for sel in inside])
        ax.bar(np.sqrt(edges[:-1] * edges[1:]), share, np.diff(edges) * 0.88, bottom=bottom, color=colour,
               edgecolor="white", linewidth=0.3, label=label)
        bottom += share
    ax.set_xscale("log")
    ax.set_xlim(edges[0] / 2, edges[-1] * 2)
    ax.set_ylim(0, 1)
    ax.set_xlabel("Hardness $H$")
    ax.set_ylabel("Fraction of histories")
    ax.set_title("(b) How the histories end", loc="left")
    handles, labels = ax.get_legend_handles_labels()
    fig.legend(handles[::-1], labels[::-1], loc="outside lower center", ncol=5, handlelength=1.0, columnspacing=1.2)

    # (c) The eccentricity the binary is born with makes no difference.
    ax = axes[2]
    one_minus_e = (1 - histories.e).to_numpy()
    for (lo, hi), colour, marker in (((0, 1), ORANGE, "s"), ((1, 10), GREEN, "o"), ((10, np.inf), BLUE, "D")):
        sel = ((histories.hardness >= lo) & (histories.hardness < hi)).to_numpy()
        centre, mean, error = binned(one_minus_e[sel], histories.passed.to_numpy()[sel], np.logspace(-2.5, 0, 6),
                                     histories.seed.to_numpy()[sel])
        label = f"$H<{hi}$" if lo == 0 else f"${lo}<H<{hi}$" if np.isfinite(hi) else f"$H>{lo}$"
        ax.errorbar(centre, mean, yerr=error, color=colour, marker=marker, mec="black", mew=0.4, capsize=1.5,
                    label=label)
    ax.set_xscale("log")
    ax.set_ylim(-0.03, 1.22)
    ax.set_yticks(np.arange(0, 1.01, 0.2))
    ax.set_xlabel(r"$1-e$ at formation")
    ax.set_ylabel("Fraction of histories that merge")
    ax.set_title("(c) Eccentricity at formation", loc="left")
    ax.legend(loc="upper center", ncol=3, handlelength=1.0, columnspacing=0.8, handletextpad=0.3)
    save(fig, out, "pipeline_mergers")


def fig_masses(t, radii, out):
    """Dependence on the masses, and the masses that come out."""
    fig, axes = plt.subplots(1, 3, figsize=(PAGE, 0.36 * PAGE))
    series = [("0", "star captured", SKY, "o"), ("A", "meets a black hole", GREEN, "s"),
              ("B", "black-hole binary forms", ORANGE, "^"), ("C", "binary merges", PURPLE, "D")]

    # (a), (b): probability of passing each stage against the two masses.
    for ax, column, title, label, ticks in (
            (axes[0], "M1", "(a) Against the 2G mass", r"2G black hole $M_1$ ($M_\odot$)", [10, 20, 30, 50, 70]),
            (axes[1], "MBH", "(b) Against the 1G mass", r"1G black hole $M_{\rm BH}$ ($M_\odot$)", [5, 10, 20, 40])):
        masses = t["0"].drop_duplicates("kk")[column]
        edges = np.unique(masses.quantile(np.linspace(0, 1, 6 if len(masses) >= 30 else 3)).to_numpy())
        edges[-1] *= 1.0001
        for s, name, colour, marker in series:
            groups = {"A": "seed", "B": "seed", "C": "seed"}.get(s)
            table = t[s]
            centre, mean, error = binned(table[column].to_numpy(), table.passed.to_numpy(), edges,
                                         table[groups].to_numpy() if groups else None)
            ax.errorbar(centre, mean, yerr=error, color=colour, marker=marker, mec="black", mew=0.4, capsize=1.5,
                        label=name)
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xticks(ticks)
        ax.xaxis.set_major_formatter(ScalarFormatter())
        ax.xaxis.set_minor_formatter(NullFormatter())
        ax.set_ylim(1e-3, 1.0)
        ax.set_xlabel(label)
        ax.set_title(title, loc="left")
    axes[0].set_ylabel("Probability of passing the stage")
    axes[1].tick_params(labelleft=False)
    fig.legend(*axes[0].get_legend_handles_labels(), loc="outside lower left", ncol=2, columnspacing=1.2)

    # (c) The masses going in and the mass coming out.
    ax = axes[2]
    pairs, mergers = t["0"].drop_duplicates("kk"), t["mergers"]
    bins = np.logspace(np.log10(4), np.log10(160), 33)
    for values, label, colour, style in ((pairs.MBH, "1G partners (input)", GREY, ":"),
                                         (pairs.M1, "2G black holes (input)", GREEN, "--"),
                                         (mergers.m2_merge, "partner at merger", ORANGE, "-"),
                                         (mergers.m_3g, "3G remnants", BLUE, "-")):
        if len(values):
            ax.hist(values, bins=bins, weights=np.full(len(values), 1 / len(values)), histtype="step", color=colour,
                    ls=style, lw=1.2, label=label)
    ax.set_xscale("log")
    ax.set_xlabel(r"Mass ($M_\odot$)")
    ax.set_ylabel("Fraction per bin")
    ax.set_title("(c) Masses in and out", loc="left")
    fig.legend(*ax.get_legend_handles_labels(), loc="outside lower right", ncol=2, handlelength=1.4,
               columnspacing=1.2)
    save(fig, out, "pipeline_masses")


def spectrum(values, groups, bins, draws=400):
    """Mass spectrum dP / dlog10(M) in `bins`, with the range that holds 68%
    of the spectra obtained when the groups (the mass pairs) are drawn again
    with replacement: values of one mass pair are not independent."""
    width = np.diff(np.log10(bins))
    keys, index = np.unique(groups, return_inverse=True)
    counts = np.zeros((len(keys), len(bins) - 1))
    np.add.at(counts, (index, np.clip(np.digitize(values, bins) - 1, 0, len(bins) - 2)), 1.0)
    density = counts.sum(axis=0) / (counts.sum() * width)
    rng = np.random.default_rng(3)
    resampled = np.array([counts[rng.integers(0, len(keys), len(keys))].sum(axis=0) for _ in range(draws)])
    resampled = resampled / (resampled.sum(axis=1, keepdims=True) * width)
    low, high = np.percentile(resampled, [16, 84], axis=0)
    return density, low, high


def fig_spectrum(t, radii, out):
    """Mass spectrum of the black holes in the mergers: 1G partner, 2G black hole, 3G remnant."""
    pairs, mergers = t["0"].drop_duplicates("kk"), t["mergers"]
    bins = np.logspace(np.log10(4), np.log10(200), 27)
    width = np.diff(np.log10(bins))
    panels = [("(a) 1G: the partner", pairs.MBH, mergers.m2_merge, ORANGE),
              ("(b) 2G: the black hole followed", pairs.M1, mergers.M1, GREEN),
              ("(c) 3G: the remnant", 0.95 * (pairs.M1 + pairs.MBH), mergers.m_3g, BLUE)]

    fig, axes = plt.subplots(1, 3, figsize=(PAGE, 0.36 * PAGE), sharey=True)
    for ax, (title, everyone, merging, colour) in zip(axes, panels):
        reference = np.histogram(everyone, bins)[0] / (len(everyone) * width)
        ax.stairs(reference, bins, color=GREY, lw=1.0, ls="--",
                  label=f"all mass pairs: median {np.median(everyone):.1f} $M_\\odot$")
        density, low, high = spectrum(merging.to_numpy(), mergers.kk.to_numpy(), bins)
        ax.stairs(high, bins, baseline=low, fill=True, color=colour, alpha=0.3, lw=0)
        ax.stairs(density, bins, color=colour, lw=1.3,
                  label=f"in the mergers: median {np.median(merging):.1f} $M_\\odot$")
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlim(bins[0], bins[-1])
        ax.set_xticks([5, 10, 20, 50, 100, 200])
        ax.xaxis.set_major_formatter(ScalarFormatter())
        ax.xaxis.set_minor_formatter(NullFormatter())
        ax.set_xlabel(r"Mass ($M_\odot$)")
        ax.set_title(title, loc="left")

    # The partners that were not there from the start but took the place of the first one in Code C.
    exchanged = mergers.n_exch > 0
    share = np.histogram(mergers.m2_merge[exchanged], bins)[0] / (len(mergers) * width)
    axes[0].stairs(share, bins, color="black", lw=0.9, ls=":",
                   label=f"gained by exchange: {exchanged.mean():.1%}")
    axes[0].set_ylabel(r"$\mathrm{d}P/\mathrm{d}\log_{10}M$")
    axes[0].set_ylim(0.01, 60)   # the top decade is for the legends
    for ax in axes:
        ax.legend(loc="upper left", handlelength=1.6)
    save(fig, out, "pipeline_spectrum")


# ----------------------------------------------------------------------------
# Numbers
# ----------------------------------------------------------------------------

def print_summary(t, radii):
    c0, cA, cC, bhb, mergers = t["0"], t["A"], t["C"], t["bhb"], t["mergers"]
    print(f"\n{c0.kk.nunique()} mass pairs at {len(radii)} radii")
    print(f"{'':28s}" + "".join(f"{'R = ' + radius_label(r) + ' pc':>16s}" for r in radii) + f"{'all radii':>16s}")
    rows = [stage_probabilities(t, lambda table, r=r: table.R_pc == r) for r in radii] + [stage_probabilities(t)]
    for key, name in (("0", "star captured"), ("A", "meets a black hole"), ("B", "black-hole binary forms"),
                      ("C", "binary merges"), ("all", "all four")):
        print(f"{name:28s}" + "".join(f"{row[key][0]:9.2e} ±{row[key][1]:5.0e}" for row in rows))
    done = [chain_complete(t, r, [HUBBLE_TIME])[0] for r in radii]
    print(f"{'3G within 13.8 Gyr':28s}" + "".join(f"{p:16.2e}" for p in done))
    wait = [np.median(capture_time(c0.M1[c0.R_pc == r].unique(), np.full(c0.kk.nunique(), r),
                                   stage_probabilities(t, lambda table, r=r: table.R_pc == r)["0"][0])) for r in radii]
    print(f"{'capture time, median (yr)':28s}" + "".join(f"{w:16.2e}" for w in wait))
    for s in "0ABC":
        table = t[s]
        shares = ", ".join(f"{label} {np.mean(table.end == rank):.1%}" for rank, (label, _, _) in enumerate(STAGES[s][2]))
        print(f"\ncode{s}: {len(table):,} {STAGES[s][1]}: {shares}")

    captured, reached = c0[c0.end <= 1], cA[cA.passed]
    print(f"\ncaptured stars: e median {captured.e.median():.2f}; pericentre median {captured.rp_rt.median():.1f} tidal"
          f" radii; {np.mean(captured.end == 1):.1%} inside the tidal radius at once")
    print(f"Code A: stellar encounters before disruption, median {cA.n_enc[cA.end == 1].median():.0f}; before the"
          f" black hole arrives, median {reached.n_enc.median():.0f}")
    growth = (reached.a / reached.a0).quantile([0.16, 0.5, 0.84]).to_numpy()
    print(f"BH/star binaries meeting a black hole: a / a0 = {growth[1]:.1f} ({growth[0]:.1f} to {growth[2]:.1f});"
          f" e median {reached.e.median():.2f}, from {reached.e0.median():.2f}")
    ratio = t["B"].speed_ratio.quantile([0.16, 0.5, 0.84]).to_numpy()
    print(f"Code B: orbital speed of the star / sigma = {ratio[1]:.2f} ({ratio[0]:.2f} to {ratio[2]:.2f});"
          f" black-hole binary forms in {t['B'].passed[t['B'].speed_ratio < 1].mean():.2%} of scatterings below 1,"
          f" {t['B'].passed[t['B'].speed_ratio >= 1].mean():.2%} above")
    print(f"black-hole binaries: {len(bhb)}; hard (hardness > 1): {np.mean(bhb.hardness > 1):.1%};"
          f" would merge alone within 13.8 Gyr: {np.mean(bhb.t_gw < HUBBLE_TIME):.1%}")
    for lo, hi in ((0, 1), (1, 10), (10, np.inf)):
        sel = cC.binary.isin(bhb.binary[(bhb.hardness >= lo) & (bhb.hardness < hi)])
        p, se = proportion(cC.passed[sel], cC.seed[sel])
        print(f"  hardness {lo} to {hi}: {sel.sum() // max(cC.groupby('binary').size().max(), 1)} binaries,"
              f" merge in {p:.1%} ± {se:.1%} of histories")
    if len(mergers):
        ratio = np.log10(mergers.tm / mergers.t_gw).quantile([0.16, 0.5, 0.84]).to_numpy()
        print(f"mergers: {len(mergers):,}; time from formation, median {mergers.tm.median():.2e} yr;"
              f" within 13.8 Gyr: {np.mean(mergers.tm < HUBBLE_TIME):.1%}")
        print(f"  the same binaries left alone: median {mergers.t_gw.median():.2e} yr; longer than 13.8 Gyr for"
              f" {np.mean(mergers.t_gw > HUBBLE_TIME):.1%}; log10(time in cluster / time alone) = {ratio[1]:.1f}"
              f" ({ratio[0]:.1f} to {ratio[2]:.1f})")
        print(f"  stellar encounters before merger, median {mergers.n_star.median():.0f}; black-hole encounters,"
              f" median {mergers.n_bh.median():.0f}; with a partner gained by exchange: {np.mean(mergers.n_exch > 0):.1%}")
        pairs = c0.drop_duplicates("kk")
        print(f"  masses in the mergers: 2G median {mergers.M1.median():.1f} Msun (all pairs {pairs.M1.median():.1f});"
              f" partner median {mergers.m2_merge.median():.1f} Msun (all pairs {pairs.MBH.median():.1f})")
        for name, values, everyone in (("1G partner", mergers.m2_merge, pairs.MBH), ("2G black hole", mergers.M1, pairs.M1),
                                       ("3G remnant", mergers.m_3g, 0.95 * (pairs.M1 + pairs.MBH))):
            q, q0 = values.quantile([0.05, 0.5, 0.95]).to_numpy(), everyone.quantile([0.05, 0.5, 0.95]).to_numpy()
            print(f"  {name:14s} in the mergers: median {q[1]:5.1f} Msun, 90% between {q[0]:5.1f} and {q[2]:5.1f};"
                  f" all mass pairs: median {q0[1]:5.1f}, {q0[0]:5.1f} to {q0[2]:5.1f}")
        ratio = (mergers.m2_merge / mergers.M1).clip(upper=None)
        ratio = np.minimum(ratio, 1 / ratio)
        print(f"  mass ratio (lighter / heavier) of the mergers: median {ratio.median():.2f}; partner heavier than"
              f" the 2G black hole in {np.mean(mergers.m2_merge > mergers.M1):.1%}; mass pairs with a merger:"
              f" {mergers.kk.nunique()}")
        print(f"  3G remnant mass: median {mergers.m_3g.median():.1f} Msun, 90% between"
              f" {mergers.m_3g.quantile(0.05):.1f} and {mergers.m_3g.quantile(0.95):.1f}")
        chain = mergers.tb + mergers.tm
        print(f"  Code A time of the BH/star binaries that form a black-hole binary: median {bhb.tb.median():.2e} yr"
              f" (all that meet a black hole: {reached.tb.median():.2e}); Code A + Code C for the mergers: median"
              f" {chain.median():.2e} yr, under 13.8 Gyr for {np.mean(chain < HUBBLE_TIME):.1%}")
        print("  median time by radius (yr): Code A " + " ".join(f"{reached.tb[reached.R_pc == r].median():.1e}" for r in radii)
              + " | Code C " + " ".join(f"{mergers.tm[mergers.R_pc == r].median():.1e}" for r in radii))


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("run", help="directory with the four tables of a pipeline run")
    parser.add_argument("-o", "--out", default="docs/figures", help="where the figures go [docs/figures]")
    parser.add_argument("--integrator", default="spacehub", help="which build wrote the tables [spacehub]")
    args = parser.parse_args()

    set_style()
    t = load(args.run, args.integrator)
    radii = np.sort(t["0"].R_pc.unique())
    print_summary(t, radii)
    print()
    for figure in (fig_outcomes, fig_funnel, fig_orbits, fig_times, fig_mergers, fig_masses, fig_spectrum):
        figure(t, radii, args.out)


if __name__ == "__main__":
    main()
