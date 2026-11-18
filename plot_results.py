#!/usr/bin/env python3
"""Plots for one run of the BH/star binary scattering code.

    python plot_results.py profile                       # cluster model (the notebook's 3 panels)
    python plot_results.py results codeA_bs.csv          # outcome plots of a Code A run
    python plot_results.py traj traj_bs.csv              # one three-body encounter

Each command writes a PNG (-o to choose the name, --show to open a window).
To compare the three integrators use compare.py.
"""
import argparse

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import Normalize

PC = 206264.8                               # AU
G = 4 * np.pi**2                            # AU^3 / (Msun yr^2)
AU_YR_TO_KMS = 1.5e8 / (365 * 24 * 3600)    # same conversion as the notebook

# Okabe-Ito colours (colour-blind safe).
BLUE, VERMILLION, GREEN, ORANGE, PURPLE, SKY, GREY = (
    "#0072B2", "#D55E00", "#009E73", "#E69F00", "#CC79A7", "#56B4E9", "#8c8c8c")

# Fate codes written by the C++ code (enum Fate), with a colour each.  1 and 2 are in tables written
# before 17 November 2026 only: such a history now ends as 6, the star being disrupted on its way in.
FATES = [
    ("BH encounter reached", BLUE),
    ("merger before BH encounter", ORANGE),
    ("merger before next star", VERMILLION),
    ("BH lost its companion", GREEN),
    ("ionised", PURPLE),
    ("unresolved", GREY),
    ("star tidally disrupted", SKY),
]


def set_style():
    plt.rcParams.update({
        "font.size": 10, "axes.labelsize": 10, "legend.fontsize": 8.5,
        "xtick.labelsize": 9, "ytick.labelsize": 9,
        "xtick.direction": "in", "ytick.direction": "in", "xtick.top": True, "ytick.right": True,
        "axes.linewidth": 0.8, "lines.linewidth": 1.4, "legend.frameon": False,
        "savefig.dpi": 200, "figure.constrained_layout.use": True,
    })


def load_csv(path):
    """Read an output file of the C++ code: returns (structured array, metadata dict).

    The first line is '# key=value key=value ...', the second the column names.
    """
    with open(path) as fh:
        first = fh.readline()
    meta = dict(kv.split("=", 1) for kv in first.lstrip("# ").split() if "=" in kv)
    data = np.genfromtxt(path, delimiter=",", names=True, skip_header=1, dtype=None, encoding="utf-8")
    return data, meta


def panel_label(ax, letter):
    ax.text(0.0, 1.03, letter, transform=ax.transAxes, fontweight="bold", fontsize=11, va="bottom", ha="left")


def finish(fig, out, show):
    fig.savefig(out)
    print("wrote", out)
    if show:
        plt.show()


# ----------------------------------------------------------------------------
def cmd_profile(args):
    msmbh, gamma, f_bh = 4e6, 1.75, 0.01
    r = np.exp(np.linspace(np.log(0.01 * PC), np.log(2 * PC), 10))
    n_star = 1.35e6 / PC**3 * 0.25**gamma * (r / PC) ** (-gamma)
    n_bh = f_bh * n_star
    sigma = np.sqrt(G * msmbh / (r * (1 + gamma)))
    rpc = r / PC

    fig, ax = plt.subplots(1, 3, figsize=(11, 3.2))
    ax[0].loglog(rpc, n_star * PC**3, "k-", label="stars")
    ax[0].loglog(rpc, n_bh * PC**3, "k-.", label="black holes")
    ax[0].set_ylabel(r"$n$ (pc$^{-3}$)")
    ax[0].legend()
    ax[1].loglog(rpc, 4 * np.pi * n_star * r**3, "k-")
    ax[1].loglog(rpc, 4 * np.pi * n_bh * r**3, "k-.")
    ax[1].set_ylabel(r"d$N$/dln$r$")
    ax[2].loglog(rpc, sigma * AU_YR_TO_KMS, "k-")
    ax[2].set_ylabel(r"$\sigma$ (km s$^{-1}$)")
    for panel, letter in zip(ax, "ABC"):
        panel.set_xlabel(r"$R$ (pc)")
        panel.set_xlim(0.01, 2)
        panel_label(panel, letter)
    finish(fig, args.out or "profile.png", args.show)


# ----------------------------------------------------------------------------
def cmd_results(args):
    d, meta = load_csv(args.file)
    radii = np.unique(d["R_pc"])
    fate = d["fate"]
    print(f"{len(d)} realisations, {len(np.unique(d['kk']))} systems, {len(radii)} radii"
          f" (integrator: {meta.get('integrator', '?')})")
    for i, (name, _) in enumerate(FATES):
        print(f"  {name:28s} {np.mean(fate == i):6.1%}")

    fig, axes = plt.subplots(2, 2, figsize=(10, 7.5))

    # A: outcome fractions against radius
    ax = axes[0, 0]
    bottom = np.zeros(len(radii))
    xs = np.arange(len(radii))
    for i, (name, col) in enumerate(FATES):
        frac = np.array([np.mean(fate[d["R_pc"] == r] == i) for r in radii])
        if frac.sum() == 0:
            continue
        ax.bar(xs, frac, bottom=bottom, color=col, label=name, width=0.8, edgecolor="white", linewidth=0.5)
        bottom += frac
    ax.set_xticks(xs)
    ax.set_xticklabels([f"{r:.2g}" for r in radii], fontsize=7.5)
    ax.set_xlabel(r"$R$ (pc)")
    ax.set_ylabel("fraction of realisations")
    ax.set_ylim(0, 1)
    ax.tick_params(axis="x", which="both", top=False)
    fig.legend(*ax.get_legend_handles_labels(), loc="outside upper center", ncol=3)

    reached = fate == 0
    cmap = plt.get_cmap("viridis")
    log_r = np.log10(d["R_pc"])
    norm = Normalize(log_r.min(), log_r.max())

    # B: a - e of the binaries that reach the BH encounter, with their starting points
    ax = axes[0, 1]
    ax.scatter(d["a0"], 1 - d["e0"], s=6, c="0.75", label="initial", rasterized=True)
    sc = ax.scatter(d["a"][reached], 1 - d["e"][reached], s=4, c=log_r[reached], cmap=cmap, norm=norm, alpha=0.5,
                    label="at BH encounter", rasterized=True)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel(r"$a$ (AU)")
    ax.set_ylabel(r"$1-e$")
    ax.legend(loc="lower left")
    fig.colorbar(sc, ax=ax, label=r"log$_{10}$ $R$ (pc)")

    # C: change in semi-major axis, one curve per radius
    ax = axes[1, 0]
    for r in radii:
        m = reached & (d["R_pc"] == r)
        if m.sum() < 5:
            continue
        x = np.sort(np.log10(d["a"][m] / d["a0"][m]))
        ax.plot(x, np.arange(1, len(x) + 1) / len(x), color=cmap(norm(np.log10(r))), lw=1.2)
    ax.axvline(0, color="0.6", lw=0.8, ls=":")
    ax.set_xlabel(r"log$_{10}$ ($a/a_0$) at BH encounter")
    ax.set_ylabel("cumulative fraction")

    # D: time to the BH encounter
    ax = axes[1, 1]
    quantiles = np.array([np.percentile(d["tb"][reached & (d["R_pc"] == r)], [16, 50, 84])
                          if np.any(reached & (d["R_pc"] == r)) else [np.nan] * 3 for r in radii])
    ax.fill_between(radii, quantiles[:, 0], quantiles[:, 2], color=BLUE, alpha=0.25, lw=0, label="16-84%")
    ax.plot(radii, quantiles[:, 1], "o-", color=BLUE, label="median")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel(r"$R$ (pc)")
    ax.set_ylabel("time to BH encounter (yr)")
    ax.legend()

    for ax, letter in zip(axes.flat, "ABCD"):
        panel_label(ax, letter)
    finish(fig, args.out or "results.png", args.show)


# ----------------------------------------------------------------------------
def cmd_traj(args):
    d, meta = load_csv(args.file)
    m1, a0 = float(meta.get("M1", 10)), float(meta.get("a0", 1))
    m2 = 1.0
    r = [np.column_stack([d[f"{c}{i}"] for c in "xyz"]) for i in (1, 2, 3)]
    com = (m1 * r[0] + m2 * r[1]) / (m1 + m2)   # centre of mass of the original binary
    r = [x - com for x in r]
    colours = [VERMILLION, BLUE, "black"]
    labels = ["m1 (black hole)", "m2 (companion star)", "m3 (incoming star)"]

    fig = plt.figure(figsize=(11, 5))
    ax3 = fig.add_subplot(1, 2, 1, projection="3d")
    for x, c, lab in zip(r, colours, labels):
        ax3.plot(x[:, 0] / a0, x[:, 1] / a0, x[:, 2] / a0, color=c, lw=0.6, label=lab)
        ax3.scatter(*(x[-1] / a0), color=c, s=15)
    lim = args.zoom
    ax3.set(xlim=(-lim, lim), ylim=(-lim, lim), zlim=(-lim, lim),
            xlabel=r"$x/a_0$", ylabel=r"$y/a_0$", zlabel=r"$z/a_0$")
    ax3.set_title("frame of the original binary", fontsize=10)
    ax3.legend(loc="upper left")

    ax = fig.add_subplot(1, 2, 2)
    for (i, j), c, lab in zip([(0, 1), (0, 2), (1, 2)], [VERMILLION, "black", BLUE], ["m1-m2", "m1-m3", "m2-m3"]):
        ax.semilogy(d["t"], np.linalg.norm(r[i] - r[j], axis=1) / a0, color=c, lw=0.8, label=lab)
    ax.set_xlabel(r"$t$ (yr)")
    ax.set_ylabel(r"separation / $a_0$")
    ax.legend()
    ax.set_title(f"{meta.get('integrator', '')}: outcome {meta.get('outcome', '?')}", fontsize=10)
    fig.suptitle(f"encounter {meta.get('encounter', '?')}:  M1 = {m1:.3g} Msun,  a0 = {a0:.4g} AU,  "
                 f"e0 = {float(meta.get('e0', 0)):.4f},  R = {float(meta.get('R_pc', 0)):.3g} pc", fontsize=10)
    finish(fig, args.out or "traj.png", args.show)


# ----------------------------------------------------------------------------
def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = parser.add_subparsers(dest="cmd", required=True)
    common = argparse.ArgumentParser(add_help=False)
    common.add_argument("-o", "--out", help="output image file")
    common.add_argument("--show", action="store_true", help="also open a window")
    sub.add_parser("profile", parents=[common]).set_defaults(func=cmd_profile)
    s = sub.add_parser("results", parents=[common])
    s.add_argument("file")
    s.set_defaults(func=cmd_results)
    s = sub.add_parser("traj", parents=[common])
    s.add_argument("file")
    s.add_argument("--zoom", type=float, default=5.0, help="half-width of the 3D box in units of a0")
    s.set_defaults(func=cmd_traj)
    args = parser.parse_args()
    set_style()
    args.func(args)


if __name__ == "__main__":
    main()
