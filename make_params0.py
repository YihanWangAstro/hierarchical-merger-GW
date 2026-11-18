#!/usr/bin/env python3
"""Write the input table of Code 0 for a list of mass pairs.

    python make_params0.py                                   # data/M12S.txt -> data/params0_1000.csv
    python make_params0.py masses.txt -o params0_mine.csv
    python make_params0.py --notebook "2 Code 0 - Yihan.txt" # the table as a Code 0 notebook holds it

The mass list has one pair per line: the 2G black hole M1 and the 1G black
hole MBH it meets in Code B, in solar masses. Every pair is placed at the six
radii of the notebook "2 Code 0.nb", with the stellar binary, its breakup
radius and the encounter time that the notebook takes at that radius:

    a_ss   = 6.875e-7 R                      semi-major axis of the stellar binary, e_ss = 0
    r_tide = a_ss (M1 / 2 m)^(1/3)           where the black hole's tide splits it
    t_a    = 10 / (2 pi G (M1 + 2 m) n(R) / sigma(R) * 2 r_tide)

with m = 1 Msun and n, sigma the density and velocity dispersion of the stars.
These reproduce the twelve rows of data/params0.csv to 16 digits.

With --notebook nothing is computed: the lists M12S, Routs, ABIs, EBIs, Rps
and TA are taken as they stand in the text copy of a Code 0 notebook (made
with nb2txt.py), which is the way to run exactly what a notebook runs.
"""
import argparse
import ast
import math
import re

PC = 206264.8                    # AU
G = 4 * math.pi**2               # AU^3 / (Msun yr^2)
STAR_MASS = 1.0
SMBH_MASS = 4.0e6
CUSP_SLOPE = 1.75
RADII = [21732.271304079746, 39153.908143096894, 70541.56932921596,     # AU: the six outer radii
         127091.09356956737, 228973.44386154477, 412529.59999999986]    # of the notebooks' grid
BINARY_SIZE = 6.875e-7           # a_ss / R
TIME_FACTOR = 10.0


def star_density(radius):
    return 1.35e6 / (STAR_MASS * PC**3) * (1 / 0.25) ** -CUSP_SLOPE * (radius / PC) ** -CUSP_SLOPE


def dispersion(radius):
    return math.sqrt(G * SMBH_MASS / (radius * (1 + CUSP_SLOPE)))


def row(kk, jj, m1, mbh, radius):
    a_ss = BINARY_SIZE * radius
    r_tide = a_ss * (m1 / (2 * STAR_MASS)) ** (1 / 3)
    rate = 2 * math.pi * G * (m1 + 2 * STAR_MASS) * star_density(radius) / dispersion(radius) * 2 * r_tide
    return [kk, jj, m1, mbh, radius, a_ss, 0, r_tide, TIME_FACTOR / rate]


def notebook_rows(path):
    """The rows of the table from the lists of a Code 0 notebook, read from its text copy."""
    with open(path) as handle:
        text = handle.read()

    def numbers(name, like=None):
        """The list assigned to `name`; if it is not written out as numbers, zeros in the shape of `like`."""
        found = re.search(rf"^{name}=(\{{.*?\}});?\s*$", text, re.M)
        if not found:
            if like is None:
                raise SystemExit(f"error: {path}: no list {name}=...")
            return [0.0] * len(like)
        literal = found.group(1).replace("`", "").replace("*^", "e").replace("{", "[").replace("}", "]")
        return ast.literal_eval(literal)

    masses, radii, a_ss = numbers("M12S"), numbers("Routs"), numbers("ABIs")
    e_ss, breakup, time = numbers("EBIs", like=a_ss), numbers("Rps"), numbers("TA")
    if not (len(breakup) == len(time) == len(masses) and all(len(r) == len(radii) for r in breakup + time)
            and len(a_ss) == len(e_ss) == len(radii)):
        raise SystemExit(f"error: {path}: the lists do not have the shapes of one grid")
    return [[kk + 1, jj + 1, float(masses[kk][0]), float(masses[kk][1]), float(radii[jj]), float(a_ss[jj]),
             float(e_ss[jj]), float(breakup[kk][jj]), float(time[kk][jj])]
            for kk in range(len(masses)) for jj in range(len(radii))], len(masses), len(radii)


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("masses", nargs="?", default="data/M12S.txt", help="mass pairs, one per line: M1 MBH")
    parser.add_argument("--notebook", help="take the table from the text copy of a Code 0 notebook instead")
    parser.add_argument("-o", "--out", default="data/params0_1000.csv")
    args = parser.parse_args()

    if args.notebook:
        rows, num_pairs, num_radii = notebook_rows(args.notebook)
    else:
        with open(args.masses) as handle:
            pairs = [tuple(float(x) for x in line.split()) for line in handle if line.strip()]
        rows = [row(kk, jj, m1, mbh, radius) for kk, (m1, mbh) in enumerate(pairs, 1)
                for jj, radius in enumerate(RADII, 1)]
        num_pairs, num_radii = len(pairs), len(RADII)
    with open(args.out, "w") as out:
        out.write("kk,jj,M1,MBH,R_AU,a_ss,e_ss,r_tide,t_a\n")
        for values in rows:
            out.write(",".join(repr(v) for v in values) + "\n")
    print(f"{args.out}: {num_pairs} mass pairs x {num_radii} radii")


if __name__ == "__main__":
    main()
