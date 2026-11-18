#!/usr/bin/env python3
"""Write the a, e, t of a Code A run on the original input as the three nested lists of `Code A.nb`.

    python export_aet.py codeA_spacehub.csv                # -> codeA_spacehub.txt
    python export_aet.py out/codeA_*.csv                   # one .txt per input
    python export_aet.py codeA_bs.csv -o mine.txt

This is for runs of `codeA --in data/bhstar_tidal_capture.csv`, where every
system has one binary per radius, as in the original notebook. A run of the
pipeline has several binaries per system and radius, which these lists cannot
hold: its a, e and tb are in the table itself, one row per history.

The file holds three lines of Mathematica input, as the notebook's analysis
cells expect them:

    AINSS={{{...},...},...};     semi-major axis at the black-hole encounter [AU]
    EINSS={{{...},...},...};     eccentricity at the black-hole encounter
    TBsS={{{...},...},...};      elapsed time at the black-hole encounter [yr]

Each list is indexed [[system]][[radius]][[realisation]] (kk, jj, ii). An entry
is 0 where the binary never reached the black-hole encounter; the analysis
drops those with DeleteCases[..., 0] and counts the rest.

Paste the file into a notebook cell, or load it with Get["codeA_spacehub.txt"].
"""
import argparse
import csv
import os
import sys
from decimal import Decimal

# CSV column -> notebook variable
COLUMNS = (("a", "AINSS"), ("e", "EINSS"), ("tb", "TBsS"))


def mathematica_number(text):
    """One CSV number written the way Mathematica writes a machine real.

    0.31735951552887776`  and  6.937361306167924`*^8 : positional below 10^6,
    scientific from there on and below 10^-5. Without the backtick Mathematica
    would read a 17-digit literal as an arbitrary-precision number. Entries
    that were never set stay the integer 0, as in the notebook.
    """
    value = float(text)
    if value == 0.0:
        return "0"
    number = Decimal(repr(value)).normalize()
    sign, digits, _ = number.as_tuple()
    power = number.adjusted()                           # value = d.ddd x 10^power
    if -5 <= power < 6:
        literal = format(abs(number), "f")
        if "." not in literal:
            literal += "."
        suffix = "`"
    else:
        literal = f"{digits[0]}." + "".join(map(str, digits[1:]))
        suffix = f"`*^{power}"
    return ("-" if sign else "") + literal + suffix


def read_results(path):
    """Return {column: {kk: {jj: {ii: literal}}}} from a codeA_*.csv file."""
    tables = {column: {} for column, _ in COLUMNS}
    with open(path, newline="") as handle:
        for row in csv.DictReader(line for line in handle if not line.startswith("#")):
            if int(row.get("h", 1)) != 1:
                sys.exit(f"error: {path}: several binaries per system and radius (column h); "
                         "the three lists hold one each")
            kk, jj, ii = int(row["kk"]), int(row["jj"]), int(row["ii"])
            for column, _ in COLUMNS:
                tables[column].setdefault(kk, {}).setdefault(jj, {})[ii] = mathematica_number(row[column])
    return tables


def nested_list(table, path):
    """{kk: {jj: {ii: literal}}} -> ('{{{...},...},...}', shape); every system must have the same shape."""
    if not table:
        sys.exit(f"error: {path}: no rows found")

    def complete(mapping, label):
        keys = sorted(mapping)
        if keys != list(range(1, len(keys) + 1)):
            sys.exit(f"error: {path}: {label} is not numbered 1..{len(keys)} without gaps")
        return keys

    first_system = table[min(table)]
    shape = (len(first_system), len(first_system[min(first_system)]))
    systems = []
    for kk in sorted(table):
        radii = []
        for jj in complete(table[kk], f"jj of system {kk}"):
            row = [table[kk][jj][ii] for ii in complete(table[kk][jj], f"ii of system {kk}, radius {jj}")]
            if (len(table[kk]), len(row)) != shape:
                sys.exit(f"error: {path}: system {kk}, radius {jj} has a different size from the rest")
            radii.append("{" + ",".join(row) + "}")
        systems.append("{" + ",".join(radii) + "}")
    return "{" + ",".join(systems) + "}", (len(systems), *shape)


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("csv", nargs="+", help="codeA_*.csv of a run on the original Code A input")
    parser.add_argument("-o", "--out", help="output file (only with a single input)")
    args = parser.parse_args()
    if args.out and len(args.csv) > 1:
        parser.error("--out needs a single input file")

    for path in args.csv:
        tables = read_results(path)
        lists = [(name, *nested_list(tables[column], path)) for column, name in COLUMNS]
        out_path = args.out or os.path.splitext(path)[0] + ".txt"
        with open(out_path, "w") as out:
            for name, text, _ in lists:
                out.write(f"{name}={text};\n")
        systems, radii, realisations = lists[0][2]
        print(f"{out_path}: AINSS, EINSS, TBsS, each {systems} systems x {radii} radii x {realisations} realisations")


if __name__ == "__main__":
    main()
