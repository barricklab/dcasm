#!/usr/bin/env python3
"""Check dcasm repository conventions. Run from the repository root.

Checks:
  * assemblies.csv has the expected columns, and one row per GenBank assembly (GCA_...)
  * biosample is a BioSample accession (SAMN/SAMEA/SAMD), runs are SRA run accessions, read_types
    are ones brefito supports, sra_source is 'reads', 'mapped-bam' or empty
  * a row with a file has references/<accession>.isescan.gbk.gz and its ISEScan table on disk, and
    every file in references/ has a row
  * status is pending, clean, residual or rejected; a clean row has n_residual 0 and reads that were
    not loaded from mapped BAM files; a residual row has n_residual above 0
  * no private lab-server paths appear in assemblies.csv
Rows still pending are listed but are not problems. Exits 1 if any problem is found.
"""

import csv
import glob
import os
import re
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from make_reference import CSV_PATH, FIELDS, REF_DIR, file_name  # noqa: E402

PRIVATE = re.compile(r"BarrickLab-Private|BarrickLab-Public|/scratch/|corral|utexas_gsaf|msu_rtsf", re.I)
ACCESSION_RE = re.compile(r"^GCA_\d+\.\d+$")
BIOSAMPLE_RE = re.compile(r"^SAM(N|EA|D)\d+$")
RUN_RE = re.compile(r"^[SED]RR\d+$")
READ_TYPES = {"illumina-PE", "illumina-SE", "nanopore"}
SRA_SOURCES = {"", "reads", "mapped-bam"}
STATUSES = {"pending", "clean", "residual", "rejected"}


def check_row(r, problems, pending):
    acc = r["accession"]
    where = f"{CSV_PATH}: {acc}"
    if not ACCESSION_RE.match(acc):
        problems.append(f"{where}: accession is not a GenBank assembly (GCA_...)")
    if r["biosample"] and not BIOSAMPLE_RE.match(r["biosample"]):
        problems.append(f"{where}: biosample {r['biosample']} is not a BioSample accession")
    for run in r["runs"].split():
        if not RUN_RE.match(run):
            problems.append(f"{where}: run {run} is not an SRA run accession")
    for t in r["read_types"].split():
        if t not in READ_TYPES:
            problems.append(f"{where}: read type {t} is not one of {sorted(READ_TYPES)}")
    if r["sra_source"] not in SRA_SOURCES:
        problems.append(f"{where}: sra_source '{r['sra_source']}' is not 'reads', 'mapped-bam' or empty")
    if any(PRIVATE.search(v or "") for v in r.values()):
        problems.append(f"{where}: mentions a private path")

    if r["file"]:
        if r["file"] != file_name(acc):
            problems.append(f"{where}: file is {r['file']}, expected {file_name(acc)}")
        elif not os.path.exists(os.path.join(REF_DIR, r["file"])):
            problems.append(f"{where}: {REF_DIR}/{r['file']} not found")
        if not os.path.exists(os.path.join(REF_DIR, acc + ".isescan.csv")):
            problems.append(f"{where}: {REF_DIR}/{acc}.isescan.csv (the ISEScan table) not found")

    status, n = r["status"], r["n_residual"]
    if status not in STATUSES:
        problems.append(f"{where}: status '{status}' is not one of {sorted(STATUSES)}")
    if n and not n.isdigit():
        problems.append(f"{where}: n_residual '{n}' is not a count")
    if status == "pending":
        pending.append(f"{where}: " + ("no reference built yet" if not r["file"] else "breseq not run yet"))
    elif status == "clean":
        if n != "0":
            problems.append(f"{where}: clean but n_residual is '{n}'")
        if r["sra_source"] != "reads":
            problems.append(f"{where}: clean but sra_source is '{r['sra_source']}', not 'reads'")
        if not r["file"]:
            problems.append(f"{where}: clean but no reference has been built")
        if not r["runs"]:
            problems.append(f"{where}: clean but no runs are listed")
    elif status == "residual" and not (n.isdigit() and int(n) > 0):
        problems.append(f"{where}: residual but n_residual is '{n}'")


def main():
    problems, pending = [], []
    rows = []
    if not os.path.exists(CSV_PATH):
        problems.append(f"{CSV_PATH}: missing")
    else:
        reader = csv.DictReader(open(CSV_PATH, newline=""))
        if reader.fieldnames != FIELDS:
            problems.append(f"{CSV_PATH}: header must be {','.join(FIELDS)}")
        else:
            rows = list(reader)

    seen = set()
    for r in rows:
        if r["accession"] in seen:
            problems.append(f"{CSV_PATH}: {r['accession']} has more than one row")
        seen.add(r["accession"])
        check_row(r, problems, pending)

    listed = {r["file"] for r in rows if r["file"]}
    for path in sorted(glob.glob(os.path.join(REF_DIR, "*"))):
        name = os.path.basename(path)
        if name.startswith("."):
            continue
        if name.endswith(".gbk.gz"):
            if name not in listed:
                problems.append(f"{path}: no {CSV_PATH} row names this file")
        elif name.endswith(".isescan.csv"):
            if name[:-len(".isescan.csv")] not in seen:
                problems.append(f"{path}: no {CSV_PATH} row for this assembly")
        else:
            problems.append(f"{path}: not a reference or an ISEScan table")

    for p in pending:
        print("pending: " + p)
    for p in problems:
        print(p)
    counts = {s: sum(r["status"] == s for r in rows) for s in sorted(STATUSES)}
    print(f"\n{len(rows)} assemblies ({', '.join(f'{v} {k}' for k, v in counts.items())}); "
          f"{len(problems)} problem(s)")
    sys.exit(1 if problems else 0)


if __name__ == "__main__":
    main()
