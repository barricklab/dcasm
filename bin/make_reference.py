#!/usr/bin/env python3
"""Build a committed reference for dcasm: a GenBank assembly with its IS elements annotated by ISEScan.

breseq needs IS element annotations to call a mobile element insertion: without them an IS insertion
shows up only as an unexplained junction (JC) rather than a MOB mutation. ISEScan finds the elements
and `breseq CONVERT-REFERENCE -s` writes them into the reference.

ISEScan runs on every assembly, including ones whose GenBank record already annotates IS elements, so
that every reference in the repository is annotated the same way. Its output depends on which ISEScan
and FragGeneScan are installed, so the annotated reference is built once, committed to references/,
and its versions and the command that produced it are recorded in assemblies.csv.

  bin/make_reference.py GCA_001683415.1          # build it
  bin/make_reference.py GCA_001683415.1 --force  # rebuild one that already exists
  bin/make_reference.py --list                   # what is in references/ already

The result is references/<accession>.isescan.gbk.gz, which says in the COMMENT of each record that its
annotations differ from GenBank's, with ISEScan's table beside it as <accession>.isescan.csv. The
sequence and sequence ids are GenBank's, so coordinates and seq_ids match the assembly.

Only the build columns of the assembly's assemblies.csv row are written; the rest of an existing row
(runs, status, notes) is kept. A new row starts with status=pending.
"""

import argparse
import csv
import datetime
import glob
import gzip
import hashlib
import os
import re
import shutil
import subprocess
import sys

BREFITO_BIN = os.path.expanduser("~/miniconda3/envs/brefito/bin")
REF_DIR = "references"
BUILD_DIR = "build"
CSV_PATH = "assemblies.csv"
# Every column of assemblies.csv, in order. check_repo.py reads the same list.
PAIR_FIELDS = ["accession", "organism", "strain", "biosample", "bioproject", "runs", "read_types",
               "sra_source", "sequencing_tech", "assembler"]
BUILD_FIELDS = ["file", "n_sequences", "seq_ids", "genome_size", "source_sha256", "isescan_version",
                "fraggenescan_version", "breseq_version", "n_mobile_element", "n_repeat_region",
                "date", "command"]
RESULT_FIELDS = ["n_residual", "status", "notes"]
FIELDS = PAIR_FIELDS + BUILD_FIELDS + RESULT_FIELDS


def run(cmd, cwd=None, env=None, check=True, log=None):
    """Run a command, keeping its output in `log` so the conda envs it used can be read back."""
    print("  $ " + " ".join(cmd), flush=True)
    with open(log, "w") if log else open(os.devnull, "w") as f:
        r = subprocess.run(cmd, cwd=cwd, env=env,
                           stdout=f if log else None,
                           stderr=subprocess.STDOUT if log else None)
    if r.returncode and log:
        print("".join(open(log, errors="replace").readlines()[-25:]))
    if check and r.returncode:
        sys.exit(f"failed: {' '.join(cmd)}")
    return r.returncode


def brefito_env():
    return dict(os.environ, PATH=BREFITO_BIN + os.pathsep + os.environ["PATH"])


def count_features(text):
    """mobile_element and repeat_region feature counts in a GenBank flat file."""
    return (len(re.findall(r"^     mobile_element ", text, re.M)),
            len(re.findall(r"^     repeat_region ", text, re.M)))


def describe(text):
    """(organism, strain, seq_ids, total length) read from a GenBank flat file."""
    loci = re.findall(r"^LOCUS\s+(\S+)\s+(\d+)\s+bp", text, re.M)
    seq_ids = [l[0] for l in loci] or re.findall(r"^LOCUS\s+(\S+)", text, re.M)
    size = sum(int(l[1]) for l in loci)
    organism = ""
    m = re.search(r"^\s+ORGANISM\s+(.+)$", text, re.M)
    if m:
        organism = m.group(1).strip()
    strain = ""
    m = re.search(r'/strain="([^"]+)"', text)
    if m:
        strain = m.group(1).strip()
    return organism, strain, seq_ids, size


def file_name(acc):
    """Name of the committed file: .isescan marks that its annotations differ from GenBank's."""
    return acc + ".isescan.gbk.gz"


def mark_modified(text, acc, isescan_version):
    """Say in each record's COMMENT that the annotations are not the ones GenBank distributes."""
    note = (f"Modified from GenBank assembly {acc}: IS elements annotated with ISEScan"
            f"{' ' + isescan_version if isescan_version else ''} by dcasm bin/make_reference.py "
            f"(mobile_element features with /source=\"isescan\"). The sequence is unchanged.")
    lines, line = [], ""
    for word in note.split():
        if line and len(line) + 1 + len(word) > 67:   # 12 columns of indent, 79 in all
            lines.append(line)
            line = word
        else:
            line = (line + " " + word).strip()
    lines.append(line)
    block = "COMMENT     " + ("\n" + " " * 12).join(lines) + "\n"

    def mark(record):
        if "Modified from GenBank assembly" in record.split("\nFEATURES", 1)[0]:
            return record
        if re.search(r"^COMMENT     ", record, re.M):
            # keep the submitters' comment, after ours and a blank comment line
            return re.sub(r"^COMMENT     ", lambda m: block + " " * 12 + "\n" + " " * 12,
                          record, count=1, flags=re.M)
        return re.sub(r"^FEATURES ", lambda m: block + "FEATURES ", record, count=1, flags=re.M)

    return "".join(mark(r) for r in re.split(r"(?m)^(?=LOCUS )", text) if r)


def read_text(path):
    opener = gzip.open if path.endswith(".gz") else open
    with opener(path, "rt", errors="replace") as f:
        return f.read()


def env_versions(env_dir, names=("isescan", "fraggenescan")):
    """Package versions actually installed in a conda env snakemake built."""
    out = {}
    for name in names:
        for meta in glob.glob(os.path.join(env_dir, "conda-meta", name + "-*.json")):
            m = re.match(rf"{re.escape(name)}-([^-]+)-", os.path.basename(meta))
            if m:
                out[name] = m.group(1)
    return tuple(out.get(n, "") for n in names)


def envs_used(build_dir):
    """The conda env directories this build's snakemake jobs actually activated.

    Read from the run logs rather than by globbing every env on the machine: several envs can hold
    a different isescan or breseq, and recording a version the run did not use would make the
    provenance wrong in exactly the way it exists to prevent."""
    dirs = []
    for log in sorted(glob.glob(os.path.join(build_dir, "*.log"))):
        for m in re.finditer(r"Activating conda environment: (\S+)", open(log, errors="replace").read()):
            d = os.path.normpath(os.path.join(build_dir, m.group(1)))
            if os.path.isdir(d) and d not in dirs:
                dirs.append(d)
    return dirs


def tool_versions(build_dir):
    """(isescan, fraggenescan, breseq) versions from the envs this build used."""
    isescan = fraggenescan = breseq = ""
    for d in envs_used(build_dir):
        i, f = env_versions(d)
        isescan, fraggenescan = isescan or i, fraggenescan or f
        if not breseq:
            exe = os.path.join(d, "bin", "breseq")
            if os.path.exists(exe):
                r = subprocess.run([exe, "--version"], capture_output=True, text=True)
                m = re.search(r"breseq\s+(\S+(?:\s+revision\s+\S+)?)", r.stdout + r.stderr)
                breseq = m.group(1) if m else env_versions(d, ("breseq-prerelease", "breseq"))[0]
    return isescan, fraggenescan, breseq


def csv_rows():
    if not os.path.exists(CSV_PATH):
        return []
    return list(csv.DictReader(open(CSV_PATH, newline="")))


def write_csv(rows):
    rows.sort(key=lambda r: r["accession"])
    with open(CSV_PATH, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=FIELDS, lineterminator="\n")
        w.writeheader()
        w.writerows({k: r.get(k, "") for k in FIELDS} for r in rows)


def cmd_list(args):
    rows = {r["file"]: r for r in csv_rows() if r.get("file")}
    files = sorted(glob.glob(os.path.join(REF_DIR, "*.gbk.gz")))
    print(f"{len(files)} reference file(s) in {REF_DIR}/, {len(rows)} built row(s) in {CSV_PATH}")
    for path in files:
        name = os.path.basename(path)
        r = rows.get(name, {})
        print(f"  {name:<32} {os.path.getsize(path) / 1e6:5.1f} MB  "
              f"IS={r.get('n_mobile_element', '?'):>4}  status={r.get('status', '?'):<8} "
              f"{r.get('organism', '(no row in ' + CSV_PATH + ')')[:40]}")
    for name in sorted(set(rows) - {os.path.basename(p) for p in files}):
        print(f"  {name:<32} in {CSV_PATH} but the file is missing")


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("accession", nargs="?", help="GenBank assembly accession, e.g. GCA_001683415.1")
    p.add_argument("--list", action="store_true", help="list what is already in references/")
    p.add_argument("--force", action="store_true", help="rebuild even if the file exists")
    p.add_argument("--keep-build", action="store_true", help="keep the build directory for inspection")
    args = p.parse_args()

    if args.list or not args.accession:
        return cmd_list(args)
    acc = args.accession
    if not re.match(r"^GCA_\d+\.\d+$", acc):
        sys.exit(f"{acc}: expected a GenBank assembly accession such as GCA_001683415.1")

    out_path = os.path.join(REF_DIR, file_name(acc))
    if os.path.exists(out_path) and not args.force:
        sys.exit(f"{out_path} exists; pass --force to rebuild it")

    build = os.path.join(BUILD_DIR, acc)
    os.makedirs(build, exist_ok=True)
    with open(os.path.join(build, "data.csv"), "w", newline="") as f:
        w = csv.writer(f, lineterminator="\n")
        w.writerow(["sample", "type", "setting"])
        w.writerow([acc, "reference", "ncbi-genome://" + acc])

    print(f"{acc}: downloading the assembly")
    run([os.path.join(BREFITO_BIN, "brefito"), "download-data", "--breseq-prerelease"],
        cwd=build, env=brefito_env(), log=os.path.join(build, "download-data.log"))
    downloaded = os.path.join(build, "references", acc)
    if not os.path.exists(downloaded):
        sys.exit(f"{downloaded}: brefito did not produce the reference")
    source_text = read_text(downloaded)
    source_sha = hashlib.sha256(open(downloaded, "rb").read()).hexdigest()
    before = count_features(source_text)
    organism, strain, seq_ids, size = describe(source_text)
    print(f"{acc}: {organism} {strain}, {len(seq_ids)} sequence(s) {' '.join(seq_ids)}, {size} bp; "
          f"GenBank has {before[0]} mobile_element, {before[1]} repeat_region")

    print(f"{acc}: running ISEScan and adding its results to the reference")
    # --notemp keeps the ISEScan output directory, which the workflow marks temp(); without it
    # snakemake deletes the table as soon as it has been folded into the reference.
    config = ["--config", "SKIP_PROKKA=1", "--config", "SKIP_PHISPY=1", "--notemp"]
    command = "brefito annotate-genomes --breseq-prerelease " + " ".join(config)
    run([os.path.join(BREFITO_BIN, "brefito"), "annotate-genomes", "--breseq-prerelease"] + config,
        cwd=build, env=brefito_env(), log=os.path.join(build, "annotate-genomes.log"))
    produced = os.path.join(build, "annotated-references", acc + ".gbk")
    if not os.path.exists(produced):
        sys.exit(f"{produced}: annotate-genomes did not produce the annotated reference")
    isescan_version, fraggenescan_version, breseq_version = tool_versions(build)
    final_text = mark_modified(read_text(produced), acc, isescan_version)
    if not isescan_version:
        print("  WARNING: could not tell which ISEScan the run used, so assemblies.csv "
              "will not record its version")
    # the ISEScan table itself, committed beside the reference as the record of what was found
    tables = glob.glob(os.path.join(build, "annotated-references-isescan", acc, "**", "*.csv"),
                       recursive=True)
    if tables:
        shutil.copyfile(tables[0], os.path.join(REF_DIR, acc + ".isescan.csv"))
        print(f"  kept the ISEScan table: {REF_DIR}/{acc}.isescan.csv")
    else:
        print("  WARNING: no ISEScan table was kept, so what it found is recorded only "
              "in the reference's own features")

    after = count_features(final_text)
    os.makedirs(REF_DIR, exist_ok=True)
    with gzip.open(out_path, "wt") as f:
        f.write(final_text)
    print(f"{acc}: wrote {out_path} ({os.path.getsize(out_path) / 1e6:.1f} MB gzipped, "
          f"{len(final_text) / 1e6:.1f} MB raw); "
          f"{after[0]} mobile_element, {after[1]} repeat_region")

    rows = csv_rows()
    row = next((r for r in rows if r["accession"] == acc), None)
    if row is None:
        row = {"accession": acc, "status": "pending"}
        rows.append(row)
    # organism and strain come from the assembly only if nothing has been filled in already
    row["organism"] = row.get("organism") or organism
    row["strain"] = row.get("strain") or strain
    row.update({
        "file": os.path.basename(out_path), "n_sequences": len(seq_ids), "seq_ids": " ".join(seq_ids),
        "genome_size": size, "source_sha256": source_sha, "isescan_version": isescan_version,
        "fraggenescan_version": fraggenescan_version, "breseq_version": breseq_version,
        "n_mobile_element": after[0], "n_repeat_region": after[1],
        "date": datetime.date.today().isoformat(), "command": command,
    })
    write_csv(rows)
    print(f"{CSV_PATH}: recorded {acc}")
    if not args.keep_build:
        shutil.rmtree(build, ignore_errors=True)


if __name__ == "__main__":
    main()
