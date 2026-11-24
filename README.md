# dcasm

Microbial genome assemblies from NCBI, each paired with the public SRA/ENA reads it was assembled from,
for which [breseq](https://barricklab.org/breseq) predicts **no residual mutations** when those reads
are mapped back to the assembly.

A pair like this is a known starting point for testing mutation prediction on real reads. Apply a set
of mutations to the assembly (`gdtools APPLY`) to make a mutated reference, and run breseq with the
original reads against it. The mutations breseq should predict are exactly the ones that turn the
mutated reference back into the real assembly. Because the reads gave no predictions against the real
assembly, anything else breseq reports is an error, and anything it misses is a false negative.

Sister repository: [dcamp](https://github.com/barricklab/dcamp), curated mutation predictions for
published evolution experiments.

## Layout

* `assemblies.csv`: one row per assembly. Its columns are described below.
* `references/`: the committed references. `<accession>.isescan.gbk.gz` is the GenBank assembly with
  ISEScan's IS element annotations added, and `<accession>.isescan.csv` is ISEScan's table of what it
  found. Other projects can use these files by their raw URLs.
* `bin/`: the tools (see Tools).
* `build/`: gitignored working directories, one per assembly, where brefito downloads the assembly and
  reads and runs ISEScan and breseq.

There are no per-dataset folders. The committed reference and its `assemblies.csv` row together hold
everything about a pair, and brefito sample sheets are written into `build/<accession>/` when needed.

### assemblies.csv

| Columns | Contents |
| --- | --- |
| `accession` | GenBank assembly accession (`GCA_…`), one row per assembly |
| `organism`, `strain` | from the assembly record |
| `biosample`, `bioproject` | the BioSample the assembly was built from and its BioProject |
| `runs` | the SRA runs (`SRR`/`ERR`/`DRR`, space-separated) used to build the assembly, all from that BioSample |
| `read_types` | brefito read types of those runs: `illumina-PE`, `illumina-SE`, `nanopore` |
| `sra_source` | `reads`, or `mapped-bam` if the SRA loaded the runs from mapped BAM files (see Conventions) |
| `sequencing_tech`, `assembler` | as given in the NCBI assembly report |
| `file` … `command` | how the committed reference was built, written by `bin/make_reference.py`: file name, sequence count and ids, genome size, SHA-256 of the GenBank download, ISEScan / FragGeneScan / breseq versions, `mobile_element` and `repeat_region` feature counts, date and command |
| `n_residual` | number of mutations breseq predicted with these reads against the committed reference |
| `status` | `pending` (not evaluated yet), `clean` (`n_residual` is 0), `residual` (breseq predicted mutations), or `rejected` (unusable for another reason, given in `notes`) |
| `notes` | free text |

## Conventions

* References are GenBank assemblies (GCA), not RefSeq (GCF). An assembly brings all of its chromosomes
  and plasmids, and its sequence ids are the GenBank names that GenomeDiff seq_ids must match.
* ISEScan runs on **every** assembly, even if its GenBank record already annotates IS elements, so all
  references are annotated the same way. Without IS annotations, breseq cannot predict a mobile element
  insertion, and an IS insertion appears only as an unexplained junction. ISEScan's output depends on
  which ISEScan and FragGeneScan are installed, so each reference is built once, committed, and recorded
  with those versions instead of being regenerated for each analysis. The committed file says in each
  record's `COMMENT` that its annotations differ from GenBank's. Its sequence and sequence ids are
  GenBank's.
* Reads are public SRA/ENA runs only, Illumina or Nanopore, and must be the runs used to build the
  assembly: runs of the assembly's own BioSample.
* Mapped-BAM uploads can't be `clean`. Some submitters deposit the BAM file from their own pipeline in
  place of the raw reads. Reads that pipeline didn't keep are missing (often those inside repeats), so
  the absence of predictions proves nothing.
* `clean` requires `n_residual` = 0, `sra_source` = `reads`, a committed reference, and listed runs.

## Adding an assembly

1. Build its reference. This downloads the assembly, runs ISEScan, writes
   `references/<accession>.isescan.gbk.gz` and its table, and adds a `pending` row to `assemblies.csv`:

   ```bash
   bin/make_reference.py GCA_001683415.1
   ```

2. Fill in the row's `biosample`, `bioproject`, `runs`, `read_types`, `sra_source`, `sequencing_tech`
   and `assembler`.
3. Run breseq with the reads against the committed reference, and record `n_residual` and `status`.
   (A driver for this step is still to come.)
4. Check everything:

   ```bash
   python3 bin/check_repo.py
   ```

## Tools

* `bin/make_reference.py <accession>`: builds a committed reference through
  `brefito annotate-genomes --config SKIP_PROKKA=1 --config SKIP_PHISPY=1`. This adds ISEScan's IS
  elements and leaves the existing gene annotations alone. It also writes the build columns of the
  assembly's `assemblies.csv` row and keeps the rest of the row. `--force` rebuilds a reference that
  already exists, `--list` shows what is in `references/`, and `--keep-build` keeps `build/<accession>/`
  for inspection.
* `bin/check_repo.py`: checks the conventions above and lists rows still `pending`.
