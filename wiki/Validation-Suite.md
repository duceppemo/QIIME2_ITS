# Validation suite

`validation/` holds small, **real** (not synthetic) fungal ITS datasets and a script,
`validation/run_validation.sh`, that runs `qiime2-its-train-ncbi`, `qiime2-its-train-fasta`, and
`qiime2-its` against them end-to-end. This is deliberately separate from the unit tests: mocked-
`subprocess` tests verify that this codebase builds the command it intends to, not that the
installed QIIME2 version still accepts that command, or that a real classifier produces sane
taxonomy. Several real QIIME2-version compatibility breaks and real bugs in this codebase were only
found by actually running the pipeline against real data — see `validation/results/*/SUMMARY.md`
for the record.

Scenarios covered by the default run:
- training a classifier from an NCBI accession list, and from a fasta + accession/taxid table
  (then actually classifying real reads with the latter)
- paired-end and single-end reads
- `--min-len`/`--max-len` read-length filtering
- `-rc`/`--reverse_complement`
- multi-sample diversity/composition/classifier stats and the PDF report, including a deliberately
  near-empty sample

## Running it

Inside an activated QIIME2 conda environment with this package installed (`pip install -e .`) and
`bbduk.sh` on `PATH` (see [Installation](Installation)):
```bash
validation/run_validation.sh
```
Exits non-zero on the first failure -- including failed *output* checks, not just crashes: one
taxonomy line per reference sequence for both trainers, a `report.pdf`/`run_metadata.json` and
every ASV classified to a fungal genus for each scenario, and group-significance results for the
multi-sample one. Output and logs go to `validation/output/` (gitignored, regenerated each run).
Needs network access (NCBI Entrez for the 20 reference sequences, and NCBI's ~70 MB
`taxdump.tar.gz`). A custom `output_dir` is wiped at the start of a run, so the script refuses an
existing non-empty directory it didn't create itself.

### `--with-unite` (heavy, opt-in)

```bash
validation/run_validation.sh --with-unite
```
Downloads a real, current UNITE release via `rescript`, trains a full classifier, and classifies
real reads with it. Not part of the default run: the naive-Bayes fit on the full release can take
an hour or more of CPU time and several GB of RAM.

## The audit trail

`validation/results/<date>_<qiime2-version>/` holds dated, committed logs plus a `SUMMARY.md` —
the permanent record of what was validated, against which QIIME2 version, and what was found. When
a bug is found and fixed, the record is corrected transparently (an erratum note, a new dated entry
cross-referencing the old one) rather than silently rewritten. See `validation/data/SOURCES.md` for
where the bundled data comes from and its licensing.
