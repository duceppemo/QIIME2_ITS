# Real-data validation

The unit tests under `tests/` mock every `subprocess.run` call to `qiime`/`biom`/`bbduk.sh`/etc.
They verify that this codebase builds the command it intends to -- they cannot verify that the
QIIME2 version actually installed still accepts that command. QIIME2's CLI has changed interface
(renamed flags, new required outputs, changed parameter semantics) between releases, and those
breaks only show up by actually running the pipeline.

This directory holds small, real (not synthetic) fungal ITS datasets and a script that runs
`qiime2-its` against them end-to-end, covering:

- paired-end reads
- single-end reads (a stand-in for full-amplicon single-end platforms like IonTorrent)
- `--min-len`/`--max-len` read-length filtering (paired-end)
- `-rc`/`--reverse_complement`

See `data/SOURCES.md` for where the bundled data comes from and its licensing.

## Running it

Inside an activated QIIME2 conda environment with this package installed (`pip install -e .`) and
`bbduk.sh` on PATH (`conda install -c bioconda -c conda-forge bbmap` -- see the main README):

```bash
validation/run_validation.sh
```

This trains a small real classifier from bundled NCBI accessions (no multi-GB taxonomy download
needed -- see `data/SOURCES.md`), then runs all four scenarios above, writing full pipeline output
and logs to `validation/output/` (gitignored -- it's regenerated each run, not a permanent record).

The script exits non-zero on the first failure (`set -euo pipefail`), so a clean exit means all
four scenarios completed without error.

## `--with-unite` (heavy, opt-in)

```bash
validation/run_validation.sh --with-unite
```

Additionally downloads a real, current UNITE release via the `rescript` plugin
(`qiime rescript get-unite-data`, bundled with QIIME2 -- see the main README's "From UNITE"
section), repackages it into a real UNITE QIIME-release archive layout, and runs it through
`qiime2-its-train-unite` and then `qiime2-its` itself, so the resulting classifier assigns real
taxonomy to the bundled `single_end` reads.

This is **not** part of the default run and nothing from it is bundled in the repo (unlike the
~100KB of data used by the four scenarios above):
- it downloads tens of MB from UNITE at runtime (UNITE data is CC BY-SA 4.0, see
  https://unite.ut.ee/cite.php -- separate from the CC0/public-domain licensing of the data in
  `data/`)
- the naive-Bayes fit on the full release can take on the order of an **hour or more of CPU time**
  and **several GB of RAM**

Override `UNITE_VERSION`/`UNITE_TAXON_GROUP`/`UNITE_CLUSTER_ID` env vars to use a different release
than the script's default (`2025-02-19`, `fungi`, `99`); see `qiime rescript get-unite-data --help`
for the currently available versions.

## `results/`

Dated, committed records of full validation runs -- the actual logs plus a `SUMMARY.md` noting the
QIIME2 version, git commit, and any bugs found. This is the audit trail: `validation/output/` is
scratch, `results/<date>_<qiime2-version>/` is the permanent record. Add a new dated subdirectory
each time a full validation run is performed against a new QIIME2 version or after a change likely
to affect runtime behavior (not for every commit).
