# Validation run: 2026-09-15

- **QIIME2**: `rachis-qiime2-2026.7` (q2cli 2026.7.0)
- **ITSxpress**: 2.2.0 (via `qiime itsxpress` plugin)
- **Package version under test**: commit `ea349a7` (the compatibility-fix commit this run's
  findings fed into) through the commit that added this validation record
- **Result**: all 4 scenarios in `run_validation.sh` completed successfully (`DONE!`, exit 0)

## What this run covers

Produced by `validation/run_validation.sh` against the bundled real data in `validation/data/`
(see `data/SOURCES.md`). Logs for each scenario are alongside this file.

| scenario | reads in | ASVs | non-chimeric retention |
|---|---|---|---|
| `paired_end` (`-pe --extract-its2`) | 116 / 119 (A/B) | 13 | 59.5% / 56.3% |
| `single_end` (`-se --extract-its2 --max-ee 4 --allow-one-off`) | 113 / 113 | 19 | 74.3% / 79.7% |
| `paired_end_size_filter` (`-pe --extract-its2 --min-len 120 --max-len 155`) | 57 / 66 (post-filter) | 9 | 50.9% / 72.7% |
| `single_end_reverse_complement` (`-se -rc --max-ee 4`) | 113 / 114 | 19 | 75.2% / 81.6% |

All numbers are plausible for real MiSeq ITS3 amplicon data (no near-zero or 100%-loss retention,
no crashes, no empty feature tables). Classifications came back `unidentified` in every scenario,
as expected: the bundled classifier is deliberately narrow (4 unrelated fungal taxa; see
`data/SOURCES.md`) and isn't meant to match these reads' actual taxonomy -- this validates the
classification *mechanism*, not biological accuracy.

Independently spot-checked (not just "didn't crash"):
- `paired_end_size_filter`: confirmed `--min-len`/`--max-len` actually filtered by length (before:
  R1 120-125bp / R2 142-191bp; after: R1 120-124bp / R2 142-155bp, the ~191bp R2 cluster dropped)
  and that R1/R2 stayed paired in sync post-filter (matching read IDs, matching counts).
- `single_end_reverse_complement`: confirmed byte-for-byte that `rc_reads/` holds the exact
  reverse complement of the input sequence and the exact reverse of the quality string, header
  unchanged.

## Bugs found and fixed by earlier runs against this same data (see commit `ea349a7`)

- `qiime tools export` writes a MANIFEST/metadata.yml that broke Casava-format re-import.
- `dada2 denoise-single`/`denoise-paired` require a new `--o-base-transition-stats` output.
- `feature-table summarize` renamed `--m-sample-metadata-file` and split its output into three.
- `feature-classifier classify-sklearn`'s `--p-n-jobs` no longer accepts `-1` for "all".
- `bbduk.sh` (needed for `--min-len`/`--max-len`) isn't installed by the QIIME2 env file and was
  undocumented; added an upfront `env_checks.check_executable()` check plus a README note.

This run is the confirmation that those fixes hold end-to-end, including for the two scenarios
(paired-end + size-filter, and reverse-complement) not covered by the runs that originally found
the bugs above.
