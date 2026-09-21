# Validation run: 2026-09-16 (corrects 2026-09-15)

> **Erratum (2026-09-21).** The two classifiers trained in this run (and in 2026-09-15's) had
> seen only 4 of their 20 reference sequences: a bug in `taxonomy.parse_id_table()` kept one
> sequence per taxid, and the suite then only checked exit codes. The lineage/rank-mapping fixes
> described below are real and stand, but "taxonomy now correct" did not hold for what those
> classifiers could resolve -- no ASV reached species level. Fixed in 0.3.1; re-validated, with
> output checks added to the suite, in [`../2026-09-21_qiime2-2026.7/SUMMARY.md`](../2026-09-21_qiime2-2026.7/SUMMARY.md).
> This record is otherwise left as written.

- **QIIME2**: `rachis-qiime2-2026.7` (q2cli 2026.7.0)
- **ITSxpress**: 2.2.0 (via `qiime itsxpress` plugin)
- **Result**: all 8 steps (2 classifier trainings + 6 pipeline runs) in `run_validation.sh`
  completed successfully (`DONE!`/`Done`, exit 0). Logs in this directory reflect the final state
  of both validation passes below (the taxid/rank-mapping fixes and the advanced-stats/report
  addition), not two separate sets.

## Why this run exists

While validating `qiime2-its-train-fasta` (not covered by the 2026-09-15 run), two real bugs came
to light that also affected every scenario in that earlier run:

1. **Corrupted `acc2taxid.tsv` fixture**: the one-off script used to build
   `data/classifier_training/acc2taxid.tsv` stringified Biopython's `Entrez.esummary` taxid field
   as `IntegerElement(559292, attributes={})` instead of `559292`. Every taxid in the training data
   was therefore invalid, so `qiime2-its-train-ncbi` silently produced a classifier where every
   reference sequence's lineage was `k__unidentified;p__unidentified;...` -- and every scenario
   using that classifier in the 2026-09-15 run correctly reported "unidentified" for a reason that
   had nothing to do with the (true, but not the actual cause) narrow-classifier explanation
   written down at the time. Fixed by extracting the real taxid from the corrupted string; verified
   by re-running `qiime2-its-train-ncbi` and inspecting `taxonomy.txt` directly (real, complete
   lineages down to species) before re-running the full suite.
2. **`qiime2_its.taxonomy`'s NCBI-rank-to-QIIME2-code mapping was wrong**, independent of bug #1:
   it mapped NCBI rank `superkingdom` to the Kingdom slot and `clade` to the Class slot. Real NCBI
   taxonomy (confirmed by walking the actual downloaded `nodes.dmp` for
   *Saccharomyces cerevisiae*, taxid 559292) uses `kingdom` (-> Fungi) and `class` (-> Saccharomycetes)
   for those; `clade` covers non-Linnean nodes like Opisthokonta, which was being mislabeled as a
   Class. Fixed by correcting `_RANK_TO_CODE` in `src/qiime2_its/taxonomy.py`; two regression tests
   added in `tests/test_taxonomy.py` using a fixture shaped like the real lineage.

Both were data/logic bugs in this project, not QIIME2 compatibility breaks. This run re-verifies
every scenario from 2026-09-15 plus the new `train_fasta`/`fasta_classifier_classification` ones,
against the corrected fixture and code.

## Results

| scenario | reads in | ASVs | non-chimeric retention |
|---|---|---|---|
| `paired_end` (`-pe --extract-its2`) | 116 / 119 (A/B) | 13 | 59.5% / 56.3% |
| `single_end` (`-se --extract-its2 --max-ee 4 --allow-one-off`) | 113 / 113 | 19 | 74.3% / 79.7% |
| `paired_end_size_filter` (`-pe --extract-its2 --min-len 120 --max-len 155`) | 57 / 66 (post-filter) | 9 | 50.9% / 72.7% |
| `single_end_reverse_complement` (`-se -rc --max-ee 4`) | 113 / 114 | 19 | 75.2% / 81.6% |
| `fasta_classifier_classification` (single_end reads, classified with the train_fasta classifier) | 113 / 113 | 19 | 74.3% / 79.7% |

DADA2/ASV numbers are identical to 2026-09-15 for the four original scenarios -- confirming the fix
only changed the classifier's taxonomy output, nothing upstream (import, ITSxpress, DADA2,
phylogeny, diversity).

**Taxonomy now correct**: e.g. `paired_end`'s `biom_table/taxonomy.tsv` assigns
`k__Fungi;p__Ascomycota` (confidence ~0.99-1.0) to most ASVs, several going deeper to
`c__Sordariomycetes;o__Sordariales` or `c__Saccharomycetes;o__Saccharomycetales;f__Saccharomycetaceae`
(confidence ~0.70-0.79). This is the expected, correct behavior of `classify-sklearn`'s
confidence-based rank truncation: the bundled classifier only knows 4 taxa unrelated to these real
reads, so it correctly stops at a shallow rank it's actually confident about rather than forcing a
specific wrong species, or (the 2026-09-15 bug's symptom) returning nothing at all.

`train_fasta`'s own `taxonomy.txt` (built directly from `data/classifier_training/seqs.fasta` +
`id_table.tsv`, no NCBI download) shows complete, correct 7-rank lineages for all 4 source taxa,
e.g. `NR_132210.1  k__Fungi;p__Ascomycota;c__Saccharomycetes;o__Saccharomycetales;f__Saccharomycetaceae;g__Saccharomyces;s__Saccharomyces_cerevisiae`.

## Bugs found in earlier runs, still confirmed fixed (see commit `ea349a7` and `8acc575`)

- `qiime tools export` MANIFEST/metadata.yml breaking Casava-format re-import.
- `dada2 denoise-single`/`denoise-paired` requiring `--o-base-transition-stats`.
- `feature-table summarize`'s renamed/restructured flags.
- `classify-sklearn`'s `--p-n-jobs` semantics (`-1` -> `0`).
- `bbduk.sh` availability check.
- `qiime2-its-train-unite`'s ambiguous file matching, leaked `.fasta` extension, dead `-t` flag.

None of these regressed in this run.

## Second pass this same day: advanced stats, sample classification, and the PDF report

Added `qiime diversity alpha-group-significance`/`beta-group-significance`, `qiime taxa collapse` +
`feature-table relative-frequency`, `qiime sample-classifier classify-samples`, and a PDF report
(`report.pdf`) built from the pipeline's own output, run by default after the core pipeline. Also
added a new bundled dataset, `data/multi_sample/` + `data/metadata_multi.tsv` (6 real samples, 5
healthy + 1 deliberately near-empty, 4 metadata columns), as `run_validation.sh`'s
`multi_sample_advanced_stats` scenario -- the 2-sample `metadata.tsv` used by the other scenarios
has only one column with one unique value per sample, too small for any of this to exercise
meaningfully.

Real bugs found and fixed while validating this against real data (all in `src/qiime2_its/`, not
in `validation/`):

1. **`qiime diversity alpha-group-significance` fails outright, for the whole metadata file, if no
   categorical column has both a repeated and a varying value** (its own message: "doesn't consist
   of unique values, and doesn't consist of exactly one value") -- exactly the case for the
   2-sample/1-column `metadata.tsv`, which crashed `paired_end` and every other small scenario the
   moment this feature was added. Added `metadata_utils.has_alpha_group_significance_column()`
   (a different, looser eligibility rule than the >=2-per-group one already used for
   beta-group-significance/classify-samples) and gated the alpha-group-significance calls on it,
   printing a skip message instead.
2. **`qiime taxa collapse --p-level 6` (genus) fails outright if the classifier never resolved
   genus for any feature in the dataset** -- routine for a classifier applied to reads unrelated to
   its training set (confirmed with the bundled toy NCBI classifier, which topped out at level 5 /
   family for `multi_sample_advanced_stats`'s real reads). Added
   `taxonomy.max_lineage_depth()` and capped the requested collapse level to
   `min(6, max_lineage_depth(...))` rather than assuming genus is always reachable.
3. **A near-empty sample DADA2 reduces to a zero-read row can break `classify-samples`'s internal
   stratified split even for a metadata column that looks balanced** (e.g. 2 samples per site on
   paper): confirmed directly by running `qiime sample-classifier classify-samples -m site` by hand
   against the real 6-sample table and getting "one or more values that match only one sample" once
   the zero-read sample effectively orphaned one group. Fixed by reading QIIME2's own
   `sample-frequencies.qza` export to determine the *actual* non-zero sample set before computing
   eligibility for beta-group-significance/classify-samples, and by wrapping `classify_samples` in
   a per-column try/except so one bad column doesn't abort the run (both already partly designed
   for this; this run is what proved the design necessary and correct).

Verified end to end with `multi_sample_advanced_stats`: 4 alpha metrics tested across 4 categorical
columns, 4 columns x 2 distance metrics for beta-group-significance (8 PERMANOVA results), genus
collapse correctly fell back to level 5, the `replicate` column's classifier trained successfully
(overall accuracy 0.33 on a 2-class, 6-sample problem -- unremarkable, but the mechanism works) while
`site`/`host-plant`/`collection-date` were correctly skipped with a clear reason, and `report.pdf`
(10 pages) rendered every section with real data -- individually confirmed by rendering each page to
PNG and inspecting it, not just checking the file exists.

## Third pass this same day: QA/provenance metadata

Added `run_metadata.json` (written unconditionally, even with `--skip-report`) and four new report
pages -- Run information (start/end/duration, user@host, platform, conda env, qiime2-its/QIIME2
versions, input/metadata/classifier/output paths, the exact command invoked), Pipeline parameters
(every CLI flag's value), Input sample files (one row per fastq file), and Installed QIIME2 plugins
(every plugin `qiime info` reports, name + version) -- for run-to-run QA/audit purposes.

One real bug found and fixed via this pass's real-data run (not present in a mocked-subprocess unit
test, since `qiime info` was never called before): the beta-group-significance table fix from
earlier today established `add_table_page()` truncates any cell too wide for its column, which is
correct for short values but wrong for this page's file-path/command-line values -- silently
dropping part of a path would defeat the point of a QA record. Added `add_keyvalue_page()` (a
label + wrapping `multi_cell()` value, never truncated) instead of reusing `add_table_page()` for
the Run information page; verified by regenerating the report against a real 8.4-minute
`paired_end` run with a long absolute output path and inspecting the rendered page directly (full
path visible, wrapped onto a second line, nothing cut off).

Re-ran the full `run_validation.sh` suite (all 6 `qiime2-its` scenarios plus `train_ncbi`/
`train_fasta`) to confirm `run_metadata.json` is written and the new pages render across every
scenario, not just the one manually inspected; refreshed this directory's `.log` files to match.
