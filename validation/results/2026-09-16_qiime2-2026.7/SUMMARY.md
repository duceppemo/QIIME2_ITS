# Validation run: 2026-09-16 (corrects 2026-09-15)

- **QIIME2**: `rachis-qiime2-2026.7` (q2cli 2026.7.0)
- **ITSxpress**: 2.2.0 (via `qiime itsxpress` plugin)
- **Result**: all 7 steps (2 classifier trainings + 5 pipeline runs) in `run_validation.sh`
  completed successfully (`DONE!`/`Done`, exit 0)

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
