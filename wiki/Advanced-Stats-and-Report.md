# Advanced stats and the PDF report

By default, after the core pipeline (import → ITSxpress → DADA2 → phylogeny → diversity →
taxonomy → barplot), `qiime2-its` also runs:

1. **Alpha diversity group significance** (`qiime diversity alpha-group-significance`, Kruskal-Wallis)
   for each of Faith's PD, observed features, Shannon, and evenness — but only if the metadata file
   has at least one categorical column with a repeated (not all-unique) and varying (not all-same)
   value; otherwise it's skipped with a message, since the action itself fails outright for the
   whole metadata file rather than degrading per-column.
2. **Beta diversity group significance** (`qiime diversity beta-group-significance`, PERMANOVA)
   against Bray-Curtis and unweighted UniFrac, for every "eligible" categorical metadata column: ≥2
   distinct values, each held by ≥2 samples.
3. **Genus-level composition**: `qiime taxa collapse` (capped to whatever taxonomic depth the
   classifier actually resolved for this dataset — a classifier commonly won't reach genus for
   reads unrelated to its training set, and requesting a deeper level than that fails outright) +
   `qiime feature-table relative-frequency`.
4. **Sample classification** (`qiime sample-classifier classify-samples`, random forest) for each
   eligible column, with cross-validation folds capped to the smallest class size. A column that's
   nominally eligible but still too small for scikit-learn's internal stratified split (e.g. one
   category has only one non-zero-read sample) is caught and skipped with a message, not fatal to
   the run.
5. **A PDF report** (`<output>/report.pdf`): run summary, QA/provenance pages (see below), DADA2
   read retention table, representative-sequence length distribution, alpha diversity
   group-significance table + boxplots, beta diversity group-significance table + PCoA plots and
   UPGMA clustering dendrograms (Bray-Curtis and unweighted UniFrac), genus-level relative-abundance
   chart, taxonomic classification confidence distribution, the rarefaction curve, and
   sample-classifier accuracy (if any column succeeded). Built with matplotlib +
   [fpdf2](https://pypi.org/project/fpdf2/) + [scipy](https://scipy.org/) — no LaTeX, no browser
   rendering, no heavy system dependencies. It discovers what to include by scanning the output
   folder for whichever of the artifacts above were actually produced, so it degrades gracefully if
   you pass `--skip-advanced-stats`.

   The alpha/beta figures aren't limited to one metadata column: every eligible column whose
   group-significance test comes back statistically significant (p < 0.05, for at least one metric)
   gets its own boxplot/PCoA/dendrogram, in addition to `--report-metadata-column`'s column (always
   shown, even when it isn't itself significant, so there's always at least one grouped view). Beta
   diversity figures are grouped one sub-section per metadata column — both distance metrics' PCoA
   and dendrogram together — rather than interleaved, and each page's intro text says why that
   column is shown (significant vs. default).

## QA / provenance

Every run writes `<output>/run_metadata.json` — regardless of `--skip-report` — recording:

- when the run started/finished and how long it took
- who ran it and on which machine (username, hostname, platform)
- the exact command invoked, and every CLI parameter's value
- the input folder, metadata file, classifier file, and output folder used
- every input fastq file, grouped by sample
- the QIIME2 framework version and every installed plugin's version (`qiime info`, captured at
  run time)

If `report.pdf` is built, this becomes four pages near the front — Run information, Pipeline
parameters, Input sample files, and Installed QIIME2 plugins — before the diversity/composition
results. A report built from an older output folder that predates this (or one missing
`run_metadata.json` for any other reason) just skips those pages rather than failing.

## A near-empty sample

A sample that DADA2 reduces to a zero-read row (its input was low enough that nothing survives
denoising/chimera removal) still exists in the metadata file, but is excluded from all of the
eligibility checks above via the sample-frequency export QIIME2 itself produces
(`sample-frequencies.qza`) — otherwise a metadata column that looks balanced on paper (e.g. 2
samples per site) can turn out to have an effectively-singleton group once the zero-read sample is
accounted for, which breaks `classify-samples`'s internal stratification even though the column
passed the simpler ≥2-per-group eligibility check.

## Turning it off

```bash
qiime2-its ... --skip-advanced-stats --skip-report
```
`--skip-advanced-stats` skips steps 1–4 above; `--skip-report` skips just the PDF (useful if you
want the group-significance/classifier artifacts but not the report, or vice versa).

## Output layout

```
<output>/
  alpha-group-significance-<metric>.qzv          # one per alpha metric tested
  beta-group-significance-<column>-<metric>.qzv  # one per (eligible column, distance metric)
  table-genus.qza, table-genus-relative.qza
  sample-classifier-<column>/                    # one per eligible, non-skipped column
  dada2_stats/, sample_frequencies/               # plain-text exports the report reads
  rep_seqs_export/dna-sequences.fasta             # plain-text FASTA the report reads for
                                                   # its sequence-length distribution
  biom_table/taxonomy.tsv                         # also the source of the report's
                                                   # classification-confidence distribution
  core-metrics-results/<metric>_pcoa_export/      # plain-text PCoA exports the report reads
  core-metrics-results/<metric>_distance_export/  # plain-text distance-matrix exports the
                                                   # report reads for its UPGMA dendrograms
  run_metadata.json                              # QA/provenance record, see above
  report.pdf
```
