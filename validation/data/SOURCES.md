# Data provenance

## `paired_end/`, `single_end/`

Real MiSeq ITS3-primed fungal amplicon reads (sample `4774-1-MSITS3`), taken from
[ITSxpress](https://github.com/USDA-ARS-GBRU/itsxpress)'s own test suite
(`tests/test_data/4774-1-MSITS3_R1.fastq.gz` / `_R2.fastq.gz` / `_merged.fastq`), retrieved
2026-09-15. ITSxpress's `LICENSE.txt` is CC0 1.0 Universal (public domain dedication), so
redistributing this test data here is unrestricted.

- `paired_end/`: the original R1/R2 pair (250 read pairs), split by read index into two
  Casava-named pseudo-samples (`sampleA`, `sampleB`) so the pipeline has >1 sample to work with.
  Purely a splitting of real reads, not synthetic data.
- `single_end/`: the same underlying reads, pre-merged into single full-length ITS reads (227
  reads) by ITSxpress's test fixtures, then split the same way. Used as a real stand-in for
  single-end/full-amplicon platforms (e.g. IonTorrent) that this pipeline also supports, since no
  small public IonTorrent ITS dataset was readily available.
- `multi_sample/` (paired with `metadata_multi.tsv`): the same 250 real read pairs re-split into 6
  samples across 3 nominal "sites" -- 5 healthy (38-60 pairs each) and one deliberately near-empty
  (`siteC-rep2`, 2 pairs, which DADA2 reduces to a zero-read row). Exercises multi-sample diversity
  statistics and the near-empty-sample edge case with real data, rather than the 2-sample
  `metadata.tsv` used elsewhere, which is too small for group-significance tests or a classifier to
  do anything meaningful.

## `classifier_training/`

`accessions.list`: 20 real NCBI RefSeq fungal ITS accessions (a bounded subset of the 66 results
for `txid4751[Organism:exp] AND internal transcribed spacer[Title] AND srcdb_refseq[PROP]`,
spanning *Saccharomyces cerevisiae*, *Eremothecium gossypii*, *Thermochaetoides thermophila*, and
*Sordaria macrospora*), retrieved via NCBI Entrez 2026-09-15. NCBI RefSeq records are U.S.
government work product / public domain.

`acc2taxid.tsv`, `dead_acc2taxid.tsv`: a small, locally-built accession-to-taxid mapping covering
only the above 20 accessions, in the same format as NCBI's `nucl_gb.accession2taxid.gz` /
`dead_nucl.accession2taxid.gz`. These stand in for those files (each several GB) so
`run_validation.sh` can pass them via `--acc2taxid`/`--dead-acc2taxid` and train a real classifier
in under two minutes without a multi-GB download. `dead_acc2taxid.tsv` is just a header row: none
of the chosen accessions are dead/merged.

`seqs.fasta`, `id_table.tsv`: the same 20 accessions/taxids, exported to the fasta + 2-column
accession/taxid format `qiime2-its-train-fasta` expects, for validating that script's own path
(a user-supplied fasta + id table, no NCBI download at all).

This is a small, deliberately narrow classifier (4 taxa) meant to validate that the training and
classification *code paths* work correctly end-to-end, not to produce biologically comprehensive
taxonomy. Classifying the bundled `paired_end`/`single_end` reads against it (which are unrelated
to these 4 taxa) correctly truncates at a shallow, high-confidence rank -- typically
`k__Fungi;p__Ascomycota` -- rather than an empty/"unidentified" result or a wrong species guess:
that's `classify-sklearn`'s confidence-based rank truncation working as intended when a read has no
close match in the reference set, and is the expected outcome here.

**Erratum (2026-09-16)**: `acc2taxid.tsv` originally had corrupted taxid values (Biopython's
`Entrez.esummary` result stringified as `IntegerElement(559292, attributes={})` instead of `559292`
in the one-off script used to build this file) that made every scenario using this classifier
produce all-"unidentified" taxonomy -- a symptom that was misdiagnosed as the expected outcome of a
narrow classifier in `results/2026-09-15_qiime2-2026.7/SUMMARY.md`. The file has since been
corrected; see `results/2026-09-16_qiime2-2026.7/SUMMARY.md` for the re-verified results. This also
surfaced a real, separate bug in `qiime2_its.taxonomy`'s NCBI-rank-to-QIIME2-code mapping, fixed in
the same pass (see that SUMMARY.md for details).
