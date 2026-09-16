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

This is a small, deliberately narrow classifier (4 taxa) meant to validate that the training and
classification *code paths* work correctly end-to-end -- not to produce biologically meaningful
taxonomy for the validation reads above (which are not related to these 4 taxa, so classification
against this classifier is expected to come back "unidentified").
