# Building a classifier

You need a QIIME2 classifier to run this pipeline. If you update your QIIME2 version, you will
generally need to recompile the classifier.

## From UNITE

**Preferred method**: the `rescript` plugin, included in the QIIME2 distribution, downloads UNITE
directly from its PlutoF REST API -- no manual download needed, and it's the way UNITE itself now
recommends fetching this data for QIIME2 (`https://unite.ut.ee/repository.php` now gates its files
behind DOI landing pages with no stable direct-download URL, which is why the URL-based
auto-download this pipeline used before generally no longer works):
```bash
conda activate rachis-qiime2-2026.7

qiime rescript get-unite-data \
    --p-version 2025-02-19 \
    --p-taxon-group fungi \
    --p-cluster-id 99 \
    --p-no-singletons \
    --o-sequences unite-sequences.qza \
    --o-taxonomy unite-taxonomy.qza

qiime feature-classifier fit-classifier-naive-bayes \
    --i-reference-reads unite-sequences.qza \
    --i-reference-taxonomy unite-taxonomy.qza \
    --o-classifier unite-classifier.qza
```
Run `qiime rescript get-unite-data --help` for the current list of available `--p-version` values.

**Fallback**: if you already have a UNITE QIIME-release archive (`.tgz`) from another source (a
manually downloaded copy, a colleague, an institutional mirror), `qiime2-its-train-unite` can
extract, import, and train from it directly:
```bash
qiime2-its-train-unite \
    -u ~/Downloads/sh_qiime_release_19.02.2025.tgz \
    -o /db/UNITE \
    -q rachis-qiime2-2026.7
```
`-u` also accepts a URL, but only if it resolves directly to the archive file -- most current UNITE
download links do not.

## From Eukaryome

[Eukaryome](https://eukaryome.org/) is a broader alternative/complement to UNITE, also fetchable
directly via `rescript` -- but unlike UNITE it isn't fungi-only: it covers ITS across *all*
eukaryotes (animals, plants, protists included), which can actually help ITS metabarcoding by
letting a classifier flag non-fungal amplification instead of forcing everything into a fungal
lineage. Confirmed by a real download (Eukaryome 2.0's `--p-rrna-gene ITS`): ~1.6 million
sequences, ~200 MB -- an order of magnitude more than UNITE, so training takes correspondingly
longer and more RAM.
```bash
qiime rescript get-eukaryome-data \
    --p-rrna-gene ITS \
    --output-dir eukaryome_its

qiime feature-classifier fit-classifier-naive-bayes \
    --i-reference-reads eukaryome_its/eukaryome_sequences/ITS_seqs.qza \
    --i-reference-taxonomy eukaryome_its/eukaryome_taxonomy/ITS_taxa.qza \
    --o-classifier eukaryome-classifier.qza
```
`--p-rrna-gene` is a Collection output, hence `--output-dir` rather than individual `--o-...`
flags -- it writes `<dir>/eukaryome_sequences/ITS_seqs.qza` and
`<dir>/eukaryome_taxonomy/ITS_taxa.qza` (naming matches whatever gene(s) you requested). If you
only want fungi, filter to it before training:
```bash
qiime taxa filter-seqs \
    --i-sequences eukaryome_its/eukaryome_sequences/ITS_seqs.qza \
    --i-taxonomy eukaryome_its/eukaryome_taxonomy/ITS_taxa.qza \
    --p-include k__Fungi \
    --o-filtered-sequences eukaryome-its-fungi-seqs.qza
```
Then train against `eukaryome-its-fungi-seqs.qza` and the original (unfiltered) taxonomy artifact
-- `fit-classifier-naive-bayes` ignores taxonomy entries that don't correspond to any reference
sequence, so there's no need to filter the taxonomy artifact too.

## From GenBank (NCBI query or accession list)

`qiime2-its-train-ncbi` downloads sequences, downloads/parses taxonomy, and trains the classifier.
It accepts either an NCBI query string or a text file with one accession number per line. Test your
query on the NCBI website first to make sure it returns the sequences you expect.
```bash
# Using a query string
qiime2-its-train-ncbi \
    -q "txid4762[Organism:exp] AND (\"internal transcribed spacer\"[Title]) NOT uncultured[Title]" \
    -t 48 \
    -o /oomycetes_DB \
    -e your_email@provider.org \
    -a your_ncbi_api_key_optional_but_recommended

# Using a text file of accessions
qiime2-its-train-ncbi \
    -q ~/Downloads/accession.list \
    -t 48 \
    -o /oomycetes_DB
```
`--acc2taxid`/`--dead-acc2taxid` let you point at already-downloaded copies of NCBI's
`nucl_gb.accession2taxid.gz`/`dead_nucl.accession2taxid.gz` (each several GB) instead of
downloading them again.

## From your own fasta file

`qiime2-its-train-fasta` builds a classifier from a fasta file and a matching "accession to taxid"
table (two tab-separated columns; accessions must match exactly between the fasta headers and the
table, up to the first whitespace):
```bash
qiime2-its-train-fasta \
    -q my_sequences.fasta \
    -i acc2taxid_table.tsv \
    -o /path/to/output_folder/
```

## Other databases (not for ITS)

`rescript` and the wider QIIME2 ecosystem also integrate several other reference databases:
[SILVA](https://www.arb-silva.de/) (`qiime rescript get-silva-data`, and SILVA also publishes
[ready-to-use, pretrained `.qza` classifiers](https://www.arb-silva.de/documentation/classifiers/qiime-2)
so you don't even need to train one yourself), GTDB (`get-gtdb-data`), PR2 (`get-pr2-data`), and a
few others. **None of these cover the ITS region** -- confirmed against each action's real
`--help` output: SILVA's `--p-target` only offers `SSURef`/`LSURef` (16S/18S and 23S/28S rRNA),
GTDB is bacterial/archaeal SSU only, PR2 is protist SSU only. This isn't a QIIME2 limitation to
work around -- the ITS region is the *non-conserved* spacer between rRNA genes, which is exactly
why these curated rRNA databases exclude it, and why UNITE/Eukaryome (both amplicon-database
projects, not rRNA-alignment projects) exist as separate resources for it in the first place.

That said, **this pipeline isn't hardcoded to ITS**: `--extract-its1`/`--extract-its2` are
optional (they just run ITSxpress before DADA2). Omit both and supply a matching classifier --
e.g. one of SILVA's pretrained SSU/LSU classifiers above -- and the rest of the pipeline (DADA2,
phylogeny, diversity, taxonomy, the advanced-stats/report features) works the same way for a
16S/18S/23S/28S rRNA amplicon dataset as it does for ITS.

Next: [Pipeline usage](Pipeline-Usage).
