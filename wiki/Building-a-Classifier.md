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

Next: [Pipeline usage](Pipeline-Usage).
