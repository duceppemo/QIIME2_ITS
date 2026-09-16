# Troubleshooting / FAQ

**"File name must be as is..." but my fastq files look right to me.**
Check for an underscore *inside* the sample identifier itself, e.g. `siteA_rep1_S1_L001_R1_001.fastq.gz`.
The naming scheme splits on every underscore and needs exactly 5 fields — use a hyphen instead
(`siteA-rep1_...`). See [Metadata and fastq naming requirements](Metadata-and-FASTQ-Requirements).

**"bbduk.sh not found on PATH."**
`--min-len`/`--max-len` need BBTools/BBMap, which the QIIME2 environment file does not install:
`conda install -c bioconda -c conda-forge bbmap`. See [Installation](Installation).

**"The following sample(s) have an empty (zero-read) fastq file."**
A sample with zero reads is rejected upfront with a clear message rather than failing several
steps into ITSxpress with a cryptic external-tool error. Remove it from the input folder and the
metadata file.

**UNITE download fails / the old `-u <url>` example doesn't work anymore.**
UNITE's site now gates release files behind DOI/PlutoF landing pages with no stable direct-download
URL. Use `qiime rescript get-unite-data` instead — see [Building a classifier](Building-a-Classifier).

**Everything classifies as "unidentified" (or stops at Kingdom/Phylum).**
This is very likely a real result, not a bug: your classifier doesn't have the actual organisms in
your reads in its training set, so `classify-sklearn`'s confidence-based rank truncation correctly
stops at whatever rank it's actually confident about instead of guessing a specific wrong species.
Check that your classifier (UNITE is the broadest, most general-purpose option) actually covers the
taxa you expect.

**"You must activate your QIIME2 conda environment to run this script."**
`qiime2-its` and friends check `$CONDA_DEFAULT_ENV` for something containing `qiime2` (matches both
legacy `qiime2-*` and current `rachis-qiime2-*` naming) before doing anything else.
`conda activate <your-qiime2-env>` first.

**A sample-classifier column got skipped ("at least one class has fewer than 2 samples").**
`qiime sample-classifier classify-samples` needs every class (value) of the target column to have
≥2 samples for its internal stratified train/test split, even if the column looked balanced in your
metadata file — a near-empty sample that DADA2 reduced to a zero-read row can turn a nominally-
2-per-group column into an effective singleton. This is expected, not a bug: see
[Advanced stats and the PDF report](Advanced-Stats-and-Report).

**"Requested level of 6 is larger than the maximum level available in taxonomy data."**
Doesn't happen anymore — `qiime2-its` caps the genus-level (`taxa collapse`) request to whatever
depth the classifier actually resolved for your dataset. If you still see this calling
`qiime taxa collapse` yourself, your classifier didn't reach genus for any feature; check the
`taxonomy.qzv`/`biom_table/taxonomy.tsv` to see how deep your actual assignments go.

**Still stuck?** [Open an issue](https://github.com/duceppemo/QIIME2_ITS/issues).
