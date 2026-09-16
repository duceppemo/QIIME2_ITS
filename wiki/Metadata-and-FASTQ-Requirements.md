# Metadata and fastq naming requirements

## Metadata file

Use a validated QIIME2 metadata TSV file. Validate it with
[Keemei](https://keemei.qiime2.org/) (a Google Sheets add-on maintained by the QIIME2 project).

Sample identifiers in the metadata file must **exactly match** the sample identifier field in the
fastq file names (see below).

## FASTQ naming

QIIME2 requires the input fastq files to follow the naming scheme used by Illumina sequencers:
```
'.+_.+_L[0-9][0-9][0-9]_R[12]_001\.fastq\.gz'
```
An example: `L2S357_15_L001_R1_001.fastq.gz`. The underscore-separated fields are:
1. the sample identifier
2. the barcode sequence or a barcode identifier
3. the lane number, starting with `L` followed by 3 digits
4. the direction of the read (`R1` or `R2`; use `R1` for single-end reads)
5. the set number (always `001`)

**The sample identifier and barcode fields must not themselves contain underscores.** File names
are split on every underscore and must yield exactly 5 fields, so e.g.
`siteA_rep1_S1_L001_R1_001.fastq.gz` (sample identifier `siteA_rep1`) is rejected — use
`siteA-rep1_S1_L001_R1_001.fastq.gz` instead. `qiime2-its` checks this upfront and fails with a
clear message before running anything, rather than partway through.

If the metadata file and fastq files aren't properly formatted, the pipeline raises an error before
doing any work.

Your samples must already be demultiplexed: one fastq file (or one pair, for paired-end) per
sample. A sample with a genuinely empty (zero-read) fastq file is also rejected upfront, by name,
before the pipeline starts — remove it from the input folder and the metadata file first.

## A note on DADA2 and IonTorrent

The DADA2 QIIME2 plugin is built for Illumina data. Standalone DADA2 has settings to compensate for
IonTorrent-specific errors
([see the DADA2 FAQ](https://benjjneb.github.io/dada2/faq.html#can-i-use-dada2-with-my-454-or-ion-torrent-data)),
but those options aren't exposed by the QIIME2 plugin. What *is* exposed and useful for noisier
single-end platforms: `--max-ee` and `--allow-one-off` (see [Pipeline usage](Pipeline-Usage)). You
may want to compare against standalone DADA2 to validate this pipeline for your application.

This pipeline was written and validated against a limited number of real datasets (see
[Validation suite](Validation-Suite)). If you run into problems, please
[report an issue](https://github.com/duceppemo/QIIME2_ITS/issues).
