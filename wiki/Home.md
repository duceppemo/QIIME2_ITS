<p align="center">
  <img src="https://raw.githubusercontent.com/duceppemo/QIIME2_ITS/master/docs/QIIME2_ITS_logo.png" alt="QIIME2_ITS logo" width="160">
</p>

# QIIME2_ITS wiki

QIIME2 pipeline for single-end (IonTorrent) and paired-end (Illumina) ITS metabarcoding data. Uses
the `qiime itsxpress` plugin to extract the ITS1/ITS2 region and DADA2 for denoising, then runs
phylogeny, diversity, taxonomic classification, and (by default) diversity/composition/classifier
statistics and a PDF summary report.

The [README](https://github.com/duceppemo/QIIME2_ITS#readme) covers a quick install and a minimal
run. Everything else lives here:

- **[Installation](Installation)** — full QIIME2 + ITSxpress + BBMap setup
- **[Building a classifier](Building-a-Classifier)** — UNITE, NCBI, or your own fasta file
- **[Pipeline usage](Pipeline-Usage)** — the full `qiime2-its` flag reference
- **[Metadata and fastq naming requirements](Metadata-and-FASTQ-Requirements)** — the Casava
  naming scheme, common gotchas
- **[Advanced stats and the PDF report](Advanced-Stats-and-Report)** — group-significance tests,
  genus-level composition, sample classification, `report.pdf`
- **[Validation suite](Validation-Suite)** — the real-data validation suite and its audit trail
- **[Development](Development)** — running the unit tests, project layout
- **[Troubleshooting / FAQ](Troubleshooting-FAQ)**

Issues and pull requests: https://github.com/duceppemo/QIIME2_ITS
