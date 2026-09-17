<p align="center">
  <img src="https://raw.githubusercontent.com/duceppemo/QIIME2_ITS/master/docs/QIIME2_ITS_logo.png" alt="QIIME2_ITS logo" width="280">
</p>

<h1 align="center">QIIME2_ITS</h1>

<p align="center">
  <a href="https://github.com/duceppemo/QIIME2_ITS/actions/workflows/tests.yml"><img src="https://github.com/duceppemo/QIIME2_ITS/actions/workflows/tests.yml/badge.svg" alt="tests"></a>
  <a href="https://github.com/duceppemo/QIIME2_ITS/blob/master/LICENSE"><img src="https://img.shields.io/github/license/duceppemo/QIIME2_ITS" alt="license"></a>
  <img src="https://img.shields.io/badge/python-3.9%2B-blue" alt="python 3.9+">
  <img src="https://img.shields.io/badge/QIIME2-2026.x-blue" alt="QIIME2 2026.x">
</p>

<p align="center">
  QIIME2 pipeline for single-end (IonTorrent) and paired-end (Illumina) ITS metabarcoding data.
  ITSxpress for ITS region extraction, DADA2 for denoising, plus diversity statistics,
  classification, and a PDF summary report.
</p>

## Quick start

```bash
# 1. Install QIIME2 (see the wiki for details/troubleshooting):
VERSION=2026.7  # check https://library.qiime2.org/quickstart/qiime2 for the current release
conda env create -n rachis-qiime2-$VERSION \
    --file https://raw.githubusercontent.com/qiime2/distributions/dev/$VERSION/qiime2/released/rachis-qiime2-linux-64-conda.yml
conda activate rachis-qiime2-$VERSION

# 2. Install ITSxpress + BBMap into it, and this package:
conda install -c bioconda -c conda-forge itsxpress bbmap
qiime dev refresh-cache

git clone https://github.com/duceppemo/QIIME2_ITS
cd QIIME2_ITS
pip install -e .

# 3. Run it:
qiime2-its \
    -q rachis-qiime2-$VERSION \
    -i /input_folder/ -o /output_folder/ \
    -m qiime2_metadata.tsv -c classifier.qza \
    -pe --extract-its2 --taxa Fungi
```

That's the short version. **[See the wiki](https://github.com/duceppemo/QIIME2_ITS/wiki) for:**
full installation, building a classifier (UNITE/NCBI/your own fasta), the complete flag reference,
fastq/metadata requirements, the diversity-stats-and-PDF-report feature, the real-data validation
suite, and troubleshooting.

## License

[MIT](https://github.com/duceppemo/QIIME2_ITS/blob/master/LICENSE)
