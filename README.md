# QIIME2_ITS

QIIME2 pipeline for single-end (IonTorrent) and paired-end (Illumina) ITS metabarcoding data.

## Description
This pipeline uses QIIME2 to process metagenomics data. It uses the `qiime itsxpress` plugin to
extract the ITS1/ITS2 region and DADA2 for denoising. Only basic analysis is performed.

## Important
Please make sure you are using a validated QIIME2 metadata TSV file. Use this tool to validate it:
https://keemei.qiime2.org/

QIIME2 also requires the input fastq files to follow the naming scheme used by Illumina sequencers:
```
'.+_.+_L[0-9][0-9][0-9]_R[12]_001\.fastq\.gz'
```
An example would be `L2S357_15_L001_R1_001.fastq.gz`. The underscore-separated fields are:
  1. the sample identifier,
  2. the barcode sequence or a barcode identifier,
  3. the lane number, starting with "L" followed by 3 digits,
  4. the direction of the read ("R1" or "R2"; use R1 for single-end reads), and
  5. the set number (always "001").

If the metadata file and fastq files are not properly formatted, this pipeline will raise an error.

Please make sure your samples are already demultiplexed, i.e. you have one fastq file (or one pair,
for paired-end) per sample.

Note that the DADA2 QIIME2 plugin is meant to process Illumina data. When running DADA2 in standalone
mode, settings can be tweaked to compensate for IonTorrent errors
(https://benjjneb.github.io/dada2/faq.html#can-i-use-dada2-with-my-454-or-ion-torrent-data), but those
options are not available in the QIIME2 plugin. You may want to compare both to validate this pipeline
for your application.

Lastly, this pipeline was written against a limited number of datasets. If you run into problems,
please report an issue.

## Installation

1. Install `conda`/`mamba`. See https://docs.conda.io/en/latest/miniconda.html for `miniconda`
   installation instructions; `mamba` speeds up environment creation and package installation.

2. Install QIIME2. See https://library.qiime2.org for current instructions — QIIME2's amplicon
   distribution is now named `rachis-qiime2` (renamed from `qiime2-amplicon`, from the 2026.4
   release onward). For example, for the current release:
   ```bash
   VERSION=2026.7  # check https://library.qiime2.org/quickstart/qiime2 for the current release
   conda env create -n rachis-qiime2-$VERSION \
       --file https://raw.githubusercontent.com/qiime2/distributions/dev/$VERSION/qiime2/released/rachis-qiime2-linux-64-conda.yml
   ```

3. Activate the environment and install ITSxpress (this also registers the `qiime itsxpress` plugin
   used by this pipeline) plus this package:
   ```bash
   conda activate rachis-qiime2-2026.7   # name may differ depending on when you installed it
   conda install -c bioconda -c conda-forge itsxpress
   qiime dev refresh-cache                # picks up the newly installed itsxpress plugin

   git clone https://github.com/duceppemo/QIIME2_ITS
   cd QIIME2_ITS
   pip install -e .

   # Test the pipeline:
   qiime2-its -h
   ```

4. Only if you plan to use `--min-len`/`--max-len` (read-length filtering): install BBTools/BBMap,
   which provides `bbduk.sh`. It is **not** installed by the QIIME2 environment file and `qiime2-its`
   will refuse to start with `--min-len`/`--max-len` until it's on `PATH`:
   ```bash
   conda install -c bioconda -c conda-forge bbmap
   ```

## Database
You need a QIIME2 classifier to run this pipeline. If you update your QIIME2 version, you will
generally need to recompile the classifier.

### From UNITE
QIIME files for UNITE are at https://unite.ut.ee/repository.php under "QIIME release". Replace the
URL (`-u`), the database location (`-o`) and the QIIME2 env name (`-q`) to suit your installation.
You can also pass an already-downloaded `.tgz` file as `-u`.
```bash
conda activate rachis-qiime2-2026.7

# Using a URL
qiime2-its-train-unite \
    -u https://files.plutof.ut.ee/public/orig/C5/54/C5547B97AAA979E45F79DC4C8C4B12113389343D7588716B5AD330F8BDB300C9.tgz \
    -o /db/UNITE \
    -q rachis-qiime2-2026.7

# Using an already-downloaded file
qiime2-its-train-unite \
    -u ~/Downloads/C5547B97AAA979E45F79DC4C8C4B12113389343D7588716B5AD330F8BDB300C9.tgz \
    -o /db/UNITE \
    -q rachis-qiime2-2026.7
```

### From GenBank (NCBI query or accession list)
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

### From your own fasta file
`qiime2-its-train-fasta` builds a classifier from a fasta file and a matching "accession to taxid"
table (two tab-separated columns; accessions must match exactly between the fasta headers and the
table, up to the first whitespace).
```bash
qiime2-its-train-fasta \
    -q my_sequences.fasta \
    -i acc2taxid_table.tsv \
    -o /path/to/output_folder/
```

## Usage — qiime2-its
Don't forget to activate your environment first.
```
usage: qiime2-its [-h] -q rachis-qiime2-2026.7 -i /input_folder/
                   -o /output_folder/ -m qiime2_metadata.tsv
                   -c unite_classifier_qiime2.qza [-t 4] [-p 1] [-rc]
                   [--min-len MIN_LEN] [--max-len MAX_LEN]
                   (-se | -pe) [--extract-its1 | --extract-its2]
                   [--taxa Fungi]

Run QIIME2 on ITS amplicon data using ITSxpress and DADA2

options:
  -h, --help            show this help message and exit
  -q, --qiime2          Name of your QIIME2 conda environment. Mandatory.
  -i, --input           Input folder where the fastq reads are located. Mandatory.
  -o, --output          Output folder for QIIME2 files. Mandatory.
  -m, --metadata        Validated QIIME2 metadata file (samples description). Mandatory.
  -c, --classifier      Classifier for QIIME2 (see qiime2-its-train-unite/-ncbi/-fasta). Mandatory.
  -t, --threads         Number of CPUs. Default is 4.
  -p, --parallel-processes
                        Samples to process in parallel. E.g. with 16 threads, 4 parallel
                        processes will run 4 samples in parallel using 4 threads each. Default 1.
  -rc, --reverse_complement
                        Use if your reads are in reverse complement (e.g. sequenced 5.8S to 18S).
  --min-len MIN_LEN     Minimum read length to keep. Default 0 (no minimum).
  --max-len MAX_LEN     Maximum read length to keep. Default 0 (no maximum).
  -se                   Reads are single-end (one fastq file per sample).
  -pe                   Reads are paired-end (two fastq files per sample). "-se" or "-pe" is mandatory.
  --extract-its1        Extract the ITS1 region with ITSxpress.
  --extract-its2        Extract the ITS2 region with ITSxpress.
  --taxa Fungi          Taxon of interest for ITSxpress. See `qiime itsxpress trim-single --help`
                         for the full list.
```

## Usage — qiime2-its-rc
Reverse-complements every fastq file (gzipped or not) found recursively under an input folder.
```
usage: qiime2-its-rc [-h] -i /input_folder/ -o /modified_fastq/ [-t N]
```

## Development
```bash
pip install -e '.[dev]'
pytest -v
```
All unit tests mock external tool invocations (`qiime`, `biom`, `bbduk.sh`), so they run with a
plain Python interpreter — no QIIME2 installation required. CI runs them on every push (see
`.github/workflows/tests.yml`).

## Migrating from the pre-0.2 flat scripts
Versions before 0.2 shipped as standalone scripts (`python qiime2_its.py ...`). As of 0.2, this is
an installable package with console-script entry points instead: `qiime2-its`, `qiime2-its-rc`,
`qiime2-its-train-unite`, `qiime2-its-train-ncbi`, and `qiime2-its-train-fasta` (`pip install -e .`
from the repository root, inside your QIIME2 conda environment). Command-line flags are otherwise
unchanged. ITS extraction now goes through the `qiime itsxpress` plugin instead of the old standalone
`itsxpress` CLI — install it with `conda install -c bioconda -c conda-forge itsxpress` followed by
`qiime dev refresh-cache` inside your QIIME2 environment.
