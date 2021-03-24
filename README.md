# IonTorrent_Fungi_Barcoding_QIIME2
QIIME2 pipeline for IonTorrent sequencing data

## Important
Please make sure you are using a validated QIIME2 metadata tsv file. Use this tool to proceed with the validation step:
https://keemei.qiime2.org/

QIIME2 also requires that the input fastq files to have a specific naming scheme, used by Illimina sequencers. Valid QIIME fastq files names are structure like this:
```
L2S357_15_L001_R1_001.fastq.gz. The underscore-separated fields in this file name are:
  1. the sample identifier,
  2. the barcode sequence or a barcode identifier,
  3. the lane number,
  4. the direction of the read (i.e. only R1, because these are single-end reads), and
  5. the set number.
```
If metadata file and fastq files are not properly formated, this pipeline will crash.

## Installation

1. Make sure you have conda installed. See https://docs.conda.io/en/latest/miniconda.html for `miniconda` installation instructions. You may prefer to install `anaconda`. `mamba` can also be installed afterward to speedup environment creation and package installation.
2. Install `qiime2`. See https://docs.qiime2.org/2021.2/install/native/#install-qiime-2-within-a-conda-environment for instructions.
3. Activate qiime2 conda environment.
```
conda activate qiime2-2021.2  # Version number may differ depending of time of installation.
```
4. Clone this repository (make sure `git` is installed in your environment. Use `conda install git` if needed):
```
git clone https://github.com/duceppemo/IonTorrent_Fungi_Barcoding_QIIME2
```
5. `cd` into cloned repository:
```
cd IonTorrent_Fungi_Barcoding_QIIME2
```
6. Test pipeline:
```
python run_qiime2.py -h
```

## Usage
```
usage: run_qimme2.py [-h] -q qiime2-2020.8 -i /input_folder/ -o
                     /output_folder/ -m qiime2_metadata.tsv -u
                     unite_classifier_qiime2.qza [-t 16] [-rc]
                     [--min_len MIN_LEN] [--max_len MAX_LEN]

Run QIIME2 on IonTorrent sequencing data using the UNITE database

optional arguments:
  -h, --help            show this help message and exit
  -q qiime2-2020.8, --qiime2 qiime2-2020.8
                        Name of your QIIME2 conda environment.Mandatory.
  -i /input_folder/, --input /input_folder/
                        Input folder where the fastq reads are
                        located.Mandatory.
  -o /output_folder/, --output /output_folder/
                        Output folder for QIIME2 files.Mandatory.
  -m qiime2_metadata.tsv, --metadata qiime2_metadata.tsv
                        Validated QIIME2 metadata file (samples
                        description).Mandatory.
  -u unite_classifier_qiime2.qza, --unite unite_classifier_qiime2.qza
                        UNITE classifier for QIIME2. See script
                        "train_unite_classifier_qiime2.py" to compile
                        it.Mandatory.
  -t 16, --threads 16   Number of CPU. Default is 16
  -rc, --reverse_complement
                        Use this flag is your reads are in reverse complement.
                        For example if you sequenced from 5.8S to 18S.
  --min_len MIN_LEN     Minimum read length to keep. Default is 0 (no min
                        length).
  --max_len MAX_LEN     Maximum read length to keep. Default is 0 (no max
                        length).

```
