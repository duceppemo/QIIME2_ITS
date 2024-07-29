# QIIME2_ITS
QIIME2 pipeline for single-end (IonTorrent) and paired-end (Illumina) sequencing data.

## Description
This pipeline uses QIIME2 to process metagenomics data. Uses ITSxpress to extract ITS1 sequences and DADA2 for denoising. Only basic analysis is performed.

## Important
Please make sure you are using a validated QIIME2 metadata tsv file. Use this tool to proceed with the validation step:
https://keemei.qiime2.org/

QIIME2 also requires that the input fastq files to have a specific naming scheme, used by Illumina sequencers. Valid QIIME fastq files names are structure like this:
```
'.+_.+_L[0-9][0-9][0-9]_R[12]_001\\.fastq\\.gz'
An example would be "L2S357_15_L001_R1_001.fastq.gz". The underscore-separated fields in this file name are:
  1. the sample identifier,
  2. the barcode sequence or a barcode identifier,.
  3. the lane number starting with "L" followed by 3 digits,
  4. the direction of the read ("R1" or "R2"; use R1 for single-end reads), and
  5. the set number (always "001").
```
If metadata file and fastq files are not properly formatted, this pipeline will crash.

Please make sure that your samples are already demultiplexed, i.e. you have one fastq file per sample.

Note that the DADA2 QIIME2 plugin is meant to process Illumina data. When running DADA2 in standalone mode, settings can be tweaked to compensate for IonTorrent errors (https://benjjneb.github.io/dada2/faq.html#can-i-use-dada2-with-my-454-or-ion-torrent-data), but those options are not available in the QIIME2 plugin. You may want to compare both to validate this pipeline for you applications.

Lastly, I only had access to limited data sets to write this pipeline. It's quite possible that you may encounter problems running it. If that is the case, please let me know by reporting an issue.
## Installation

1. Make sure you have conda installed. See https://docs.conda.io/en/latest/miniconda.html for `miniconda` installation instructions. You may prefer to install `anaconda`. `mamba` can also be installed afterward to speedup environment creation and package installation.

2. Install `qiime2`. See https://docs.qiime2.org/2021.2/install/native/#install-qiime-2-within-a-conda-environment for instructions. If you're using a newer version, make sure to install the "QIIME 2 Amplicon Distribution" and follow the updated instructions for "Install QIIME 2 within a `conda` environment"
3. Run the following command lines:
```bash
# Activate qiime2 conda environment
conda activate qiime2-2021.2  # Version number may differ depending of time of installation.

# Install ITSxpress and NCBI's entrez
conda install -c bioconda itsxpress entrez-direct

# Clone this repository (make sure `git` is installed in your environment. Use `conda install git` if needed):
git clone https://github.com/duceppemo/QIIME2_ITS

# Go to cloned repository
cd QIIME2_ITS

# Test pipeline:
python qiime2_its.py -h
```
## Database
In order to use this pipeline, you have to have a qiime2 classifier available. Note that if you update your QIIME2 version, you might need to recompile the classifier.

Here are the instruction to build the UNITE qiime2 classifier using the included helper script `train_unite_classifier_qiime2.py`. Use `python train_unite_classifier_qiime2.py -h` for more options.
```bash
# QIIME files for unite are located at https://unite.ut.ee/repository.php under the "QIIME release" drop down menu.
# Replace the url ("-u"), the datase location ("-o") and the qiime2 ("-q") to suit your installation. If you have
# already downloaded the QIIME release file, you can substitute the URL with the actual ".tgz" file.

conda activate qiime2-2021.2  # If you didn't alread activated it

# Using the URL
python train_unite_classifier_qiime2.py \
    -u https://files.plutof.ut.ee/public/orig/C5/54/C5547B97AAA979E45F79DC4C8C4B12113389343D7588716B5AD330F8BDB300C9.tgz \
    -o '/db/UNITE' \
    -q qiime2-2021.2

# Using the downloaded file
python train_unite_classifier_qiime2.py \
    -u ~/Download/C5547B97AAA979E45F79DC4C8C4B12113389343D7588716B5AD330F8BDB300C9.tgz \
    -o '/db/UNITE' \
    -q qiime2-2021.2
```
Here's how to build a classifier with sequences hosted on GenBank with the script `train_ncbi_classifier_qiime2.py`. That script takes care of downloading the sequences, the taxonomy and train the QIIME2 classifier. It takes as input a NCBI query input or a text file with one accession number per line. It is strongly recommended testing your query string on NCBI's website first to make sure you get the right sequences. Use `python train_ncbi_classifier_qiime2.py -h` for more options.  
```bash
# Using query string:
python train_ncbi_classifier_qiime2.py \
    -q "txid4762[Organism:exp] AND (\"internal transcribed spacer\"[Title]) NOT uncultured[Title]" \
    -t 48 \
    -o /oomycetes_DB \
    -e your_email@provider.org \
    -a ncbiapikeyisoptionalebutrecommended0

# Using a text file with accessions:
python train_ncbi_classifier_qiime2.py \
    -q ~/Download/accession.list \
    -t 48 \
    -o /oomycetes_DB \
```
ere's how to build a classifier with sequences hosted on GenBank with the script `train_fasta_classifier_qiime2.py`. You need to supply your own fasta file and "Accession to taxid" table. The script will expand the taxonomy and train the QIIME2 classifier. Use `python train_fasta_classifier_qiime2.py -h` for more options.  

## Usage - qiime2.py
Don't forget to activate your environment, if not already done.
```commandline
usage: python qiime2_its.py [-h] -q qiime2-2020.8 -i /input_folder/ -o /output_folder/ -m qiime2_metadata.tsv -c unite_classifier_qiime2.qza [-t 4] [-p 1] [-rc] [--min-len MIN_LEN] [--max-len MAX_LEN] [-se] [-pe]
                            [--extract-its1] [--extract-its2] [--taxa Fungi]

Run QIIME2 on IonTorrent sequencing data using the UNITE database

optional arguments:
  -h, --help            show this help message and exit
  -q qiime2-2020.8, --qiime2 qiime2-2020.8
                        Name of your QIIME2 conda environment. Mandatory.
  -i /input_folder/, --input /input_folder/
                        Input folder where the fastq reads are located. Mandatory.
  -o /output_folder/, --output /output_folder/
                        Output folder for QIIME2 files. Mandatory.
  -m qiime2_metadata.tsv, --metadata qiime2_metadata.tsv
                        Validated QIIME2 metadata file (samples description). Mandatory.
  -c unite_classifier_qiime2.qza, --classifier unite_classifier_qiime2.qza
                        Classifier for QIIME2. See script "train_unite_classifier_qiime2.py" to compile it. Mandatory.
  -t 4, --threads 4     Number of CPU. Default is 4.
  -p 1, --parallel-processes 1
                        Processes to run in parallel. Adjust according the number of threads. For example, if 16 threads, using 4 parallel processes will run 4 samples in parallel using 4 threads each (16/4). Default
                        is 1.
  -rc, --reverse_complement
                        Use this flag is your reads are in reverse complement. For example if you sequenced from 5.8S to 18S. Optional.
  --min-len MIN_LEN     Minimum read length to keep. Default is 0 (no min length). Optional.
  --max-len MAX_LEN     Maximum read length to keep. Default is 0 (no max length). Optional.
  -se                   Reads are single-end. One fastq file per sample. "-se" or "-pe" mandatory.
  -pe                   Reads are paired-end. Two fastq file per sample. "-se" or "-pe" mandatory.
  --extract-its1        Extract ITS1 sequence from reads with ITSxpress. Cannot be used with "--extract-its2".
  --extract-its2        Extract ITS2 sequence from reads with ITSxpress. Cannot be used with "--extract-its1".
  --taxa Fungi          Select taxa of interest for ITSxpress:
                        {Alveolata,Bryophyta,Bacillariophyta,Amoebozoa,Euglenozoa,Fungi,Chlorophyta,Rhodophyta,Phaeophyceae,Marchantiophyta,Metazoa,Oomycota,Haptophyceae,Raphidophyceae,
                        Rhizaria,Synurophyceae,Tracheophyta,Eustigmatophyceae,All}
```

## Usage - train_ncbi_classifier_qiime2.py
```commandline
usage: python train_ncbi_classifier_qiime2.py [-h] -q "txid4762[Organism:exp] AND "internal transcribed spacer"[Title] NOT uncultured[Title]" -o /output_folder/ [-t 4] [-e your.email@example.org] [-a]
                                              [--taxdump /path/to/taxdump.tar.gz] [--acc2taxid /path/to/nucl_gb.accession2taxid.gz] [--dead-acc2taxid /path/to/dead_nucl.accession2taxid.gz]

Download DNA sequence from NCBI and add taxonomy for QIIME2.

optional arguments:
  -h, --help            show this help message and exit
  -q "txid4762[Organism:exp] AND ("internal transcribed spacer"[Title]) NOT uncultured[Title]", --query "txid4762[Organism:exp] AND ("internal transcribed spacer"[Title]) NOT uncultured[Title]"
                        NCBI query string OR a text file with one accession number per line. Mandatory.
  -o /output_folder/, --output /output_folder/
                        Output folder. Mandatory.
  -t 4, --threads 4     Number of CPU. Default is 4. Optional.
  -e your.email@example.org, --email your.email@example.org
                        Your email address. Optional.
  -a , --api-key        Your NCBI API key. Allows up to 10 requests per second instead of 3. Optional.
  --taxdump /path/to/taxdump.tar.gz
                        Path to downloaded taxdump.tar.gz. Optional.
  --acc2taxid /path/to/nucl_gb.accession2taxid.gz
                        Path to downloaded nucl_gb.accession2taxid.gz. Optional.
  --dead-acc2taxid /path/to/dead_nucl.accession2taxid.gz
                        Path to downloaded dead_nucl.accession2taxid.gz. Optional.
```
## Usage - train_unite_classifier_qiime2.py
```commandline
usage: python train_unite_classifier_qiime2.py [-h] -u http://www.fileserver.com/file.txt -o /unite_folder/ -q qiime2-2020.8 [-t 16]

Train Unite classifier for QIIME2

optional arguments:
  -h, --help            show this help message and exit
  -u http://www.fileserver.com/file.txt, --url http://www.fileserver.com/file.txt
                        URL of UNITE QIIME2 database (compressed tar file with sequences and taxonomy) OR the path the the file already downloaded. Mandatory.
  -o /unite_folder/, --output_folder /unite_folder/
                        Output folder for classifier. Mandatory.
  -q qiime2-2020.8, --qiime2 qiime2-2020.8
                        Name of your QIIME2 conda environment. Mandatory.
  -t 16, --threads 16   Number of CPU. Default is 16

```
## Usage - train_fasta_classifier_qiime2.py
```commandline
usage: python train_fasta_classifier_qiime2.py [-h] -q my_sequences.fasta -i /path/to/acc2taxid_table.tsv -o /path/to/output_folder/ [--taxdump /path/to/taxdump.tar.gz]

Prep a QIIME2 classifier from a fasta file and corresponding "acc to taxid" table.

optional arguments:
  -h, --help            show this help message and exit
  -q my_sequences.fasta, --query my_sequences.fasta
                        A fasta file. Mandatory.
  -i /path/to/acc2taxid_table.tsv, --id-table /path/to/acc2taxid_table.tsv
                        Tab-separated text file with 2 columns (accession + taxid) matching your input fasta. The accession numbers (i.e. everything befor the first whitespace in the header of sequencers, must match
                        exactly between the fasta and the table.Mandatory.
  -o /path/to/output_folder/, --output /path/to/output_folder/
                        Output folder. Mandatory.
  --taxdump /path/to/taxdump.tar.gz
                        Path to downloaded taxdump.tar.gz. Will be downloaded otherwise. Optional.
```
