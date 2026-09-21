# Worked example: a public SRA dataset, start to finish

Want to see what `qiime2-its` produces before installing anything? This page walks through a real,
production-scale run on public data -- 53 soil-fungi samples from NCBI BioProject
[PRJNA767765](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA767765) -- and publishes everything:
the scripts, the metadata, and the resulting report.

**[Open the full PDF report (37 pages)](https://github.com/duceppemo/QIIME2_ITS/blob/master/examples/PRJNA767765/report.pdf)** ·
[all the files](https://github.com/duceppemo/QIIME2_ITS/tree/master/examples/PRJNA767765)

## A look inside the report

| | |
|---|---|
| ![Alpha diversity boxplots](https://raw.githubusercontent.com/duceppemo/QIIME2_ITS/master/examples/PRJNA767765/screenshots/alpha_boxplots.png) | ![PCoA](https://raw.githubusercontent.com/duceppemo/QIIME2_ITS/master/examples/PRJNA767765/screenshots/pcoa.png) |
| **Alpha diversity** per group, one page for each metadata column that tested significant. | **PCoA** for Bray-Curtis and unweighted UniFrac. Colorblind-safe colors, plus marker shapes when a column has more groups than the palette has colors. |
| ![Sample clustering](https://raw.githubusercontent.com/duceppemo/QIIME2_ITS/master/examples/PRJNA767765/screenshots/dendrogram.png) | ![Genus composition](https://raw.githubusercontent.com/duceppemo/QIIME2_ITS/master/examples/PRJNA767765/screenshots/genus_composition.png) |
| **UPGMA sample clustering** from the same distance matrix, leaves colored like the PCoA. | **Genus-level composition** of every sample. |

![PERMANOVA table](https://raw.githubusercontent.com/duceppemo/QIIME2_ITS/master/examples/PRJNA767765/screenshots/permanova_table.png)

The report also carries the run's full provenance (exact command, parameters, tool and plugin
versions, input files), DADA2 read retention per sample, sequence-length and
classification-confidence distributions, rarefaction curves and the sample-classifier results --
see [Advanced stats and the PDF report](Advanced-Stats-and-Report) for what each page means.

## The dataset

| | |
|---|---|
| Source | BioProject PRJNA767765 -- fungal ITS amplicons, Illumina MiSeq, paired-end |
| Samples | the 53 field samples from Poland and Iran (the BioProject's laboratory microcosm samples are left out) |
| Size | 106 fastq files, 3.3 GB |
| Metadata | country, environmental medium, local environment, collection date, elevation, replicate -- built from the samples' public BioSample records |
| Classifier | UNITE 2025-02-19, fungi, 99% clusters, no singletons |

## Reproduce it

Inside an activated QIIME2 environment with `qiime2-its` installed ([Installation](Installation)),
from an empty working folder with ~15 GB free:

```bash
EX=/path/to/QIIME2_ITS/examples/PRJNA767765

python3 $EX/01_fetch_metadata.py .   # optional: rebuilds metadata.tsv + download_manifest.tsv
bash $EX/02_download_reads.sh        # 3.3 GB from ENA into ./raw_reads, MD5-verified, resumable
bash $EX/03_train_classifier.sh      # UNITE classifier into ./classifier (slow; skip if you have one)
THREADS=40 PARALLEL=10 bash $EX/04_run_pipeline.sh   # -> ./output/report.pdf
```

The analysis itself is a single command, everything else left at its default:

```bash
qiime2-its -q rachis-qiime2-2026.7 \
    -i raw_reads -o output -m metadata.tsv -c classifier/unite-classifier.qza \
    -pe --extract-its2 --taxa Fungi -t 40 -p 10
```

The published run used `qiime2-its` 0.3.2 on QIIME2 2026.7 and took **3 h 28 min** with 40 threads
on a 64-core workstation. Expect small numerical differences if you re-run it (rarefaction and
PERMANOVA permutations are random), not different conclusions.

| Script | What it does |
|---|---|
| [`01_fetch_metadata.py`](https://github.com/duceppemo/QIIME2_ITS/blob/master/examples/PRJNA767765/01_fetch_metadata.py) | Builds the QIIME2 metadata file and the download manifest from two public APIs (ENA file report, NCBI BioSample); standard library only. Its output is committed, so this step is optional -- it exists so nothing in the example is hand-made. |
| [`02_download_reads.sh`](https://github.com/duceppemo/QIIME2_ITS/blob/master/examples/PRJNA767765/02_download_reads.sh) | Downloads each run's fastq pair from ENA under the [Casava-style names](Metadata-and-FASTQ-Requirements) QIIME2 requires, checking every file against ENA's MD5. |
| [`03_train_classifier.sh`](https://github.com/duceppemo/QIIME2_ITS/blob/master/examples/PRJNA767765/03_train_classifier.sh) | Fetches UNITE with QIIME2's `rescript` plugin and trains the naive-Bayes classifier (see also [Building a classifier](Building-a-Classifier)). |
| [`04_run_pipeline.sh`](https://github.com/duceppemo/QIIME2_ITS/blob/master/examples/PRJNA767765/04_run_pipeline.sh) | The `qiime2-its` command above. |

[`metadata.tsv`](https://github.com/duceppemo/QIIME2_ITS/blob/master/examples/PRJNA767765/metadata.tsv) is also a realistic template for your own metadata file.

## One disclosed edit

Before publishing, the user name, host name and the home-directory part of the program path were
blanked in the run's `run_metadata.json`, and the PDF was rebuilt from the same output folder
(`qiime2_its.report.build_report()`), so "Run by" reads "(not recorded)". Nothing else was changed.
You can do the same before sharing your own reports.

## Data and citation

The sequencing data belongs to its authors; the repository only holds scripts and metadata derived
from the public archive records. If you use it, cite the study registered on the BioProject:

> Okrasińska A. *et al.* (2022). Marginal lands and fungi -- linking the type of soil contamination
> with fungal community composition. *Environmental Microbiology*.
> <https://doi.org/10.1111/1462-2920.16007>

UNITE is distributed under CC BY-SA 4.0 -- see <https://unite.ut.ee/cite.php>.
