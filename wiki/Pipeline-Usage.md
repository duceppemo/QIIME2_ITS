# Pipeline usage — `qiime2-its`

Don't forget to activate your QIIME2 environment first. See
[Metadata and fastq naming requirements](Metadata-and-FASTQ-Requirements) before your first run.

```
usage: qiime2-its [-h] -q rachis-qiime2-2026.7 -i /input_folder/ -o
                  /output_folder/ -m qiime2_metadata.tsv -c
                  unite_classifier_qiime2.qza [-t 4] [-p 1] [-rc]
                  [--min-len MIN_LEN] [--max-len MAX_LEN] (-se | -pe)
                  [--extract-its1 | --extract-its2] [--taxa Fungi]
                  [--max-ee 2.0] [--max-ee-r 2.0] [--trunc-q 2]
                  [--pooling-method {independent,pseudo}]
                  [--chimera-method {consensus,none}]
                  [--min-fold-parent-over-abundance 1.0] [--allow-one-off]
                  [--n-reads-learn 1000000] [--sampling-depth 1000]
                  [--max-rarefaction-depth 4000] [--skip-advanced-stats]
                  [--skip-report] [--report-metadata-column COLUMN]
                  [--version]
```

## Core options

| Flag | Description |
|---|---|
| `-q`, `--qiime2` | Name of your QIIME2 conda environment. Mandatory. |
| `-i`, `--input` | Input folder with the fastq reads. Mandatory. |
| `-o`, `--output` | Output folder for QIIME2 files. Mandatory. |
| `-m`, `--metadata` | Validated QIIME2 metadata file. Mandatory. |
| `-c`, `--classifier` | Classifier `.qza` — see [Building a classifier](Building-a-Classifier). Mandatory. |
| `-t`, `--threads` | Number of CPUs. Default 4. |
| `-p`, `--parallel-processes` | Samples to process in parallel (e.g. 16 threads / 4 processes = 4 threads/sample). Default 1. |
| `-rc`, `--reverse_complement` | Reads are in reverse complement (e.g. sequenced 5.8S to 18S). |
| `--min-len` / `--max-len` | Read-length filtering (needs `bbduk.sh` — see [Installation](Installation)). Default 0 (off). |
| `-se` / `-pe` | Single-end or paired-end. Exactly one required. |
| `--extract-its1` / `--extract-its2` | Extract ITS1 or ITS2 with ITSxpress. At most one. Both optional — omit both to run the pipeline on a non-ITS marker (16S/18S/etc.) with a matching classifier, see [Building a classifier](Building-a-Classifier#other-databases-not-for-its). |
| `--taxa` | Taxon for ITSxpress. Default `Fungi`. Full list: `qiime itsxpress trim-single --help`. |

## DADA2 denoising

QIIME2's DADA2 defaults are tuned for Illumina data. For noisier single-end platforms (e.g.
IonTorrent), consider raising `--max-ee` and passing `--allow-one-off`: their higher per-base error
rate means more real reads get discarded by the default max-expected-errors threshold, and
homopolymer-driven indels get miscalled as one-off bimeras by the default chimera search.

| Flag | Description |
|---|---|
| `--max-ee` | Max expected errors (forward read, for paired-end) before a read is discarded. Default 2.0. |
| `--max-ee-r` | Max expected errors for the reverse read (paired-end only). Defaults to `--max-ee`. |
| `--trunc-q` | Truncate reads at the first quality score at or below this value. Default 2. |
| `--pooling-method` | `independent` or `pseudo`. Default `independent`. |
| `--chimera-method` | `consensus` or `none`. Default `consensus`. |
| `--min-fold-parent-over-abundance` | Chimera parent abundance fold-change threshold. Default 1.0. |
| `--allow-one-off` | Also flag one-mismatch/indel-from-exact bimeras as chimeric. Off by default. |
| `--n-reads-learn` | Minimum reads used to train the DADA2 error model. Default 1000000 (lower for faster runs on small datasets). |

### Recommended starting point for IonTorrent (single-end) data

IonTorrent's per-base error profile is dominated by homopolymer indels rather than Illumina-style
substitutions, so the Illumina-tuned defaults above throw away a lot of otherwise-good reads and
over-call chimeras. A reasonable starting point, and the one this pipeline's own
[validation suite](Validation-Suite) uses against real single-end data:

```bash
qiime2-its \
    -q <your-qiime2-env> \
    -i /input_folder/ -o /output_folder/ \
    -m qiime2_metadata.tsv -c classifier.qza \
    -se --extract-its2 --taxa Fungi \
    --max-ee 4 --allow-one-off
```

- **`--max-ee 4`** (up from the default 2.0): IonTorrent's higher raw error rate means the default
  threshold discards a large fraction of otherwise-usable reads before denoising ever sees them.
  Treat 4 as a starting point, not a universal value — check the `denoising-stats.qzv` retention
  numbers after a first run and raise further if too many reads are being dropped at the filtering
  step specifically (as opposed to the denoising/chimera steps).
- **`--allow-one-off`**: without it, a true ASV that differs from a more abundant parent by exactly
  one homopolymer-length indel — routine on IonTorrent — is *not* checked as a one-off bimera and so
  is more likely to survive as a spurious extra ASV. Turning it on makes the chimera search also
  flag those one-mismatch/indel bimeras, trading a little sensitivity for a lot fewer homopolymer-
  indel artifacts being called as real variants.
- Leave `--trunc-q`, `--pooling-method`, `--chimera-method`, and
  `--min-fold-parent-over-abundance` at their defaults unless you have a specific reason: they're not
  platform-specific in the same way.
- There is no length-truncation flag for single-end data in this pipeline (`qiime dada2
  denoise-single` truncates by position, and IonTorrent read lengths are naturally variable, so
  truncating by position is generally the wrong tool here). Use `--min-len`/`--max-len` instead if
  you need to filter by read length — see the table above.

### IonTorrent-specific DADA2 options that QIIME2 does not expose

Standalone DADA2 (the R package) has two additional options, set via `setDadaOpt()`, that its own
["Big Data: Pyrosequencing"](https://benjjneb.github.io/dada2/bigdata_paired.html) workflow
specifically recommends for 454/IonTorrent data:

- `HOMOPOLYMER_GAP_PENALTY` — a separate, more lenient gap penalty specifically for homopolymer
  gaps during alignment (DADA2's tutorial suggests `-1`).
- `BAND_SIZE` — widens the alignment band used when comparing reads (DADA2's tutorial suggests `32`),
  needed because IonTorrent's indel-driven errors shift alignments further than Illumina's
  substitution-driven ones.

**`qiime dada2 denoise-single`/`denoise-paired` do not expose either option** (confirmed against
`qiime dada2 denoise-single --help` in `rachis-qiime2-2026.7`) — QIIME2's q2-dada2 plugin only
surfaces the parameter subset listed in the table above, and neither
`--p-homopolymer-gap-penalty` nor `--p-band-size` (or any equivalent) exists as a CLI flag. There is
no workaround within this pipeline or within `qiime dada2` itself; getting that level of control
means running DADA2's R package directly outside QIIME2 (`setDadaOpt(HOMOPOLYMER_GAP_PENALTY=-1,
BAND_SIZE=32)` before `dada()`) and importing the resulting ASV table back into QIIME2, which is
outside what this pipeline automates. In practice, `--max-ee`/`--allow-one-off` recover most of the
same benefit for typical ITS-amplicon-scale datasets; reach for the R package directly only if you
still see excess homopolymer-driven spurious ASVs after tuning those.

## Diversity analysis

| Flag | Description |
|---|---|
| `--sampling-depth` | Rarefaction depth for `core-metrics-phylogenetic`. Samples with fewer reads are excluded. Default 1000 — lower this for small/pilot datasets. |
| `--max-rarefaction-depth` | Maximum depth for the alpha-rarefaction plot. Default 4000. |

## Advanced stats and report

Runs by default after the core pipeline — see [Advanced stats and the PDF report](Advanced-Stats-and-Report)
for what these actually produce.

| Flag | Description |
|---|---|
| `--skip-advanced-stats` | Skip group-significance tests, taxonomy collapse, and sample classification. |
| `--skip-report` | Skip building the PDF summary report. |
| `--report-metadata-column` | Categorical column the report's alpha/beta plots are grouped by. Defaults to the first eligible (≥2 distinct values, ≥2 samples each) column. Group-significance/classification still run against every eligible column regardless of this choice. |

## `qiime2-its-rc`

Reverse-complements every fastq file (gzipped or not) found recursively under an input folder:
```
usage: qiime2-its-rc [-h] -i /input_folder/ -o /modified_fastq/ [-t N]
```

Next: [Advanced stats and the PDF report](Advanced-Stats-and-Report).
