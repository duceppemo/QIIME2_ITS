# Installation

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

4. Only if you plan to use `--min-len`/`--max-len` (read-length filtering) -- or want to run
   `validation/run_validation.sh` (see [Validation suite](Validation-Suite)), which exercises that
   path -- install BBTools/BBMap, which provides `bbduk.sh`. It is **not** installed by the QIIME2
   environment file and `qiime2-its` will refuse to start with `--min-len`/`--max-len` until it's on
   `PATH`:
   ```bash
   conda install -c bioconda -c conda-forge bbmap
   ```

Next: [build a classifier](Building-a-Classifier), then see [Pipeline usage](Pipeline-Usage).
