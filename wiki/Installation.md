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

3. Activate the environment and install ITSxpress (this also registers the `qiime itsxpress`
   plugin used by this pipeline) and BBTools/BBMap (provides `bbduk.sh`, used for
   `--min-len`/`--max-len` read-length filtering and by `validation/run_validation.sh` -- see
   [Validation suite](Validation-Suite)). Neither is installed by the QIIME2 environment file
   itself, so both are worth installing upfront rather than hitting a missing-tool error later:
   ```bash
   conda activate rachis-qiime2-2026.7   # name may differ depending on when you installed it
   conda install -c bioconda -c conda-forge itsxpress bbmap
   qiime dev refresh-cache                # picks up the newly installed itsxpress plugin

   pip install qiime2-its                 # from PyPI: https://pypi.org/project/qiime2-its/

   # Test the pipeline:
   qiime2-its -h
   ```
   `pip install qiime2-its` must run *inside* the activated QIIME2 environment: the package
   drives that environment's `qiime` command and checks for it at start-up. To upgrade later,
   `pip install -U qiime2-its` in the same environment. To work from a clone instead (for the
   validation suite, or to contribute), see [Development](Development).

Next: [build a classifier](Building-a-Classifier), then see [Pipeline usage](Pipeline-Usage).
