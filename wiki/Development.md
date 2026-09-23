# Development

```bash
git clone https://github.com/duceppemo/QIIME2_ITS
cd QIIME2_ITS
pip install -e '.[dev]'
pytest -v
```

All unit tests mock external tool invocations (`qiime`, `biom`, `bbduk.sh`, etc.), so they run with
a plain Python interpreter — no QIIME2 installation required. CI runs them on every push (see
`.github/workflows/tests.yml`).

## Releasing

Bump the version in `pyproject.toml`, `src/qiime2_its/_version.py` and `CITATION.cff` (plus its
`date-released`), commit, tag `vX.Y.Z`, push both, and create the GitHub release. Then publish to
PyPI from a clean export of the tag, so no untracked file can leak into the package:
```bash
mkdir /tmp/release && git archive vX.Y.Z | tar -x -C /tmp/release && cd /tmp/release
python -m build && python -m twine check dist/* && python -m twine upload dist/*
```
Finally `scripts/sync_wiki.sh`.

## Project layout

```
src/qiime2_its/
  cli/                  pipeline.py, fastq_rc.py, train_unite.py, train_ncbi.py, train_fasta.py
  fastq_utils.py         fastq listing/parsing/cleaning (pure Python)
  env_checks.py           conda-env/CPU/executable sanity checks
  qiime_wrapper.py         `qiime ...` subprocess command builders
  biom_utils.py            `biom ...` subprocess command builders + taxonomy-in-BIOM merging
  itsxpress_wrapper.py     `qiime itsxpress ...` command builders
  size_filter.py           bbduk.sh command builders
  taxonomy.py              NCBI taxdump parsing + lineage-string construction
  downloader.py            generic download + tar.gz extraction
  metadata_utils.py        QIIME2 metadata TSV parsing, column eligibility
  provenance.py            run-provenance/QA metadata (run_metadata.json)
  report_data.py           parses pipeline output (.qzv/.tsv) into plain data structures
  report.py                assembles report.pdf from report_data.py's output
  timing.py                compact elapsed-time formatting ("1d2h3m4s")
tests/                    unit tests, one file per module above
validation/                real-data end-to-end validation (see Validation-Suite)
```

The wrapper modules (`qiime_wrapper.py`, `itsxpress_wrapper.py`, `size_filter.py`) each build one
external-tool command per function and run it with `subprocess.run(cmd, check=True)`. Their tests
mock `subprocess.run` and assert the constructed argv — this catches *our* bugs, not whether the
installed QIIME2 version still accepts that command; that's what the
[validation suite](Validation-Suite) is for.

## Editing this wiki

The GitHub wiki you're reading is a separate git repo
(`github.com/duceppemo/QIIME2_ITS.wiki.git`) that GitHub maintains alongside this one, but its
source of truth is `wiki/*.md` in *this* repo, not that repo directly — that way wiki changes go
through the same review/diff/history as code changes instead of being editable, unreviewed, straight
from the GitHub web UI. To change a page: edit the file under `wiki/`, commit it here as usual, then
run `scripts/sync_wiki.sh` to push the current contents of `wiki/` live to the GitHub wiki (it clones
the wiki repo to a temp dir, mirrors `wiki/*.md` onto it, and pushes — nothing to set up beforehand).

## Migrating from the pre-0.2 flat scripts

Versions before 0.2 shipped as standalone scripts (`python qiime2_its.py ...`). As of 0.2, this is
an installable package with console-script entry points instead: `qiime2-its`, `qiime2-its-rc`,
`qiime2-its-train-unite`, `qiime2-its-train-ncbi`, and `qiime2-its-train-fasta` (`pip install -e .`
from the repository root, inside your QIIME2 conda environment). Command-line flags are otherwise
unchanged. ITS extraction now goes through the `qiime itsxpress` plugin instead of the old
standalone `itsxpress` CLI — install it with `conda install -c bioconda -c conda-forge itsxpress`
followed by `qiime dev refresh-cache` inside your QIIME2 environment.
