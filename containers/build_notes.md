# Container rebuild provenance (2026-10-03)

## Why this exists

The paper's Table 1 lists exact verified software versions, recorded directly
from the HPC analysis server's `ciriquant` conda environment (the environment
that actually produced all 19 published datasets). `choukaihsuan/circrna-pipeline:1.0.1`
predates that verification and has several version mismatches against Table 1
(see `audit/code_freeze.md` and `CLAUDE.md`'s "資料集特性稽核校正" /
Table 2 reproducibility sections). This work rebuilds the image to match the
verified versions, without touching the existing `1.0.1` tag.

## Why `envs/circrna.yaml` could not be used to build this image

`envs/circrna.yaml` was edited to pin Table 1's exact versions (python=3.7.12,
star=2.7.11a, samtools=1.18, r-base=4.2.1, ciriquant=1.1.2, etc.). Running
`mamba env create -f envs/circrna.yaml` against today's live conda-forge/
bioconda package index fails with `Could not solve for environment specs`:
samtools/star now require htslib builds that link against libcurl>=8.x,
which in turn requires openssl>=3.0 and krb5>=1.22 -- but python=3.7.12
requires an old-style openssl (1.1.1-series) build. No combination of
currently-indexed package builds satisfies both simultaneously. This is a
real infeasibility in today's package index, not a mistake in the pins.

## The approach actually used: explicit spec from the live server

Rather than ask a solver to reconstruct a combination it can no longer find,
this build reproduces the **exact package set already running** on the
verified server -- which is known to be mutually compatible because it is
literally what executed the 19 published analyses.

1. `containers/spec-file.linux-64.txt` -- captured via, on the server:
   ```
   export PATH=/home/choukaihsuan/miniconda3/envs/ciriquant/bin:$PATH
   conda list --explicit -n ciriquant
   ```
   450 packages (the file is 455 lines: 5 header lines plus 450 package URLs), all from `conda.anaconda.org` (conda-forge/bioconda) or
   `repo.anaconda.com/pkgs/main` (pandas is a `defaults`/`pkgs/main` build).
   **All 450 download URLs were checked with a HEAD request on 2026-10-03
   and returned HTTP 200 -- none have been pruned from their channel.**
   `conda create --name ciriquant --file spec-file.linux-64.txt` installs
   this set directly, with no solving involved (the versions are already
   fixed in the file).

2. `containers/pip-freeze.txt` -- captured via, on the server:
   ```
   /home/choukaihsuan/miniconda3/envs/ciriquant/bin/pip freeze
   ```
   165 lines. Most of these are conda-managed packages that also show up
   in `pip freeze` (conda-forge/bioconda builds register pip-visible
   dist-info); those lines carry an `@ file:///home/conda/feedstock_root/...`
   or `@ file:///opt/conda/conda-bld/...` build-artifact path and are
   **not** re-installed via pip -- those paths only existed on conda-forge's
   own CI and don't exist in this build.

3. `containers/pip-only-requirements.txt` -- the 14 lines from
   `pip-freeze.txt` that are **not** matched by name in
   `spec-file.linux-64.txt` (cross-checked both by name-matching and, for
   every ambiguous case, by directly running `import X; X.__version__` on
   the server to confirm which version is actually importable at runtime):
   `Flask==2.2.5`, `Flask-WTF==1.1.1`, `itsdangerous==2.1.2`,
   `kaleido==0.2.1`, `lxml==5.4.0`, `Markdown==3.3.7`, `plotly==5.18.0`,
   `pyliftover==0.4.1`, `pysradb==2.2.2`, `requests-ftp==0.3.1`,
   `tenacity==8.2.3`, `tqdm==4.68.2`, `Werkzeug==2.2.3`, `WTForms==3.0.1`.
   Installed with `pip install --no-deps` so pip cannot silently pull in
   different dependency versions than what's already pinned by step 1.

   **Important caveat found during this cross-check**: `pip freeze`'s
   reported version is *stale/unreliable* for several conda-managed
   packages where dist-info from an earlier install evidently lingered
   after a later conda-managed upgrade overwrote the actual module files.
   Confirmed via live `import X; X.__version__` on the server:

   | Package | pip freeze says | actually importable (ground truth) | source of truth |
   |---|---|---|---|
   | scipy | 1.2.2 | **1.7.3** | conda-forge, confirmed by runtime import |
   | scikit-learn | 0.20.3 | **1.0.2** | conda-forge, confirmed by runtime import |
   | pysam | 0.15.2 | **0.21.0** | bioconda, confirmed by runtime import |
   | PyYAML | 5.4 | **6.0** | conda-forge, confirmed by runtime import |
   | numexpr | 2.6.9 | **2.8.3** | conda-forge, confirmed by runtime import |

   These five are conda-managed (present in `spec-file.linux-64.txt`) and
   were **not** added to `pip-only-requirements.txt` -- the conda explicit
   spec already installs the correct version. The stale pip freeze numbers
   are left as-is in `containers/pip-freeze.txt` (an unedited raw capture)
   since editing a captured log to "fix" it would defeat its purpose as a
   provenance record.

   For the 14 genuinely pip-only packages, pip freeze's reported version
   **did** match live `import` output wherever checked (Flask, Markdown,
   plotly, pysradb, pyliftover all confirmed) -- the staleness problem is
   specific to packages that were later re-installed via conda, not to
   pip-only packages.

## `envs/circrna.yaml` vs the files under `containers/`

`envs/circrna.yaml` remains in the repo as a readable summary of intent (see
its header comment). It is **not** used by the Dockerfile. The files in this
directory (`spec-file.linux-64.txt`, `pip-only-requirements.txt`,
`pip-freeze.txt`) are the actual build input and the more precise
environment record -- citable directly for the paper's §2.14/§2.15.

## CIRIquant version: three-way check (2026-10-03)

| Source | Reports | Note |
|---|---|---|
| `conda list -n ciriquant` on server | `1.1.2 pypi_0` | conda registry entry is stale |
| `pip show CIRIquant` on server | `1.1.3` | pip metadata |
| `CIRIquant --version` on server and in container | `CIRIquant 1.1.3` | runtime binary |
| bioconda package in spec file | `ciriquant-1.1.3-pyhdfd78af_0` | what the build installs |

The authoritative version is **1.1.3**. Table 1's "1.1.2" and the earlier
CLAUDE.md edit (c7e0559) relied only on the stale conda registry entry.

## Post-build verification (candidate image `choukaihsuan/circrna-pipeline:1.1.0-candidate`, local, 2026-10-03)

Done. Nothing has been pushed or tagged as 1.1.0 yet.

- [x] Versions read from inside the container. Conda package versions match the
      spec file. Binary self-reports were checked for CIRIquant (1.1.3), DCC
      (0.5.0), STAR (2.7.11a), HISAT2 (2.2.1), BWA (0.7.18), samtools (1.18),
      Snakemake (7.24.0), FastQC (0.12.1), Python (3.7.12), and R (4.2.1).
- [x] Compared against Table 1. All rows match. CIRIquant was the only
      mismatch: Table 1 said 1.1.2, the real version is 1.1.3. Per the
      author's decision, Table 1 was corrected to 1.1.3. The container was not
      adjusted to fit the old value.
- [x] Plotly interactivity. Under the container's Python 3.7.12, numpy 1.21.6
      and plotly 5.18.0 import cleanly. A test figure renders to HTML with a
      `Plotly.newPlot` call, which is the interactive path. This tests the
      exact `_PLOTLY` import check in `generate_report.py`. It is not a full
      end-to-end report run; there is no synthetic dataset for one.
- [x] Caveat (present on the live server too, so it is not a container
      artifact): `Rscript` warns that edgeR, limma, and DESeq2 were built under
      R 4.2.3, while the installed r-base is 4.2.1.
- [ ] find_circ is not in the container. Disclose this in paper §2.14 (see
      the scope statement in the commit message).
