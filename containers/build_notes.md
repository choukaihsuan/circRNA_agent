# Container rebuild provenance (2026-10-03)

## Two images, two different jobs (added 2026-10-31)

This file originally described one rebuild effort. It has since split into
two images with different jobs and different obligations, and conflating
them would commit the platform to never upgrading past Table 3's versions:

| | Job | Built from | Must match Table 3? |
|---|---|---|---|
| **Image A** | carries the paper's reproducibility claim | `envs/ciriquant_published.explicit.txt` **+** `containers/pip-only-requirements.txt` (two files -- see below; `envs/ciriquant_published.from-history.yml` is an intent record alongside these, not build input) | Yes — that is its entire purpose |
| **Image B** (`1.1.0-candidate` below, and whatever follows it) | the platform's own, forward-moving container | the Dockerfile in this repo | No — it iterates normally; `find_circ` and version drift from Table 3 are expected here and not a defect |

Image A is the only one that gets frozen and bound to a Zenodo DOI once it
exists. Image B is disclosed as "not the paper's environment -- see Image A
for that" rather than being held to Table 3's pins.

**Image A's recipe is two files, not one -- this was missed when
`envs/ciriquant_published.explicit.txt` was first written and called "the
real reproducibility asset."** `conda list --explicit` (what that file is)
never captures pip-installed packages at all. 14 of the packages that
actually ran on the server -- Flask among them -- only exist via `pip
install`, and are not in that file in any form. Proof this isn't
theoretical: the canary build below only got 27/27 on the Table 3 check
*after* `containers/pip-only-requirements.txt` was also installed into it;
with only the conda explicit spec, Flask was absent and the check would
have read 26/27. **Both files together are Image A's recipe.** Re-verified
2026-10-31: every one of `pip-only-requirements.txt`'s 14 pinned versions
still matches the server's current `pip list --format=freeze` exactly (no
new file was needed -- there was no drift to correct).

## Image A: canary test of `envs/ciriquant_published.explicit.txt` (2026-10-31)

Before building Image A, the explicit spec it depends on was tested for
whether it still resolves -- every line in a `conda list --explicit` capture
names one exact package build from a specific channel, and bioconda/
conda-forge do prune old builds over time. This was the actual risk being
tested, not a formality:

```
conda create --name ciriquant_canary_A --file envs/ciriquant_published.explicit.txt
```

Result: **all 450 packages resolved and installed with zero errors** (no
404s, no "package not found," checked by grepping the full build log for
error/404/not-found/fail -- none). `containers/pip-only-requirements.txt`'s
14 pip packages were then installed into the same canary env
(`pip install --no-deps -r containers/pip-only-requirements.txt`) and also
resolved with zero errors.

Full version check against Table 3 (same method as the earlier Docker-image
and server-env checks in `audit/`): **27/27 rows match**, including Flask
2.2.5 (only present after the pip-only-requirements step -- `conda list
--explicit` does not capture pip-installed packages at all, which is why
that file is paired with the pip-only requirements rather than being
sufficient on its own) and Bowtie2/CIRCexplorer2 (both conda packages,
came in with the rest of the explicit spec with no extra step needed).
BiocManager is absent (same as every other check this project has done --
not installed, not inferred). The canary env was removed after this check
(`conda env remove -n ciriquant_canary_A`); nothing from it was kept.

**Verification level, stated explicitly so it is not read as more than it
is**: this confirms the package *versions* this file lists still resolve and
match Table 3 today. It does **not** run Snakemake, does **not** run any
pipeline stage, and does **not** claim that output from this environment
would be bit-identical to what originally produced the 19 published
datasets. "Same versions" and "same output" are different claims; only the
first one has been checked here. Find_circ (see `patches/find_circ_py3.patch`)
is not yet part of this canary test -- it is not a conda package and is
handled separately (see that patch file's own header for why it is a patch,
not a vendored copy, and for the license reason).

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
- [x] find_circ is not in this (Image B) container, and per the
      two-image split above (added 2026-10-31) it does not need to be --
      Image B does not carry the reproducibility claim Table 3 describes.
      find_circ is handled for **Image A** instead, via
      `patches/find_circ_py3.patch` (fetch upstream commit `8655dca5`,
      verify its sha256, apply the patch at build time) -- never vendored
      into this repo, since find_circ is GPLv3 and this repo is MIT.

## Image A + find_circ: wired in and re-checked (2026-10-31)

`containers/build_image_a_find_circ.sh` fetches upstream find_circ.py at
commit `8655dca5`, verifies its sha256, applies
`patches/find_circ_py3.patch`, and verifies the RESULT's sha256 against the
exact file running on the server -- both checks are hard failures (non-zero
exit), not warnings. Tested three ways: (1) the real upstream file passes
both checks and the script's output is byte-identical to the server's
find_circ.py; (2) a deliberately wrong "upstream" sha256 hard-fails before
ever reaching the patch step; (3) a deliberately wrong "expected patched"
sha256 hard-fails after patching and deletes the bad output rather than
leaving it on disk.

**A real bug was caught by check (1) on the first actual run, not a
contrived test**: the first version of `find_circ_py3.patch` (committed
earlier today) was generated against the server's local
`find_circ_orig.py`, trusting its filename to mean "genuine unmodified
upstream." It wasn't -- diffing a fresh fetch of the real commit `8655dca5`
against that "orig" file turned up 8 further differences (Python 2 `print`
statements, `izip_longest`, three `dict.values()/.items()` calls) that
predate this project's own patch entirely; `find_circ_orig.py` had already
been through a separate, unrecorded Python 3 pass before anyone here used
it as a baseline. The build script's own upstream-sha256 check caught this
immediately (hard-failed on the real mismatch) rather than silently
building an image around a wrong assumption. `patches/find_circ_py3.patch`
has been corrected: it is now a diff against the genuine upstream blob
(273 lines, not 222), the upstream sha256 in its header is the real one
(`b67a0439...`), and the resulting patched file is unchanged (`28a3cb9e...`
-- still bit-identical to the server). See that patch file's own
"CORRECTION" note for the full account.

Full Table 3 re-check with find_circ wired in (conda create from
`envs/ciriquant_published.explicit.txt`, pip install from
`containers/pip-only-requirements.txt`, then `build_image_a_find_circ.sh`
into the same env's `bin/`): **28/28 rows now match**, with find_circ's row
now asserted at the sha256 level -- not just "commit 8655dca5 + has a
patch that looks plausible," but "produces the exact bytes that generated
the 19 published datasets." Canary env removed after the check.

## find_circ license disclosure (GPLv3)

find_circ (Memczak et al., github.com/marvin-jens/find_circ) is licensed
GPLv3. circDEX itself is MIT. Building Image A fetches find_circ's source
and applies a patch to it at build time (see above) -- the resulting image
therefore **contains GPLv3-licensed code**, which carries disclosure and
distribution obligations under GPLv3 for that component specifically. This
does not make circDEX itself GPL, and MIT code elsewhere in the image is
unaffected -- but it is also not adequately covered by circDEX's own MIT
notice alone.

**Disclosure** (to be carried in both this file and the image's own
`LABEL`, for any image that bundles find_circ):

> This image contains find_circ (GPLv3), fetched from upstream commit
> 8655dca5 and patched for Python 3 compatibility. circDEX itself is
> licensed MIT; the find_circ component is distributed under its own
> GPLv3 license, and the modifications applied to it are published as
> `patches/find_circ_py3.patch` in this repository as the corresponding
> source for that component.

Currently only Image A is expected to bundle find_circ (Table 3 lists it;
Image B does not need to, per the two-image split above). If Image B, or
any future image, also bundles find_circ, the same disclosure applies to
it. Not yet added to an actual Docker `LABEL` -- Image A has no Dockerfile
yet, only the conda-env equivalent tested above; this is noted here so it
isn't forgotten when one is written. If Image A is ever pushed publicly
and bound to a Zenodo DOI, this wording is worth a second read by someone
with more licensing background before that happens -- not a blocker, but
worth asking.
