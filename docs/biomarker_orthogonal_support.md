# Orthogonal support for top biomarker candidates

Every dimension of the composite biomarker score (`rank_biomarkers.py`) is computed
from a single pipeline run, so a high rank is not *independently* corroborated. This
document collects **orthogonal evidence** from external, independently-generated data
for the top-ranked candidates:

1. **Literature** — is the specific circRNA (circBase id) already reported as a
   biomarker in the matching or a related cancer?
2. **Host-gene clinical/survival** — does the circRNA's host gene have an established
   TCGA/cohort survival association (the circRNA itself is often novel, but its host
   locus's clinical relevance is orthogonal support)?

All citations below come from grounded literature searches (July 2026); each row notes
whether the **direction** (up/down in tumor) is consistent between our result and the
external evidence. Direction conflicts are stated plainly, not hidden.

> Method note for the paper: this is corroborative, not confirmatory. "Reported in
> literature" and "host-gene survival association" reduce the chance a top candidate is
> a pipeline artefact; they do not substitute for experimental validation (qRT-PCR /
> RNase R) of the specific circRNA.

---

## Tier 1 — circRNA isoform itself is an established biomarker (strongest)

### hsa_circ_0004771 (circ_NRIP1) — NRIP1, chr21
- **Our result:** GSE248612 gastric cancer, **up-regulated** (log2FC = +5.58, Type_I,
  score 0.72; a top-3 candidate).
- **Literature:** circ_NRIP1 (= hsa_circ_0004771, from NRIP1 exon 2–3) is a
  well-characterised **oncogenic** circRNA:
  - Gastric cancer — plasma **diagnostic & dynamic-monitoring biomarker**
    ([PMC7549879](https://pmc.ncbi.nlm.nih.gov/articles/PMC7549879/)).
  - Breast cancer — up-regulated, proposed **diagnostic biomarker** (with CircCSPP1 &
    CircSMAD2) ([PMC12396114](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC12396114/)).
  - ESCC — up-regulated, correlates with tumour burden & poor prognosis
    ([PMID 32050790](https://pubmed.ncbi.nlm.nih.gov/32050790/);
    [Cancer Cell Int PMC8101145](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8101145/)).
  - Cervical cancer — oncogenic, promotes invasion
    ([Cell Death Dis](https://www.nature.com/articles/s41419-020-2607-9)).
- **Verdict:** ✅ **Direct validation** — same circRNA, matching cancer type (gastric),
  matching direction (up / oncogenic).

---

## Tier 2 — same host gene yields a reported oncogenic circRNA in a matching cancer

### PTK2 (Protein Tyrosine Kinase 2 / FAK), chr8
- **Our results:** GSE121842 CRC — hsa_circ_0002483, **up** (log2FC = +9.28, Type_II,
  top-1, score 0.85); GSE133998 breast — hsa_circ_0003221 (PTK2), a top candidate.
- **Literature:** circPTK2 (**hsa_circ_0005273**, a different isoform of the same host
  gene) is an established **oncogenic** circRNA:
  - Colorectal cancer — elevated, **poor prognosis**, promotes proliferation &
    metastasis ([PMC6977296](https://pmc.ncbi.nlm.nih.gov/articles/PMC6977296/)).
  - Breast cancer — circPTK2/0005273 regulates YAP1–Hippo signalling
    ([PMC7802350](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7802350/)).
- **Verdict:** ✅ **Host-gene circRNA support** — the PTK2 locus produces oncogenic
  circRNAs in CRC (matching cancer, matching up-direction) and breast; our circRNAs are
  different isoforms of the same validated locus.

---

## Tier 3 — circRNA appears novel, but host gene has a strong survival association

### SLC38A1 (amino-acid transporter SNAT1), chr12  ← strongest recurrent candidate
- **Our results:** **up-regulated in two independent breast cohorts** —
  GSE113230 (log2FC +6.94, top-1) and SRP156355 (log2FC +5.35, top-3); circRNA
  hsa_circ_0000397.
- **Host-gene survival:** High SLC38A1 expression → **poor prognosis in breast cancer**,
  higher in tumour vs. adjacent normal, an independent prognostic factor
  ([PMC11656554](https://pmc.ncbi.nlm.nih.gov/articles/PMC11656554/)); also unfavourable
  prognosis / pro-metastatic in lung adenocarcinoma
  ([BMC Cancer s15329-9](https://link.springer.com/article/10.1186/s12885-025-15329-9)).
- **circRNA literature:** hsa_circ_0000397 itself is **not** previously reported in
  breast cancer → potentially **novel**.
- **Verdict:** ✅ Host-gene up-regulation and poor-prognosis association match our
  up-regulated circRNA, in the same cancer, reproduced across two cohorts. A strong,
  paper-worthy candidate (novel circRNA on a clinically-validated locus).

### BACH1 (BTB domain transcription factor), chr21
- **Our result:** GSE77509 HCC, **up-regulated** (hsa_circ_0001181, log2FC +5.49, Type_I,
  top-2).
- **Host-gene survival:** BACH1 is up-regulated in HCC, associated with **poor overall
  survival** and high recurrence, an **independent predictor**; drives HCC growth and
  metastasis via IGF1R and **PTK2**
  ([Theranostics PMC8771560](https://pmc.ncbi.nlm.nih.gov/articles/PMC8771560/)). Also a
  survival predictor in early-stage lung adenocarcinoma
  ([PMC8758312](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8758312/)).
- **Verdict:** ✅ Host-gene direction (up in HCC, poor prognosis) matches our
  up-regulated circRNA. Note the BACH1→PTK2 axis dovetails with the recurrent PTK2
  finding above.

---

## Direction conflicts / cautions (reported honestly)

### XPO1 (Exportin-1 / CRM1), chr2 — hsa_circ_0001017
- **Our results:** GSE248612 gastric **up** (log2FC +6.61, top-1) **but** GSE130078 ESCC
  **down** (log2FC −1.71) — our own two datasets disagree on direction.
- **Host gene:** XPO1 is a validated pan-cancer poor-prognosis marker and drug target
  (selinexor/CRM1) ([pan-cancer PMC8797940](https://pmc.ncbi.nlm.nih.gov/articles/PMC8797940/);
  [GI cancers PMC11680630](https://pmc.ncbi.nlm.nih.gov/articles/PMC11680630/)); in gastric
  its prognostic effect is stage/sex-dependent (nuanced).
- **circRNA literature:** hsa_circ_0001017 has been reported as **tumour-suppressive /
  chemo-sensitising** in gastric cancer (miR-543/PHLPP2 axis,
  [Biochem Genet](https://link.springer.com/article/10.1007/s10528-021-10110-6)) — the
  **opposite** direction to our gastric up-regulation.
- **Verdict:** ⚠️ **Conflicting** at the circRNA level and internally inconsistent across
  our cohorts. Host-gene XPO1 clinical relevance is solid, but this specific circRNA
  should not be presented as a confident up-regulated gastric biomarker without wet-lab
  confirmation.

---

## Convergence note (worth a sentence in Discussion)

A published breast-cancer circRNA biomarker panel highlights **CircCSPP1 + CircNRIP1 +
CircSMAD2** ([PMC12396114](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC12396114/)). Our
pipeline independently surfaced circRNAs from **all three host loci** as top candidates
(NRIP1 → gastric GSE248612; CSPP1 → hsa_circ_0003388, LUAD GSE148036 top-1;
SMAD2 → hsa_circ_0000847, prostate GSE221107 top-3). Independent recovery of the same
biomarker-producing loci across unrelated cancer types corroborates that these loci
generate disease-relevant circRNAs, and that the ranking is capturing real signal.

---

## Cross-dataset recurrence (from `cross_dataset_biomarkers.py`, real data)

Top candidates recurring across independent cohorts (orthogonal by construction):

| circBase id | host gene | cohorts (direction) | note |
|---|---|---|---|
| hsa_circ_0000397 | SLC38A1 | GSE113230 breast (↑), SRP156355 breast (↑) | direction-consistent; host gene poor-prognosis in breast |
| hsa_circ_0000471 | N4BP2L2 | GSE133998 breast (↑), SRP156355 breast (↓), PRJNA553289 SCLC (↓) | recurs in 3 cohorts; direction mixed → investigate |
| hsa_circ_0031584 | ARHGAP5 | GSE77509 HCC (↑), PRJNA553289 SCLC (↑) | direction-consistent |
| hsa_circ_0001359 | PHC3 | GSE77509 HCC (↓), GSE130078 ESCC (↑) | mixed direction |
| hsa_circ_0001017 | XPO1 | GSE248612 gastric (↑), GSE130078 ESCC (↓) | mixed; see caution above |
| — (PTK2 locus) | PTK2 | GSE121842 CRC (↑), GSE133998 breast | different isoforms, same oncogenic locus |

`hsa_circ_0000471 / N4BP2L2` recurs in the most cohorts (3) but with mixed direction —
a good target for focused follow-up; it is comparatively understudied in the circRNA
literature.

---

---

## Focused follow-up: hsa_circ_0000471 / N4BP2L2 (chr13:33091994|33101669)

The most-recurrent candidate in the panel, and an instructive case for why recurrence
alone is not enough — direction matters.

**Our data (top-ranked in three independent cohorts, all Type_I):**

| cohort | cancer | log2FC | direction | Type | score | rank |
|---|---|---|---|---|---|---|
| GSE133998 | breast | +6.65 | **up** | Type_I | 0.728 | 1 |
| SRP156355 | breast (early IDC) | −6.89 | **down** | Type_I | 0.741 | 1 |
| PRJNA553289 | SCLC | −8.71 | **down** | Type_I | 0.750 | 2 |

**Host gene N4BP2L2 (a.k.a. PFAAP5):** a nuclear transcriptional regulator,
phosphorylated on Ser199 in response to DNA damage (ATM/ATR), interacts with the Gfi1
repressor and neutrophil elastase, and functions in hematopoietic stem-cell
differentiation/proliferation; associated with congenital neutropenia
([PMC2725743](https://pmc.ncbi.nlm.nih.gov/articles/PMC2725743/);
[GeneCards](https://www.genecards.org/cgi-bin/carddisp.pl?gene=N4BP2L2);
[OMIM 615788](https://omim.org/entry/615788)). It is **not** an established cancer gene —
no TCGA/cohort survival association was found; the DNA-damage-response link is the only
cancer-adjacent hook, and it is speculative.

**circRNA literature:** hsa_circ_0000471 itself has **no** published cancer report —
genuinely understudied / novel.

**Assessment (honest):**
- ✅ Strongest *recurrence* signal in the panel — a top-1/2 candidate in 3 independent
  cohorts spanning 2 cancer types. Unlikely to be pure noise.
- ⚠️ **Direction is unstable**: the two independent breast cohorts *disagree*
  (GSE133998 up vs SRP156355 down). Because `edgeR_ciriquant` tests the BSJ/FSJ ratio,
  the sign reflects a change in circularisation efficiency relative to the linear host
  transcript, not absolute abundance — but a sign flip between two same-cancer cohorts is
  still a red flag (subtype/stage difference, paired-vs-unpaired design, or a few outlier
  samples).
- ⚠️ No prior biological grounding for either the circRNA or the host gene in cancer.

**Recommended before treating it as a biomarker:**
1. Inspect the raw BSJ and FSJ matrices per cohort — is the sign driven by BSJ change or
   by FSJ (linear) change? Is it a few outlier samples?
2. Compare cohort composition (GSE133998 general breast / unpaired vs SRP156355 early IDC
   / paired) to see whether the flip tracks a clinical variable.
3. If it survives (1)–(2), it is a strong *novel* target for wet-lab validation
   (qRT-PCR + RNase R); if the flip is an artefact, drop it. Do **not** present it as a
   confident biomarker on recurrence alone.

This candidate is the clearest example of why `cross_dataset_biomarkers.py` reports a
`dir_consistent` flag: recurrence with consistent direction (e.g. SLC38A1) is support;
recurrence with a sign flip is a lead that needs investigation, not a result.

---

*Generated as orthogonal validation for the biomarker ranking. Regenerate the recurrence
table with `scripts/cross_dataset_biomarkers.py`; refresh citations before submission.*
