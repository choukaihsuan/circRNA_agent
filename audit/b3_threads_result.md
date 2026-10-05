# B3: does config.threads change results? (subsample test)

Status: result recorded. Decision: `threads` is Tier 2 (changes results), not Tier 1.

## Setup

- Sample: SRR444655 (SRP011042), first 10,000,000 read pairs (fastq-dump -X), scratch only.
- Runs: threads=8 and threads=24, same repo code (HEAD + B2), same config except `threads`.
- Steps: fastp → STAR (paired, mate1, mate2) → CIRIquant. DCC and consensus excluded.
- Each run: `--rerun-triggers mtime`, scratch results directory, scratch repo copy.
- Gate: CIRIquant GTF circRNA records = 445 (≥ 100 required), so the comparison is not against an empty result.

## Endpoint: CIRIquant GTF (feeds consensus_filter)

- Circ IDs: 445 in each run, all shared, same coordinates and strand.
- 9 of 445 circRNAs differ:
  - `bsj`: 1 locus (12 vs 11)
  - `fsj`: 2 loci (10 vs 11; 151268 vs 151267)
  - `junc_ratio`: 1 locus (derived from fsj)
  - `score` (abundance): 8 loci, differences of 0.0001 to 0.0004
- Conclusion: the endpoint is NOT identical across thread counts → Tier 2.

## Diagnostics (not conclusions)

- fastp trimmed FASTQ (decompressed content): identical, both mates.
- STAR Chimeric.out.junction: raw differs; sorted identical → order only.
- STAR alignments (SAM body, no @ header): unsorted differs; sorted identical → order only.
- So the content of the fastp and STAR outputs is identical. The GTF difference arises downstream, in CIRIquant (its internal HISAT2/BWA alignments were not diagnosed here).

## Implications

- Published datasets used several thread counts (config.threads: 8, 12, 18, 24 across the 20 projects, see audit/parameter_table). Their circRNA counts may differ slightly from a rerun with another thread count.
- Appendix wording: "In this subsample, CIRIquant output differed for 9 of 445 circRNAs (1–2 read counts or abundance differences in the 4th decimal) when threads changed from 8 to 24; the thread count is therefore recorded per dataset and exposed as a Tier 2 parameter."

## Limits

- One 10M-pair subsample of one sample; not a proof for full-depth data.
- DCC and consensus were not run, so the consensus-level effect is not measured.
- The 2M-pair attempt yielded 78 circRNAs (below the 100 gate) and DCC failed on zero passing circRNAs (see audit/dcc_zero_output_pipeline_failure.md); that run is not used.
