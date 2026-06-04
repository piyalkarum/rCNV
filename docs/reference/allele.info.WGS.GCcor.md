# Allele Information with GC-Bias Correction for WGS Data

Detects per-SNP duplication (and optionally deletion) signals in
whole-genome sequencing (WGS) data using a joint genotype +
sequencing-depth likelihood-ratio framework, with GC-bias correction
applied at the single-SNP level. Sequencing depth is modelled separately
within each 1-% GC bin using SNP-centred windows, so the GC and depth
signals share the same spatial footprint (avoiding the geometry mismatch
that arises when GC is estimated from sliding windows but depth from SNP
positions).

## Usage

``` r
allele.info.WGS.GCcor(
  vcf = NULL,
  fis = NULL,
  parallel = FALSE,
  ref = NULL,
  gc_win_size = 1000,
  numCores = NULL,
  nb_para = NULL,
  sig = 0.95,
  deletion = FALSE,
  perm_test = TRUE,
  chroms = NULL,
  bed = NULL,
  ...
)
```

## Arguments

- vcf:

  VCF input. Either (a) a file path (character) that will be read via
  [`readVCF()`](readVCF.md), (b) a list returned by
  [`readVCF()`](readVCF.md) with a `$vcf` element, or (c) a `data.frame`
  whose first column holds chromosome names.

- fis:

  Numeric scalar, global inbreeding coefficient. If `NULL` (default), it
  is estimated internally with [`h.zygosity()`](h.zygosity.md) and
  reported via message.

- parallel:

  Logical. If `TRUE`, parallelise the per-sample and permutation loops
  with a `parallel` cluster.

- ref:

  Character path to the reference FASTA (required). Must be indexable by
  [`Rsamtools::indexFa()`](https://rdrr.io/pkg/Rsamtools/man/scanFa.html);
  the index will be created on the fly if missing.

- gc_win_size:

  Integer window size in bp used to compute GC content around each SNP
  (default `1000`). The window is centred on `POS`:
  `[POS - floor((win-1)/2), POS + floor(win/2)]`.

- numCores:

  Number of cores to use when `parallel = TRUE`. Defaults to
  `parallel::detectCores() - 1`.

- nb_para:

  Optional pre-computed negative-binomial parameter table (one row per
  sample, with columns `mean`, `size`, `mean2`, `size2`). If `NULL`, the
  table is fit internally.

- sig:

  Numeric in `(0, 1)`, quantile threshold for the permutation-based
  significance cut-off (default `0.95`).

- deletion:

  Logical. If `TRUE`, also test for deletions (N = 1 vs. N = 2) in
  addition to duplications (N = 4 vs. N = 2).

- perm_test:

  Logical. If `TRUE` (default), estimate per-Nhet significance
  thresholds by permutation; if `FALSE`, the function returns raw
  likelihood ratios only (no `*.sig` or `*.stat` columns).

- chroms:

  Optional character vector of chromosome names to process (e.g.
  `c("chr1", "chr2")`). Mutually exclusive with `bed`.

- bed:

  Optional path to a BED file. Only SNPs overlapping at least one BED
  interval are processed. Mutually exclusive with `chroms`.

- ...:

  Currently unused; reserved for future arguments.

## Value

A single `data.frame` containing one row per processed SNP (merged
across all chromosomes / regions). Columns include:

- **`CHROM, POS, win_start, win_end, win_width, win_clipped`**: SNP
  coordinates and the GC window actually used (clipped at chromosome
  boundaries).

- **`gc_percent, nonN_fraction, gc_idx`**: GC content of the window,
  fraction of non-`N` bases, and the 1-% GC bin.

- **`mdep, Nhet, Nmiss, p_alt, delta`**: Mean depth across samples,
  heterozygous-genotype count, missing-genotype count, alt-allele
  frequency, and the deviation from expected Nhet given \\F\_{is}\\.

- **`lhr.dup, dup.sig, dup.stat`**: Duplication likelihood ratio,
  permutation threshold for the SNP's Nhet, and classification
  (`"duplicated"` / `"non-duplicated"`).

- **`lhr.del, del.sig, del.stat`**: (only when `deletion = TRUE`)
  deletion analogues of the above.

When `perm_test = FALSE`, the `*.sig` and `*.stat` columns are omitted
and only the raw `lhr.*` columns are returned.

## Pipeline overview

1.  Optional region filter (`chroms` or `bed`).

2.  GC content estimated for a `gc_win_size`-bp window centred on each
    retained SNP, then binned at 1-\\

3.  Global inbreeding coefficient (\\F\_{is}\\), negative-binomial depth
    parameters, and beta-binomial allelic-ratio parameters are fit on
    the full filtered set so per-chromosome likelihoods remain directly
    comparable.

4.  Per-chromosome loop: depth mean/SD is re-fit within each 1-\\ bin,
    depth-based and allelic-ratio likelihoods are computed for copy
    numbers N = 2, N = 4 (and optionally N = 1).

5.  Likelihood ratios for duplication (and deletion) are computed per
    SNP; if `perm_test = TRUE`, significance thresholds are estimated
    from a permutation pool drawn across all chromosomes.

6.  All per-chromosome results are merged into one `data.frame`.

## Region-filter behaviour

`chroms` and `bed` cannot be supplied together. With no filter, every
SNP in the VCF is processed. With `chroms`, SNPs whose `CHROM` is in the
vector are kept (missing chromosomes generate a warning). With `bed`,
SNPs overlapping any interval are kept. In all cases, the per-chromosome
loop iterates over the chromosomes that actually appear in the retained
SNPs, and outputs are merged before return.

## Statistical caveat

The per-chromosome GC-depth model uses only SNPs on that chromosome. If
a chromosome carries fewer than ~500 SNPs, the 1-\\ and the fitted depth
mean/SD may be unstable; the function emits a warning in that case but
still returns the result.

## Memory

The dominant memory term is the per-chromosome likelihood matrix
(`nSNP_chr` rows by `4 * nsample` or `6 * nsample` columns). It is
allocated and freed inside each iteration of the chromosome loop, so
peak memory scales with the largest chromosome, not the full genome.

## See also

[`allele.info.WGS3`](allele.info.WGS3.md) (non-GC-corrected variant
using the negative-binomial depth model directly), `run_gc_estimation`,
`parse_region_filter` (internal helpers).

## Examples

``` r
if (FALSE) { # \dontrun{
# Whole genome with default 1 kb GC window
res <- allele.info.WGS.GCcor(
  vcf = "sample.vcf.gz",
  ref = "genome.fa"
)

# Restrict to a few chromosomes, also test deletions
res <- allele.info.WGS.GCcor(
  vcf      = "sample.vcf.gz",
  ref      = "genome.fa",
  chroms   = c("chr1", "chr2", "chr3"),
  deletion = TRUE
)

# Restrict to a BED file of candidate regions, in parallel
res <- allele.info.WGS.GCcor(
  vcf       = "sample.vcf.gz",
  ref       = "genome.fa",
  bed       = "candidate_regions.bed",
  parallel  = TRUE,
  numCores  = 8,
  perm_test = TRUE
)
} # }
```
