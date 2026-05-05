# concordR: Concordance-Aware QC for Proteomic Interpretation

**concordR** provides an atlas of RNA–protein concordance across 11,154 genes, 24 human tissues, and 53 cell types. One job: given a gene list from a proteomics experiment, tell you what you can and cannot conclude from it.

## Why

Standard proteomic interpretation assumes RNA predicts protein. For ~34% of genes it doesn't, they are consistently suppressed at the protein level despite high RNA. Cell-type enrichment via scRNA, biomarker nomination from intracellular proteins in plasma, and pathway enrichment driven by discordant genes all produce results that do not survive concordance correction.

## Install

```r
remotes::install_github("Art83/concordR")
```

## Functions

| Function | Question it answers |
|---|---|
| `query_atlas()` | What is this gene's concordance class, suppression rate, mechanistic tier? |
| `flag_genes()` | Which genes in my list have reliable RNA–protein concordance in this context? |
| `audit_enrichment()` | Do my enrichment results survive concordance correction? |
| `triage()` | Is this gene defensible as a drug target? (concordance + compartment + tissue gates) |
| `plot_gene()` | Show me this gene's tissue-resolved concordance profile. |
| `plot_audit()` | Show me standard vs corrected enrichment side by side. |
| `load_atlas()` | Load the atlas data (auto-cached, called internally). |

## Quick start

```r
library(concordR)

# 1. What do I have?
flag_genes(c("NEFL", "GFAP", "NRGN", "HDAC1", "TP53"), sample_type = "plasma")
#>   gene_symbol                gene_class  flag
#> 1        NEFL  consistently_concordant  reliable
#> 2        GFAP  consistently_concordant  reliable
#> 3        NRGN  consistently_suppressed  unreliable
#> 4       HDAC1  consistently_concordant  unreliable   # concordant but nuclear → implausible in plasma
#> 5        TP53  consistently_suppressed  unreliable

# 2. Do my enrichment results hold up?
audited <- audit_enrichment(
  enrichment_result = my_clusterProfiler_output,
  gene_list         = my_genes,
  universe          = all_measured_genes,
  sample_type       = "csf"
)
audited[, c("Description", "pvalue", "padj_corrected", "status_change")]

# 3. Is my candidate defensible as a drug target?
triage(c("AR", "ESR1", "CD274"), sample_type = "tissue")
#>   gene_symbol                gene_class           verdict
#> 1          AR  consistently_suppressed    RNA-unreliable
#> 2        ESR1  consistently_suppressed    RNA-unreliable
#> 3       CD274  consistently_suppressed    RNA-unreliable
# All three are established drug targets — but RNA-based evidence
# for their abundance is unreliable. Protein-level measurement required.

# Compartment-implausible nominations from a CSF proteomics study:
triage(c("HDAC1", "TP53", "NEFL"), sample_type = "csf")
#>   gene_symbol  compartment_rescue                   verdict
#> 1       HDAC1               FALSE   compartment_implausible
#> 2        TP53               FALSE   compartment_implausible
#> 3        NEFL                TRUE             RNA-supported
# HDAC1 / TP53 are nuclear with no rescue route into CSF — the proteomic
# nomination is suspect regardless of how significant the differential
# abundance was. NEFL is also non-secreted/non-membrane, but as a
# projection protein under-detected by HPA microscopy it gets a damage-
# release rescue (axonal injury is a well-characterised CSF route).

# Tissue-of-origin gate: catches "this brain protein" claims that don't
# survive the per-tissue protein-level data. The atlas resolves brain
# into four regions; pass one explicitly.
triage(c("NEFL", "GFAP", "SPANXA1"),
       sample_type    = "csf",
       context_tissue = "Cerebellum")
#>   gene_symbol  context_tissue_class  context_tissue_detected      verdict
#> 1        NEFL  concordant_in_tissue                     TRUE  RNA-supported
#> 2        GFAP  concordant_in_tissue                     TRUE  RNA-supported
#> 3     SPANXA1        low_expression                    FALSE  tissue_origin_mismatch
# SPANXA1 fails: testis-restricted, not detected in cerebellum.
# Whatever cell-type enrichment said it was "neuronal" was wrong about
# the protein level, and so a target nomination on that basis is hollow.
```

## Companion to

> Shvetcov A (2026). *An integrated atlas of RNA–protein concordance
> reveals the determinants of translational fate across human tissues.*

The atlas paper provides the foundation. This package provides the tool.

## See also

[**xEnrich**](https://github.com/Art83/xEnrich) — information-theoretic enrichment with CMI redundancy selection and Wilson confidence intervals. Independent package, complementary scope.

## License

MIT © Artur Shvetcov
