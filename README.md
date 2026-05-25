# concordR: Concordance-Aware QC for Proteomic Interpretation

**concordR** provides an atlas of RNA–protein concordance across 11,154 genes, 24 human tissues, and 53 cell types. One job: given a gene list from a biofluid proteomics experiment, tell you what you can and cannot conclude from it.

## Why

Standard proteomic interpretation assumes RNA predicts protein, which is often not true. Some genes are either consistently suppressed despite high RNA. Drug-target nomination from intracellular proteins in plasma/csf, or based on transcriptomic cell-type enrichment, and tissue-of-origin assertions that don't survive protein-level tissue data all produce results that do not hold up under scrutiny. concordR makes those checks systematic and reproducible.

## Install

```r
remotes::install_github("Art83/concordR")
```

## Functions

| Function | What it does |
|---|---|
| `load_atlas()` | Load the atlas (auto-cached, called internally). |
| `query_atlas()` | Look up concordance class, suppression rate, mechanistic tier, localisation. |
| `query_tissue()` | Tissue-resolved lookup: per-tissue concordance, protein detection, cross-platform evidence. |
| `triage()` | Drug-target defensibility audit: compartment, tissue-of-origin, concordance, and specificity gates. |
| `plot_gene()` | Tissue-resolved concordance profile for a single gene. |

## Quick start

```r
library(concordR)



# 1. Triage drug-target claims from a CSF study
triage(c("NEFL", "GFAP", "SPANXA1"),
       sample_type    = "csf")


# 2. Narrow to a specific brain region
triage(c("NEFL", "GFAP", "SPANXA1"), sample_type = "csf", context_tissue = "Cerebellum")


# 3. Visualise a gene's tissue profile
plot_gene("TREM2")
```

## Scope

concordR audits inferential claims made from biofluid proteomics: plasma, CSF, urine. It catches the gap between where you measured and what you're claiming. The four gates in `triage()` are:

1. **Compartment plausibility** - is the protein's subcellular localisation consistent with the biofluid?
2. **Tissue-of-origin** - is the protein detected at protein level in the tissue implied by the biofluid?
3. **Concordance** - does RNA reliably predict protein for this gene?
4. **Tissue specificity** - is the protein broadly expressed, undermining tissue-of-origin claims?

## Companion to

> Shvetcov A (2026). *[An integrated atlas of RNA–protein concordance
> reveals the determinants of translational fate across human tissues.] preliminary name*

The atlas paper provides the foundation. This package provides the tool.

## License

MIT © Artur Shvetcov
