# concordR: QC for biofluid proteomics in NDD


## Motivation
For drug nomination, standard proteomic practice (particularly in neurodegenerative disorders) often relies on a chain of assumptions that are rarely tested at the protein level.

**1.** RNA predicts protein abundance, despite substantial discordance between transcript and protein expression. This leads to drug targets being nominated from transcriptomic evidence: a differentially abundant protein (DAP) in plasma -> corresponding gene is found exprerssed in brain -> the DAP is reported as a drug candidate. This assumption doesn't account for the fact that reported correlation between mRNA and protein, replicated in several studies, is only **ρ ~ 0.4**, i.e. mRNA ≠ protein.

**2.** Any DAP discovered in plasma/csf has a natural path to plasma/csf and has a direct causal role in disease. Detecting a protein in plasma or CSF doesn't establish either. Intracellular proteins found in a biofluid may reflect passive release from damaged cells reacting _to_ disease rather than causing it. Such proteins can be solid biomarkers, yet they contribute little to a understanding causes of the disease and selecting druggable targets. Their utility for drug discovery is at best questionable.

**3.** The final assumption is that any DAP in the plasma or CSF of NDD patients originates from brain. This simply ignores the fact that plasma collects proteins from all tissues, and even a protein reliably detected in brain is rarely brain-specific. For CSF the situation is similar: roughly 80% of its proteins originate from plasma, making it an unreliable proxy for brain tissue.

In short, a protein can be a useful circulating biomarker without being a mechanistically relevant drug target, and transcriptomic enrichment of a gene in a brain region does not by itself establish that the corresponding protein is present, brain-derived, or disease-relevant.

**concordR** makes these checks systematic and reproducible. Given a set of differentially abundant proteins from plasma or CSF, it evaluates whether the corresponding proteins are actually detected in the relevant tissues, how specifically they are expressed, whether their subcellular localisation is compatible with their presence in the biofluid, and whether the protein-level evidence supports the proposed tissue of origin. It integrates protein and tissue resources including GTEx (Jiang et al., 2020), PaxDB, HPA, and UniProt to provide an evidence-based assessment of protein provenance, specificity, localisation, and biological plausibility.

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
       sample_type = "csf")


# 2. Narrow to a specific brain region
triage(c("NEFL", "GFAP", "SPANXA1"), sample_type = "csf", context_tissue = "Cerebellum")


# 3. Visualise a gene's tissue profile
plot_gene("TREM2")
```

## Scope

concordR audits inferential claims made from biofluid proteomics: plasma, CSF, urine. It catches the gap between where you measured and what you're claiming. The four gates in `triage()` are:

1. **Compartment plausibility** - is the protein's subcellular localisation consistent with the biofluid? 
2. **Tissue-of-origin** - is the protein detected at protein level in the tissue implied by research?
3. **Concordance** - does RNA reliably predict protein for this gene?
4. **Tissue specificity** - is the protein broadly expressed, undermining tissue-of-origin claims?

## Companion to

> Finney & Shvetcov (2026). An integrated atlas of RNA to protein concordance across human tissues and cell types. https://www.biorxiv.org/content/10.64898/2026.07.13.738336v1

The atlas paper provides the foundation. This package provides the assessment tool.


