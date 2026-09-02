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
| `triage()` | Score a gene list on four independent gates. One column per gate, no verdict. |
| `audit_call()` | Accept, caution, or decline from the gate columns. |
| `funnel_class()` | First failing gate in a stated order, with an optional survivor table. |
| `query_atlas()` | Gene-level lookup: concordance class, suppression rate, mechanistic tier. |
| `query_tissue()` | Tissue-resolved lookup: per-tissue concordance and protein evidence. |
| `plot_gene()` | Tissue-resolved concordance profile for one gene. |
| `load_atlas()` | Load the atlas objects (cached, called internally). |

## Quick start

```r
library(concordR)

# 1. Gate columns for proteins nominated from CSF, no tissue implied
triage(c("NEFL", "GFAP", "SPANXA1"), sample_type = "csf")

# 2. The same list against a named brain region
res <- triage(c("NEFL", "GFAP", "SPANXA1"), sample_type = "csf",
              context_tissue = "Cerebral cortex")

# 3. Combine the gates
res$call <- audit_call(res)
funnel_class(res, table = TRUE)

# 4. Evidence behind each gate
triage("NEFL", sample_type = "csf", details = TRUE)

# 5. Visualise a gene's tissue profile
plot_gene("TREM2")
```

## Scope

concordR audits inferential claims made from biofluid proteomics: plasma, CSF, urine. It catches the gap between where you measured and what you're claiming. The four gates in `triage()` are:

1. **Concordance** - does RNA reliably predict protein for this gene?
2. **Tissue** - is the protein detected in the tissue the claim is about, by a vote across HPA IHC, GTEx MS, and PaxDb?
3. **Route** - how does the protein reach the sample? Secreted and shed proteins can carry a mechanism claim. Intracellular proteins reach a fluid by damage release, so their level reports that cells broke.
4. **Specificity** - how tissue-specific is the protein at the cell-type level, and where is it made?

Each gate is reported on its own. `audit_call()` and `funnel_class()` are the two supplied ways to combine them, and a user is free to write another.

## Companion to

> Finney & Shvetcov (2026). An integrated atlas of RNA to protein concordance across human tissues and cell types. https://www.biorxiv.org/content/10.64898/2026.07.13.738336v1

The atlas paper provides the foundation. This package provides the assessment tool.


