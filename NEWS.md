# concordR 0.1.0

* Initial release.
* Gene-level atlas: 11,154 genes with concordance class, protein confidence,
  detection rate, mechanistic tier, tau specificity index, and subcellular
  localisation flags.
* Tissue-resolved atlas: per-(gene, tissue) concordance with IHC, GTEx MS,
  and PaxDb protein evidence.
* `triage()`: four-gate defensibility audit for drug-target claims from
  biofluid proteomics (compartment, tissue-of-origin, concordance,
  specificity).
* `flag_genes()`: pre-analysis sanity check with traffic-light flags.
* `plot_gene()`: tissue-resolved concordance profile visualisation.
