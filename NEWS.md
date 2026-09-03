# concordR 0.2.0

The output of `triage()` changes shape and the specificity
gate changes definition.

## triage()

* Returns one column per gate and no verdict. Six columns by default:
  `gene_symbol`, `concordance`, `tissue`, `route`, `specificity`,
  `primary_tissue`. `details = TRUE` appends the evidence behind each gate.
* `verdict` and `reason` columns removed.
* Each gate is a function of the atlas data and the sample type only. No
  gate reads another gate's result.
* The tissue gate now reads `hpa_tissue`, the three-source presence vote
  (HPA IHC, GTEx MS, PaxDb) over 59 tissues. `context_tissue` is validated
  against that vocabulary. In 0.1.0 the vote was shipped but not read.
* Localisation flags are joined at run time by `ensembl_id` with symbol
  fallback. In 0.1.0 the join was missing, so the localisation gate never
  fired.
* The localisation gate is `route`: how the protein reaches the sample,
  from UniProt alone. `secreted`, `membrane`, `intracellular`, `unknown`
  for biofluids; `matched`, `unmatched`, `unknown` for subcellular
  fractions; `na` for tissue and lysate. The 0.1.0 gate asked whether a
  protein was "plausible" in a fluid, and the CSF data answered that
  everything is: a thousand intracellular proteins are in the CSF
  proteome. "Implausible" was intracellular plus not present in the
  tissue, a conjunction of two gates, which is why the gate had to read
  tissue data and needed a rescue for projection proteins. The route gate
  is independent, the rescue column is gone, and `audit_call()` handles
  the conjunction: intracellular is caution (the fluid level reports cell
  damage, which neither supports nor contradicts a target claim),
  intracellular with a tissue mismatch is a decline on the mismatch.
* The specificity gate is protein-level. It reads a tau computed over the
  highest IHC level per tissue across the HPA panel. RNA tau and protein
  breadth are shipped as evidence only. In
  0.1.0 specificity mixed a rank-based RNA tau with an any-platform breadth
  whose denominator was the atlas total; both are documented in
  `data-raw/derive_thresholds.R` as unfit to gate.
* `primary_tissue` is the set of HPA tissues tied at the gene's highest
  IHC level, comma-joined, `"multiple"` for broad genes. A single winner
  was tried twice and rejected both times: in 0.1.0 an ordinal that tied
  across many tissues, and in development a fraction-weighted score under
  which a quarter of the atlas had Adrenal gland as primary tissue because
  HPA scores it on one cell type and the tiebreak was alphabetical. Median
  six tissues share a gene's maximum level, so the set is reported and
  `audit_call()` passes a gene if any member is an accepted tissue.
* Specificity has a `"sparse"` value for genes IHC-detected in too few
  tissues for tau to describe a profile. Tau on one detection is 1.0 by
  construction, and a quarter of the atlas's specific calls rested on a
  single antibody hit. The bulk platforms decide the outcome: if they
  detect the protein in at least half of covered tissues, the antibody is
  uninformative and the gene reads `"na"`; otherwise the restriction is
  consistent across platforms but thin and the gene reads `"sparse"`.
  Both get `primary_tissue = NA`, so a single antibody hit can no longer
  produce a wrong-tissue decline. Later revised: sparse genes report their
  tied-max set in `primary_tissue` so the caution is readable (UMOD in
  urine reads `sparse | Kidney`), and `audit_call()` declines on the set
  only for specific genes.
* The word "compartment" no longer appears in the interface.

## New functions

* `audit_call()`: accept/caution/decline from the gate columns, with
  `accept_tissues` derived from the context tissue.
* `funnel_class()`: first-failure partition in a stated gate order, with an
  optional survivor table. `funnel_class(x) == "pass"` is identical to
  `audit_call(x) == "accept"`.
* `query_presence()`: the three-source presence vote for a gene list
  across the 59 HPA tissues, which is what the tissue gate reads.

## Vectorised context

* `triage()` accepts several `context_tissue` values. The output is then
  long, one row per gene per tissue, with a `context_tissue` column after
  `gene_symbol`. `audit_call()` and `funnel_class()` read each row's own
  context. Only the tissue gate and the cross-platform clause of
  concordance vary with context; the rest repeat.

## Data

* `atlas_thresholds` shipped as internal data. Every threshold `triage()`
  uses is read from it, and `data-raw/derive_thresholds.R` is the single
  place to change one. `cross_platform_agree_max` from 0.1.0 is removed:
  with a single context tissue the quantity it compared against was a
  boolean, so the clause now reads `platforms_agree == FALSE` directly.
* `concordance_atlas` gains `tau`, `protein_breadth`, `protein_tau`,
  `primary_tissue_ihc`, `n_primary_tissues`, `n_tissues_ihc`,
  `n_tissues_ihc_detected`,
  `n_tissues_bulk_detected`, `n_tissues_bulk_covered`,
  `n_tissues_profiled`. Per-gene summaries are grouped by `ensembl_id`
  rather than symbol.
* `hpa_tissue` gains `ihc_reliability_max`.

## Fixed

* `context_tissue` errored for every value in 0.1.0 because the validator
  was handed a data frame instead of a vocabulary.
* `load_atlas("both")` returns three elements; the test expecting two is
  updated.
* `stats` and `utils` declared in Imports.
* No `toupper()` over a full table inside `triage()`; row-wise loops
  removed.

# concordR 0.1.0

* Initial release. State used for the companion atlas paper, tagged
  `v0.1.0`.
