# CNV Reversion For Pfmdr1 And Pfplasmepsin

## Summary

Add a separate daily CNV reversion path for copy-number genes so copy number
can decrement outside the current forward drug-mutation flow. The first
version supports `Pfmdr1` and `Pfplasmepsin` through per-gene reversion
multipliers, with optional fallback to `default_cnv_reversion_multiplier`,
uses `mutation_probability_per_locus * multiplier` as the reversion rate, and
only attempts reversion when there are no active drugs or when current drugs
do not select for the higher copy.

## Implementation Changes

- Add gene-level optional `cnv_reversion_multiplier` to genotype configuration.
- Add optional `default_cnv_reversion_multiplier` under `genotype_parameters`
  as a fallback for CNV genes that do not specify their own multiplier.
- Precompute the set of CNV-capable genes (`max_copies > 1`) inside
  `PfGenotypeInfo` via `refresh_cnv_gene_indices()` so the daily reversion
  step does not need to rescan every gene.
- Add a genotype/helper path for CNV reversion that:
  - only inspects genes in the precomputed CNV gene index list
  - decrements copy number by one (for example `2 -> 1`, `3 -> 2`, ...) with
    a floor of `1`
  - resolves the effective multiplier from the gene-level value first, then
    `default_cnv_reversion_multiplier`
  - computes the daily reversion probability as
    `mutation_probability_per_locus * effective_multiplier`
- Add a separate within-host daily function, distinct from `update_by_drugs`,
  for example `apply_cnv_reversion()`. Call it from `Person::update()` as its
  own step.
- Reversion eligibility rule for v1:
  - allow reversion if the host has no active drugs in blood, or
  - allow reversion if active drugs do not select for a higher copy number for
    that gene
- Use a simple gene-local selection heuristic for "drug selects for a higher
  copy":
  - if any active drug with concentration `> 0` has a CNV EC50 factor pair
    where `factors[1] > factors[0]` and `factors[1] > 1.0`, treat copy 2 (or
    higher) as selected
  - otherwise treat the gene as not selected and allow reversion attempts
- Refactor the forward mutation in `update_by_drugs()` to use a uniform
  `proposed = old +/- 1` rule clamped to `[1, max_copies]`. The clamping
  semantics from the previous boundary-only branching are preserved.
- Add validation in `Config::validate_all_cross_field_validations()` so both
  reversion multipliers must lie in `[0, 10]`. The upper bound is well above
  `1.0` because the reverse mutation rate is expected to sometimes exceed
  the forward `mutation_probability_per_locus`; the ceiling is a guardrail
  against typos and absurd values, not a physical limit. The runtime
  probability is not re-clamped, so `multiplier > 1.0` is honored (and will
  make the per-gene roll always succeed on days the reversion gate passes).
- Record successful reversions through the existing mutation tracking path so
  they contribute to mutation event counts and genotype-transition logs.
- Represent the higher reverse fitness penalty through existing
  `cnv_daily_crs` values rather than adding a new fitness mechanism in code.

## Key Code Touchpoints

- `GenotypeParameters` / config decoding: add `cnv_reversion_multiplier`,
  `default_cnv_reversion_multiplier`, and the `PfGenotypeInfo::CnvGeneIndex`
  precomputed index.
- `Config::validate_all_cross_field_validations()`: validate the new
  multiplier fields against `[0, 10]` and the `max_copies > 1` constraint.
- `Genotype` (and a small free helper in `Genotype.cpp`): implement the
  copy-number reversion candidate generation and drug-selection check.
- `SingleHostClonalParasitePopulations` plus `Person::update()`: add and call
  the separate daily CNV reversion step, and record the resulting mutation
  events.

## Test Plan

- Config parsing test for both `cnv_reversion_multiplier` and
  `default_cnv_reversion_multiplier`.
- Reversion helper test: copy 2 with multiplier set can become copy 1.
- Higher-copy test: copy 3 with `max_copies = 3` can become copy 2.
- Fallback test: CNV gene without its own multiplier can still revert via
  `default_cnv_reversion_multiplier`.
- Negative test: non-CNV genes are ignored.
- Daily host update behavior:
  - reversion can occur with no drugs present
  - reversion can occur when active drugs do not select for copy 2 for that
    gene
  - reversion does not occur when an active drug selects for copy 2 for that
    gene
- Forward mutation behavior preserved:
  - copy `2` with `max_copies = 3` can step up to `3` or down to `1` under
    `perform_mutation_by_drug`
  - copy `3` with `max_copies = 3` steps down to `2`
- Coverage for both current genes: `Pfmdr1` and `Pfplasmepsin`.
- Mutation logging test: successful reversion is recorded as a mutation
  event.
- Validation tests: multipliers above `10` are rejected; multipliers above
  `1` but at or below `10` are accepted (smoke-test at, e.g., `2.0` and
  `5.0`).

## Assumptions And Defaults

- Only `Pfmdr1` and `Pfplasmepsin` receive explicit `cnv_reversion_multiplier`
  in current example config, but the code supports any future CNV gene the
  same way.
- v1 supports copy-number reversion as `n -> n - 1` with a floor of `1`,
  generalized from the original `2 -> 1` plan once the implementation
  converged.
- Reversion rate uses `mutation_probability_per_locus * multiplier`, not a
  separate hard-coded baseline such as `0.001983`, unless config is later
  changed to make that necessary.
- The reverse rate is allowed to exceed the forward rate, so the multiplier
  range is `[0, 10]` rather than `[0, 1]`. The exact upper bound is a
  guardrail; raise it in `Config.cpp` if a real configuration needs to go
  higher.
- "Drug does not select for copy 2" is intentionally implemented as a simple
  per-gene CNV EC50-factor check, not a full phenotype-validation model.
- The precomputed CNV gene index keeps the daily cost bounded by the number
  of CNV genes, which is two in the current configuration.
