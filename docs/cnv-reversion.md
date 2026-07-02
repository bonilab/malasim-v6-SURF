# CNV Reversion

## Summary

Malasim supports a daily copy-number reversion step for CNV genes, run
independently of the drug-driven forward mutation path. The first shipped
configuration targets `Pfmdr1` and `Pfplasmepsin`, but the code is generic for
any gene with `max_copies > 1`.

In this version reversion decrements the copy number by one per successful
event (for example `2 -> 1`, `3 -> 2`, ...). The daily probability is
`mutation_probability_per_locus * effective_cnv_reversion_multiplier`, and the
step is only attempted when active drugs do not currently select for a higher
copy number on the affected gene.

## Configuration

CNV reversion is configured under `genotype_parameters`.

```yaml
genotype_parameters:
  mutation_probability_per_locus: 0.001
  default_cnv_reversion_multiplier: 0.5

  pf_genotype_info:
    - chromosome: 5
      genes:
        - name: "Pfmdr1"
          max_copies: 2
          cnv_reversion_multiplier: 0.5

    - chromosome: 14
      genes:
        - name: "Pfplasmepsin"
          max_copies: 2
          cnv_reversion_multiplier: 0.5
```

### Multiplier Precedence

The effective reversion multiplier is resolved in this order:

1. `gene.cnv_reversion_multiplier` if present and `>= 0`
2. `genotype_parameters.default_cnv_reversion_multiplier` if present and `>= 0`
3. disabled for that gene if neither is set

### Validation

- Both multipliers must lie in `[0, 10]`. Values outside this range cause
  `Config::validate_all_cross_field_validations()` to throw.
- The upper bound is intentionally above `1.0` because the reverse mutation
  rate can be substantially higher than the forward
  `mutation_probability_per_locus`. A multiplier of `1.0` keeps the reversion
  rate equal to the forward mutation rate; values above `1.0` make reversion
  faster than the corresponding forward step. The chosen ceiling of `10` is
  a guardrail to catch typos and absurd values; raise it in code if a real
  configuration needs to exceed it.
- `cnv_reversion_multiplier` may only be set on genes with `max_copies > 1`.
  Setting it on a non-CNV gene causes a validation error.
- The configured `multiplier * mutation_probability_per_locus` is not
  re-clamped to `[0, 1]` at runtime; values above `1.0` will cause the
  per-gene roll to always succeed on the days where the reversion gate
  passes. This is intentional.

## Behavior

- Reversion is implemented as `n -> n - 1` (clamped at `1`). It is not
  restricted to `2 -> 1` even though all shipped CNV genes use `max_copies = 2`.
- Daily reversion probability per CNV gene is
  `mutation_probability_per_locus * effective_cnv_reversion_multiplier`.
- A gene is considered protected from reversion on a given day if any active
  drug with concentration `> 0` has a CNV EC50 factor pair for that gene in
  which `factors[1] > factors[0]` and `factors[1] > 1.0`. The check is
  gene-local and only inspects the gene's own `cnv_multiplicative_effect_on_EC50`
  entries.
- Only the precomputed list of CNV-capable genes (those with
  `max_copies > 1`) is iterated. The list is built once when
  `PfGenotypeInfo` is set via `set_pf_genotype_info()`.
- Successful reversions are recorded through the existing mutation tracking
  path (`Model::get_mdc()->record_1_mutation(...)`) and contribute to mutation
  event counts and genotype-transition logs.

## Call Sites

- `Person::update()` calls
  `all_clonal_parasite_populations_->apply_cnv_reversion(drugs_in_blood_.get())`
  immediately after `update_by_drugs()`.
- `SingleHostClonalParasitePopulations::apply_cnv_reversion()` iterates the
  current blood-stage parasite populations and delegates the actual mutation
  decision to `Genotype::perform_cnv_reversion(...)`.
- `Genotype::perform_cnv_reversion()` returns the same genotype instance when
  no gene reverts, so the caller only records a mutation event when the
  genotype actually changes.
- Reverted sequences are re-resolved through `Model::get_genotype_db()` so the
  canonical, interned genotype instance is used.

## Notes

- The forward drug-driven mutation path in `update_by_drugs()` now uses a
  uniform `proposed = old +/- 1` rule that is clamped to `[1, max_copies]`,
  replacing the previous boundary-only branching. The clamping semantics are
  preserved.
- Reverse fitness pressure continues to be applied through the existing
  `cnv_daily_crs` values; no new fitness mechanism is added.
- The current implementation is intentionally gene-local and lightweight, so
  it scales cleanly as more CNV genes are added later.
