# ViennaRNA patches

This repository contains patches for specific ViennaRNA releases.
Patches are not officially supported and are maintained by our group.

Each patch is designed to be applied to the corresponding upstream release and
contains all modifications introduced up to that patch version.

## Patch

**Base version:** ViennaRNA 2.7.2

**Patch:** `patch-viennarna-2.7.2-1`

## Summary

This patch removes the discretization of soft constraints in partition function
calculations while preserving the original integer arithmetic for MFE
calculations.

In the original ViennaRNA implementation, soft constraints added through

- `vrna_sc_add_up()`
- `vrna_sc_add_bp()`

are immediately rounded to 0.01 kcal/mol and stored as integers. This affects
partition function calculations and base-pair probabilities, which therefore
become discontinuous functions of the applied soft constraints.

This patch changes the internal storage of soft constraints to floating-point
values while keeping the same internal energy units (0.01 kcal/mol).

The consequences are:

- **Partition function (PF):**
  - soft constraints are represented with full floating-point precision;
  - Boltzmann factors are computed from the exact values;
  - partition functions and base-pair probabilities become continuous functions
    of the supplied soft constraints.

- **Minimum Free Energy (MFE):**
  - unchanged with respect to upstream ViennaRNA;
  - energies are rounded only when preparing the MFE dynamic-programming tables.

Therefore, MFE predictions remain fully compatible with the original ViennaRNA
implementation, while PF calculations gain substantially improved numerical
precision.

## Files modified

- `src/ViennaRNA/constraints/soft.c`
- `src/ViennaRNA/constraints/soft.h`

## Applying the patch

```bash
git clone https://github.com/ViennaRNA/ViennaRNA.git
cd ViennaRNA
git checkout v2.7.2
git apply patch-viennarna-2.7.2-1.patch
```

## Limitations

This patch affects only soft constraints handled through

- `vrna_sc_add_up()`
- `vrna_sc_set_up()`
- `vrna_sc_add_bp()`
- `vrna_sc_set_bp()`

(and their comparative counterparts).

Stack soft constraints (`vrna_sc_add_stack()` / `vrna_sc_set_stack()`) still use
the original integer discretization.

Suboptimal structure enumeration (`subopt`) still uses the rounded MFE energies,
as in the original ViennaRNA implementation.
