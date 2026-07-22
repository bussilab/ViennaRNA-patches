# ViennaRNA patches

This repository contains patches for specific ViennaRNA releases.
Patches are not officially supported and are maintained by our group.

Each patch is independent and applies directly to the corresponding upstream
release. By default, patches are applied in numerical order, but either patch
may also be applied alone.

## Patches

**Base version:** ViennaRNA 2.7.2

- `patch-viennarna-2.7.2-1.patch`: continuous soft-constraint values
- `patch-viennarna-2.7.2-2.patch`: correct `subopt()` backtracking with
  unpaired soft constraints

## Summary

### Patch 1

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

### Patch 2

This patch fixes suboptimal-structure backtracking when unpaired soft
constraints contribute to multiloop decompositions. Without the correction,
`subopt()` can omit valid structures, return structures outside the requested
energy band, or produce duplicate structures.

## Files modified

Patch 1:

- `src/ViennaRNA/constraints/soft.c`
- `src/ViennaRNA/constraints/soft.h`

Patch 2:

- `src/ViennaRNA/subopt/subopt.c`

## Applying the patch

Download the corresponding ViennaRNA release tarball from

https://www.tbi.univie.ac.at/RNA/

For example:

```bash
wget https://www.tbi.univie.ac.at/RNA/download/sourcecode/2_7_x/ViennaRNA-2.7.2.tar.gz
tar xzf ViennaRNA-2.7.2.tar.gz
cd ViennaRNA-2.7.2

patch -p1 < ../patch-viennarna-2.7.2-1.patch
patch -p1 < ../patch-viennarna-2.7.2-2.patch

./configure \
    --prefix=$HOME/viennarna \
    --without-doc \
    --without-doc-html \
    --without-doc-pdf

make -j$(nproc)
make install
```

To use only one fix, apply only its corresponding patch command. Neither patch
depends on the other.

## Limitations

Patch 1 affects only soft constraints handled through

- `vrna_sc_add_up()`
- `vrna_sc_set_up()`
- `vrna_sc_add_bp()`
- `vrna_sc_set_bp()`

(and their comparative counterparts).

Stack soft constraints (`vrna_sc_add_stack()` / `vrna_sc_set_stack()`) still use
the original integer discretization.

Suboptimal structure enumeration (`subopt`) still reports rounded MFE energies,
as in the original ViennaRNA implementation. Patch 2 corrects which structures
are reached during backtracking; it does not change that reporting convention.
