# ViennaRNA patches

This repository contains patches for specific ViennaRNA releases.
Patches are not officially supported and are maintained by our group.

Each patch is designed to be applied to the corresponding upstream release and
contains all modifications introduced up to that patch version.

## Patch

**Base version:** ViennaRNA 2.7.2

**Patch series:**

1. `patch-viennarna-2.7.2-1`
2. `patch-viennarna-2.7.2-2`

Apply them consecutively. Advanced users who only need the continuous
soft-constraint changes can stop after patch 1.

## Summary

Patch 1 removes the discretization of soft constraints in partition function
calculations while preserving the original integer arithmetic for MFE
calculations.

Patch 2 fixes suboptimal structure enumeration with unpaired soft constraints
added through `vrna_sc_add_up()`.

In the original ViennaRNA implementation, soft constraints added through

- `vrna_sc_add_up()`
- `vrna_sc_add_bp()`

are immediately rounded to 0.01 kcal/mol and stored as integers. This affects
partition function calculations and base-pair probabilities, which therefore
become discontinuous functions of the applied soft constraints.

Patch 1 changes the internal storage of soft constraints to floating-point
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

- **Suboptimal structure enumeration (`subopt`):**
  - completed structures are re-evaluated when soft constraints are present;
  - final output is filtered by the evaluated energy, so structures outside the
    requested energy band are not returned;
  - duplicate structures caused by this soft-constraint backtracking issue are
    avoided.

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

## Testing

After installation, the tests in this repository can be run with

```bash
python -m unittest discover -s tests -v
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

Suboptimal structure enumeration (`subopt`) still generates candidates from the
rounded MFE dynamic-programming tables. Patch 2 re-evaluates and filters final
structures when soft constraints are active; it does not make `subopt` itself a
continuous-energy enumerator.
