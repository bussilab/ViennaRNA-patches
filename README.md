# ViennaRNA patches

This repository contains patches for specific ViennaRNA releases.
Patches are not officially supported and are maintained by our group.

Patches 1 and 2 apply directly to the corresponding upstream release and may
be used independently. Patch 3 applies after patches 1 and 2.

## Patches

**Base version:** ViennaRNA 2.7.2

- `patch-viennarna-2.7.2-1.patch`: continuous soft-constraint values
- `patch-viennarna-2.7.2-2.patch`: correct `subopt()` backtracking with
  unpaired soft constraints
- `patch-viennarna-2.7.2-3.patch`: fix follow-up issues in the continuous
  soft-constraint implementation in patch 1

## Summary

### Patch 1

This patch retains floating-point soft-constraint values until they are used
for partition function (PF) or minimum free energy (MFE) calculations.

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
  - the dynamic-programming tables still use integer energies in units of
    0.01 kcal/mol;
  - contributions to the same unpaired nucleotide or base pair are accumulated
    before being rounded for those tables;
  - results can differ from upstream ViennaRNA, which rounds each contribution
    when it is added. For example, two contributions of 0.003 kcal/mol to one
    unpaired nucleotide sum to 0.006 kcal/mol and give an MFE contribution of
    0.01 kcal/mol; upstream rounds each contribution to zero.

Thus, PF calculations retain the supplied precision, while MFE calculations
round the accumulated contributions to the existing integer energy scale.

**If you use patch 1, apply patch 3 as well.** Patch 1 alone leaves the
comparative base-pair add path rounded to 0.01 kcal/mol. On a fresh comparative
fold compound, that path can also crash because it checks the wrong constraint
pointer during initialization. Its MFE dimer wrapper truncates stored
fractional base-pair energies. These issues are addressed in patch 3.

### Patch 2

This patch fixes suboptimal-structure backtracking when unpaired soft
constraints contribute to multiloop decompositions. Without the correction,
`subopt()` can omit valid structures, return structures outside the requested
energy band, or produce duplicate structures.

The official ViennaRNA repository has fixed this bug on its `test_dev` branch
in [commit de9d96ba](https://github.com/ViennaRNA/ViennaRNA/commit/de9d96ba93b7d0e7df304865f2ddb8618c8e771c).
Patch 2 remains useful for the 2.7.2 release, which does not contain that fix.

### Patch 3

This source-only follow-up makes the raw storage type consistently `double`,
including builds configured with `--enable-floatpf`. It retains fractional
base-pair contributions in comparative folding, initializes comparative
constraints correctly, and preserves fractional pair contributions when the
MFE dimer wrapper copies constraints to a second fold compound. The patch also
exposes the comparative base-pair add function to the Python bindings so its
behavior can be tested there. Regression tests live in this repository's
`tests/` directory and run in its GitHub Actions workflow using the Python
bindings.

Patch 3 changes the public base-pair storage layout in `--enable-floatpf`
builds. Rebuild ViennaRNA and any compiled bindings or extensions together
when upgrading an installation that uses patch 1.

The `sc_add_bp_comparative()` Python method is a small SWIG wrapper added to
exercise the comparative base-pair path in the regression tests. It exposes
the existing C function to Python; the tests in this repository use it to
check comparative initialization and repeated fractional contributions. If
the release archive's generated Python wrapper is stale, remove its generated
files after applying the patches and configuring ViennaRNA:

```bash
rm -f interfaces/Python/RNA_wrap.cpp \
      interfaces/Python/RNA.py \
      interfaces/Python/RNA/RNA.py
```

Then run the regular `make` and `make install` commands from the ViennaRNA
source directory. The build regenerates these files and compiles the Python
extension with the wrapper.

## Files modified

Patch 1:

- `src/ViennaRNA/constraints/soft.c`
- `src/ViennaRNA/constraints/soft.h`

Patch 2:

- `src/ViennaRNA/subopt/subopt.c`

Patch 3:

- `src/ViennaRNA/constraints/soft.c`
- `src/ViennaRNA/constraints/soft.h`
- `src/ViennaRNA/mfe/mfe_wrappers.c`
- `interfaces/constraints_soft.i`

Regression tests are maintained separately in `tests/` in this repository.

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
patch -p1 < ../patch-viennarna-2.7.2-3.patch

./configure \
    --prefix=$HOME/viennarna \
    --without-doc \
    --without-doc-html \
    --without-doc-pdf

make -j$(nproc)
make install
```

Patch 2 may be used alone. Patch 3 requires patches 1 and 2 in this order.

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
