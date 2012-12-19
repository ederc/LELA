This folder contains legacy versions of the Faugère-Lachartre algorithm (directory `LEGACY_COMPLETE_VERSIONS`) and experimental data structures
for the block version (directory `BLOCK_DATA_TYPES`).

The legacy versions are completely functionnal and are mainly simpler than the block version and could be used to understand the overall parts of the algorithm
used in the block version.

* `LEGACY_COMPLETE_VERSIONS/FGL-sparse-vector-version` contains the earliest and simplest implementation of the **standard algorithm**, using only sparse vectors (no blocks).
* `LEGACY_COMPLETE_VERSIONS/FGL-sparse-vector-new-method`: the same as the above version, but implements the **new method** (new ordering of operations).
* `LEGACY_COMPLETE_VERSIONS/FG-multiline: In this version`, the data structure **multiline** is used for both the standard and the new Faugère-Lachartre.