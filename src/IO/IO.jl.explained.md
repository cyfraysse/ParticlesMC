# IO.jl — Plain-language explanation

## What this file does

This module handles reading and writing particle configurations to disk. It
defines two main actions:

- **`load_configuration` / `load_chains`** — read a file (XYZ, EXYZ, or LAMMPS
  format) and reconstruct the particle system from it.
- **`store_trajectory` / `store_lastframe` / `store_backup`** — write the
  current system state to a file in one of those formats.

---

## The three write functions and when each is used

| Function | When called | Includes bonds? |
|---|---|---|
| `store_trajectory` | Every scheduled step (trajectory output) | No |
| `store_lastframe` | At the end of a run | Yes (for Molecules) |
| `store_backup` | At each checkpoint (restart files) | Yes (for Molecules) |

**Why do some include bonds and others don't?**

Bonds never change during a simulation — they are a fixed property of the
molecule topology. Trajectory files are written thousands of times and can
be huge, so skipping bonds keeps them compact. But any file that is meant
to be *loaded back* to resume a simulation must include bonds, because
`load_configuration` needs them to reconstruct the molecular system.

`store_lastframe` was always correct (it included bonds). `store_backup` was
added to match: it writes positions **and** bonds so checkpoint files are
fully self-contained restart files.

---

## What `load_configuration` does

1. Reads the header line to get `N` (number of particles) and box dimensions.
2. Reads the metadata line to find which columns contain species, position,
   and (for molecules) molecule index.
3. Reads `N` particle lines.
4. **For molecular systems**: detects the `molecule` column in the header and
   calls `read_bonds` to read the bond list that follows the particle data.
5. Returns a dictionary with all of this.

**Important**: if the file has a `molecule` column in the header but no bonds
section after the particle data, `read_bonds` raises `"No bonds found in the
file"`. This is what happened before `store_backup` was fixed.

---

## What `load_chains` does

Takes a file or directory, calls `load_configuration` on each file found,
then calls `System(...)` to build a fully initialised simulation system
(including energy computation and neighbour list). You never need to build
the energy or neighbour list yourself — `System` does it from the positions.

---

## The `store_bonds` helper

Counts the total number of bonds (each bond `i–j` is stored once, as `i j`),
writes that count, then writes each pair on its own line. This is the format
`read_bonds` expects when loading.
