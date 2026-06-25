# rotation.jl — Plain-language explanation

## What this file does

This file tracks how much each molecule has rotated during the simulation.
It is only relevant for molecular systems (ortho-terphenyl, trimers, …) — for
atomic systems (`Atoms`) it is never used.

---

## Key concept: the rotation accumulator

Imagine you freeze a molecule at time 0 and take a photo of it.
Then, every time the molecule rotates past a threshold angle θ_T, you update
the "reference photo" and accumulate the total rotation so far.
At any later time, you can report the total angle the molecule has rotated
since the beginning of the simulation.

This is done with two objects per molecule:

| Variable | Meaning |
|---|---|
| `R_ref[k][m]` | Reference orientation of molecule `m` under threshold `k` (a 3×3 rotation matrix) |
| `Φ_acc[k][m]` | Accumulated rotation vector so far (a 3D vector whose length = total angle rotated, in radians) |
| `system.Φ[k][m]` | Current total rotation = `Φ_acc + rotation since last R_ref update` |

The threshold `θ_T` controls how often you update the reference frame. A small
threshold gives very precise tracking but updates more often.

---

## What a rotation matrix is

A 3×3 matrix `R` represents a rotation in 3D space. If you multiply `R` by a
vector, you get the rotated vector. The matrix `R_ref' * R_current` gives you
the rotation that happened *between* the reference and now (the apostrophe means
"transpose", which for rotation matrices equals "inverse").

`rotation_vector(R)` converts that relative rotation matrix into a single 3D
vector: the direction is the rotation axis, and the length is the angle in
radians (this is called the "axis-angle" or "Rodrigues" representation).

---

## ComputeRotation

```
ComputeRotation(chains; theta_T=[π/4])
```

`theta_T` is a list of thresholds (one per "measurement channel"). Using
multiple thresholds lets you study rotation at different angular resolutions
at once.

**`initialise`** (called once at the start of each job):
- Computes the current body frame for every molecule from its atomic positions
  (loaded either from the original config or from `lastframe.exyz`).
- If `simulation.t_start > 0` (restart): looks for `lastphiframe.dat` in
  `simulation.path/chains/c/` and loads `R_ref`, `Φ_acc`, `Φ` from it.
- Otherwise (fresh start): sets `R_ref = current frame`, `Φ_acc = 0`, `Φ = 0`.

The restart detection is done inside `initialise` by reading `simulation.t_start`
directly — no `initial_state` field is needed on the struct.

**`make_step!`** (called at every scheduled MC step):
- Recomputes every molecule's current body frame.
- For each (molecule, threshold): computes the rotation since `R_ref`, checks
  if it exceeds `θ_T`. If yes, adds to `Φ_acc` and updates `R_ref`.
- Writes the current total rotation `Φ = Φ_acc + Φ_current` into `system.Φ`.

---

## StorePhiTrajectories

Writes `phitrajectories_k.dat` (one file per threshold `k`, one per chain).
Each line: `N_mol` header then one line per molecule with `m Φx Φy Φz`.

**On restart:** `initialise` checks `simulation.t_start > 0` and opens files
in append mode. The `store_first` write is also skipped (the value at `t_start`
is already in the file from the previous job).

---

## StoreLastPhiFrame and load_phi_frame

`StoreLastPhiFrame` saves the *complete internal state* of `ComputeRotation`
to `lastphiframe.dat` in `finalise` — at the end of each job. This is the
rotation checkpoint: everything needed to resume tracking from exactly where
it left off.

`load_phi_frame(path)` is the inverse: it reads `lastphiframe.dat` and
reconstructs `R_ref`, `Φ_acc`, and `Φ` as Julia arrays. This is called
inside `ComputeRotation.initialise` when a restart is detected.

**Why only `finalise` and not `make_step!`?**
Because each job runs exactly `t_restart` steps and always ends cleanly
(the loop stops at `t_start + t_restart`, then `finalise` fires). There is
no need to write mid-job since the job never gets killed partway through.

**Why do we need to save `R_ref` and not just `Φ_acc`?**
Because `R_ref` is the reference orientation used to detect the *next* crossing
of the threshold. Without it, the first step after a restart would compute the
rotation relative to the molecule's orientation at t=0, not at the last
crossing — giving a wrong accumulated angle.

---

## Body frame computation

`body_frame(r1, r2, r3, L)` builds a right-handed orthonormal basis from the
positions of three atoms within a molecule, applying the minimum-image
convention (periodic boundary conditions). This gives a coordinate frame that
is fixed to the molecule regardless of its absolute position in the box.

---

## Numerical note

`rotation_vector(R)` has three special cases to avoid division by zero:
- When `sin(θ) ≈ 0` and `cos(θ) > 0`: the rotation is nearly identity → return zero vector.
- When `sin(θ) ≈ 0` and `cos(θ) < 0`: the rotation is nearly 180° → extract axis from diagonal of `R`.
- Otherwise: standard Rodrigues formula.
