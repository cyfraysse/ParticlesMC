# ParticlesMC.jl — Plain-language explanation

## What this file does

This is the entry point: it reads a TOML configuration file, builds the
simulation object, and runs it. The `@main` macro (from Comonicon) turns
`particlesmc(params)` into a command-line program you call as:

```
julia --project -e 'using ParticlesMC; particlesmc("params.toml")'
```

---

## How a simulation is built

1. **Read the TOML** — extract system properties (temperature, density, model)
   and simulation settings (steps, moves, outputs).
2. **Load chains** — each configuration file becomes one independent
   replica (called a "chain") of the system.
3. **Build the algorithm list** — each entry in the list is a named tuple that
   describes one algorithm (Metropolis, StoreCallbacks, …). The `Simulation`
   constructor reads these and instantiates the algorithm objects.
4. **Create `Simulation`** and call `run!`.

---

## The restart/checkpoint feature

### How it is activated

Add `trestart = N` to the `[simulation]` section of your TOML. That is all.
The code then:

- At every `trestart` MC steps, saves a checkpoint.
- At startup, looks for existing checkpoints and restarts automatically if
  any are found.
- Optionally calls a shell command (`submit_command = "sbatch job.sh"`) at
  each checkpoint to queue the next job on the cluster.

### Auto-detection logic (`find_latest_backup`)

At startup (when `trestart` is present), the code scans
`output_path/chains/1/` for files named `restart_t{number}.exyz` or
`restart_t{number}.xyz`. It picks the one with the largest step number.
If none exists → fresh start. If found → `t_start` is set to that number.

### What happens on a restart

1. Each chain is loaded from its own checkpoint file
   (`output_path/chains/c/restart_t{t_start}.exyz`).
2. If `lastphiframe.dat` exists (molecular system with rotation tracking),
   it is loaded with `load_phi_frame` and passed to `ComputeRotation` so
   the rotation accumulators resume from exactly where they stopped.
3. All output algorithms open their files in **append** mode (`"a"`) instead of
   overwrite mode (`"w"`), so previously accumulated data is preserved.
4. The `Simulation` is created with `t_start` so the main loop starts at
   `t_start + 1` instead of `1`.
5. The scheduler counters for every algorithm are fast-forwarded to the first
   scheduled event after `t_start`, so no measurement is missed or doubled.

### The `multi_origins` schedule is preserved automatically

The `multi_origins` schedule is fully determined by `(tw, N, steps)` — these
never change between runs. When restarting from `t_start`, the same schedule
is rebuilt and the counters skip past `t_start`. Origins that started before
`t_start` still have their future observation times in the schedule, and they
will be hit as the loop advances.

### Job resubmission (`maybe_resubmit`)

If `submit_command` is set in the TOML and `t < steps`, Julia calls:
```julia
run(Cmd(split(submit_command)))
```
This submits the next job to the cluster. The submitted job runs the same
TOML — it will find the latest checkpoint and resume from there.

**Each job runs exactly `trestart` steps, then exits cleanly.** This is
enforced by passing `job_steps = min(steps, t_start + trestart)` to the
`Simulation` constructor instead of the full `steps`. The simulation reaches
`job_steps`, writes one checkpoint, submits the next job, and terminates
normally — no wall-time kill needed, no overlap between jobs.

The chain looks like:
- Job 001: t = 0 → trestart, checkpoint, submit Job 002, exit
- Job 002: t = trestart → 2×trestart, checkpoint, submit Job 003, exit
- …
- Last job: t = (n-1)×trestart → steps, checkpoint (t = steps → no submit), exit

### What is NOT saved (and why it is OK)

The RNG (random number generator) state is not saved. This means the exact
sequence of random numbers after a restart will differ from a hypothetical
single long run. The physics is still correct: the particles start from the
right positions, and the Markov chain property of Metropolis is preserved.
The only thing lost is *bit-for-bit reproducibility from a fixed seed*,
which is rarely required in statistical physics.

---

## Helper functions added for restart

| Function | What it does |
|---|---|
| `find_latest_backup(output_path, n_chains)` | Scans disk for the most recent checkpoint, returns `(t, extension)` |
| `maybe_resubmit(t, steps, submit_command)` | Calls `sbatch` (or any command) if `t < steps` and a command is configured |

---

## parse_schedule

Translates the `scheduler_params` block from the TOML into a concrete sorted
list of timesteps. Three modes:

| TOML key | Schedule type |
|---|---|
| `linear_interval = N` | Every N steps |
| `log_base = b` | Log-spaced with base b |
| `multi_origins = {tw=…, N=…}` | Multi-origin scheme (see below) |

### Multi-origins schedule

Used for computing time-correlation functions efficiently. Instead of measuring
at fixed intervals, you place "origins" every `tw` steps and, from each origin,
measure at a log-spaced set of delays. This gives good statistics at both short
and long times with far fewer total measurements than a purely linear schedule.

The schedule is fully deterministic from `(tw, N, steps)` and is rebuilt
identically on restart.
