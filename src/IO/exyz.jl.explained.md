# exyz.jl — Plain-language explanation

## What is EXYZ?

EXYZ (Extended XYZ) is a text format for particle configurations. Each frame
has two header lines followed by one line per particle:

```
1290
Lattice="32.9 0.0 0.0 0.0 32.9 0.0 0.0 0.0 0.0" Properties=molecule:I:1:species:S:1:pos:R:3 Time=400
1 A  0.123  4.567  8.901
1 A  1.234  5.678  9.012
...
```

Line 1: number of particles.
Line 2: key=value pairs describing the box, column layout, and metadata.
Lines 3+: one particle per line, columns defined by `Properties=`.

---

## How columns are described

The `Properties=` field uses a compact notation:

```
name:type:dimension
```

For example `pos:R:3` means "a column group called `pos`, type Real, 3 numbers
wide". Multiple groups are separated by `:` with no extra separator, so the
full string looks like `molecule:I:1:species:S:1:pos:R:3`.

`parse_column_string` splits this on `:` and groups it into triples
`(name, type, dimension)`, building an ordered dictionary that maps each
name to its column index and width.

---

## The regex bug that was fixed

The metadata line can have extra key=value pairs after `Properties=`, like
`Time=400`. The original regex was:

```julia
match(r"Properties=(.*)", metadata_line)
```

`(.*)` is greedy and matches everything to the end of the line, so it would
capture `molecule:I:1:species:S:1:pos:R:3 Time=400`. Then `parse_column_string`
would receive `"3 Time=400"` as the last token and fail trying to
`parse(Int, "3 Time=400")`.

The fix:

```julia
match(r"Properties=(\S+)", metadata_line)
```

`\S+` matches only non-whitespace characters, stopping at the space before
`Time=400`. Now only the column description is captured.

**When does this matter?** Only when *loading* an EXYZ file that has extra
metadata after `Properties=`. Checkpoint files written by `StoreBackups`
always include `Time=t` at the end, which triggered this bug on restart.

---

## 2D vs 3D systems

The `Lattice=` field always has 9 numbers (a 3×3 matrix), even for 2D systems
where the third diagonal entry is 0. The reader extracts the diagonal
`[Lx, Ly, Lz]` and then truncates to match the position dimensionality:
if `pos:R:2`, box becomes `[Lx, Ly]`.

---

## read_t_from_lastframe — the restart hook

```julia
function Arianna.read_t_from_lastframe(path::String, ::EXYZ)
    open(path) do f
        readline(f)      # skip N
        comment = readline(f)
        m = match(r"Time=(\d+)", comment)
        return parse(Int, m.captures[1])
    end
end
```

This overrides the generic Arianna function for EXYZ files. Arianna's `detect_restart` calls
it when looking for `lastframe.exyz` at the start of a new job.

`Time=t` in the comment line is the MC step at which the job ended (written by
`StoreLastFrames.finalise`). Reading it back gives `t_start` for the next job.

The default Arianna version expects `t` as the first comma-separated field (DAT format).
EXYZ stores it differently — embedded in a metadata string — so it needs its own parser.
