# Pell equations — experimental C++/GMP solver

This standalone implementation solves

\[
X^2-DY^2=1,\qquad D>1\text{ nonsquare},\quad X,Y>0,
\]

using a Pell-specialized reduced-form/continued-fraction computation, optional
square-factor reduction, and modular lifting. It was developed from the
`conway_full.cpp` experiments in [`octonion/puzzles`](https://github.com/octonion/puzzles/tree/master/continued).

**Dependencies:** a C++17 compiler, GMP (`gmp` and `gmpxx`), and Python 3 for
CLI tests. Square-factor detection uses only standard C++ arithmetic and an
included table of primes up to 10,000; **PARI and FLINT are not dependencies**.
The fast arithmetic path uses `unsigned __int128` and the build assumes a
64-bit `unsigned long` (LP64). This is an experimental solver, not yet a
published, proof-audited library.

## Build and use

```sh
make
make test
./conway_fast 61
./conway_fast 410286423278424 --no-print
./conway_fast 410286423278424 --no-print --expand
```

By default, solutions with at most 1,000 decimal digits in each coordinate
are printed. Larger solutions are summarized without decimal conversion; if
a square-factor reduction is available, they can be represented exactly by
a smaller unit and a lifting exponent, without constructing the large integers.

| Option | Meaning |
| --- | --- |
| `--print` | Force expansion and decimal display |
| `--no-print` | Suppress display (and avoid expansion when a compact representation exists) |
| `--print-limit N` | Set the automatic printing cutoff in decimal digits |
| `--expand` | Construct full coordinates, but do not force decimal display |
| `--verify` | Expand and verify the norm equation |
| `--direct` or `--no-discover-square-factor` | Bypass automatic square-factor discovery |
| `--square-factor S` | Use a supplied square factor; require `S*S` to divide `D` |
| `--discover-square-factor` | Explicitly request automatic detection (64-bit inputs only) |
| `--repeat N` | Repeat the computation and report average microseconds per run |

Timing ends **before** summary formatting or decimal printing. The `--verify`
norm check occurs after the timed interval. The direct path currently expands
the fundamental unit even under `--no-print`; the compact representation is
available when a useful square factor is extracted.

## Square-factor detection

For positive input `D`, seek

\[
D=d s^2,\quad s\ge1.
\]

The **maximal** choice of `s` makes `d` squarefree. We use a bounded,
dependency-free detector rather than fully factoring every input:

1. Visit primes from the included, compile-time table (at most 10,000),
   extracting their entire powers. For each exponent `e`, multiply `s` by
   `p^(e/2)` (integer quotient) and keep the odd part in the remainder.
2. Immediately after removing a prime power, test whether the remaining
   cofactor is a perfect square. If it is `r*r`, multiply `s` by `r` and stop.
3. Stop at the configured prime limit, or earlier when `p*p > remainder`.
   Finally, perform one more exact square test on the remaining cofactor.

Both square and cube root routines use floating point **only for an initial
guess**; 128-bit integer products correct the guess and make all acceptance
decisions exact. The detector always returns a mathematically valid identity
`D == d*s*s`, even when `d` has not been proved squarefree.

### When is the decomposition complete?

Let the trial-prime ceiling be `B` (10,000 by default). For
`D <= B^3`, the returned factor `s` is maximal: every prime up to
`cuberoot(D)` is tested, unless the remaining cofactor is already a square.
After those primes have been removed, a remaining nonsquare number has at most
two prime factors, counted with multiplicity. A repeated prime would make
that remainder a square, which the final test detects. Therefore no
undetected squared prime factor remains.

For the default bound, **every `D <= 10^12` is fully decomposed**. For larger
64-bit `D`, square factors with large prime bases might be missed; then `d`
can still contain a square factor. This is deliberate: expensive general
factorization would make some easy Pell instances much slower. The direct
Pell solver works for nonsquarefree `d` as well. If `D` does not fit in an
unsigned 64-bit integer, the CLI automatically uses the direct GMP path
instead of running this detector. A user-supplied square factor is checked
for exact divisibility, but it need not be maximal.

### Why the early square test matters for Archimedes

For the cattle-problem Pell parameter,

\[
D=410286423278424=4729494\cdot9314^2.
\]

Removing small prime factors through 353 leaves the square `4657^2`.
The detector can stop at that point, rather than trial-dividing through
4,657. The compact Pell result then records the reduced parameter and the
least required unit exponent, `2329`.

### Relation to existing libraries

- **PARI/GP** provides `core2` (complete squarefree decomposition) and
  `core2partial` (bounded extraction). Useful for validation or an optional
  backend, but not needed for this implementation.
- **FLINT** provides fast 64-bit integer factorization via `n_factor`, among
  other factorization routines. Full factorization is not necessarily faster
  than our bounded detector on typical Pell parameters.

The small-prime and early-square approach is chosen for low overhead, not
because it is proven fastest for every integer. Large inputs that require
complete squarefree decomposition may benefit from a different algorithm.

## Performance and testing

On the **previous experiment's local container CPU**, the isolated Archimedes
square-factor detection took approximately 0.18 µs with the early-square
shortcut, versus 3.08 µs for the older trial-division detector. A separate
5000-run compact-solver benchmark measured approximately 2.21 µs per call
versus 5.08 µs for the previous version. These measurements were **not made
on a Ryzen 9 5950X**. They are not a controlled comparison with PARI or FLINT;
benchmark on your target machine to establish practical performance.

`make test` runs the detector against a simple independent reference for all
integers from 2 through 100,000, structured inputs through `10^12`, and larger
bounded cases checking the exact identity. It also compares direct and
square-reduced Pell solutions on small inputs, and exercises the print/compact
CLI behavior. A correct Pell norm alone does not establish *minimality*;
minimality also depends on the reduced-form termination and modular-order
arguments, which should be documented and audited separately before a
production release.

This standalone project contains **no cattle-specific CLI flag**. Archimedes'
problem is a useful benchmark; its herd-count computation belongs in a
separate example rather than the general Pell solver.
