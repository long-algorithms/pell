# Pell — fast C++/GMP solver

A standalone solver for the **fundamental positive solution** of Pell's equation

$$
X^2-DY^2=1,\qquad D>1\text{ a nonsquare integer}.
$$

The algorithm is a Pell-specialized **reduction of indefinite binary quadratic forms**, implemented through a closely related continued-fraction recurrence. Its mathematical viewpoint comes from Conway's *topograph* and the classical reduction theory of Lagrange and Gauss; Will Jagy's exposition of reduced-form cycles is especially useful for connecting that viewpoint to an explicit algorithm [1–3]. The code is **not** a literal transcription of Conway's drawings or a general binary-quadratic-form solver.

It derives from the author's earlier [`conway_full.cpp`](https://github.com/octonion/puzzles/blob/master/continued/conway_full.cpp), with specialized arithmetic, earlier stopping, and an optional square-factor/lifting strategy.

## Algorithm

### Reduced forms and continued fractions

For a form $(A,B,C)$ of discriminant $\Delta=B^2-4AC>0$, the right-neighbor change of variables described by Jagy [2] has the form

$$
(A,B,C)\longmapsto
(C,\,-B+2C\delta,\,A-B\delta+C\delta^2).
$$

The associated change-of-basis matrix is sparse, so updating the accumulated transformation does not require general matrix multiplication. For Pell's equation the principal form is $X^2-DY^2$, of discriminant **$4D$** (not necessarily the fundamental discriminant of the field).

The solver uses a normalized, integer-only state for the equivalent continued-fraction walk. Writing $r=\lfloor\sqrt D\rfloor$, its successive steps can be expressed as

$$
q=\left\lfloor\frac{h+r}{k}\right\rfloor,\qquad
h'=qk-h,\qquad k'=\frac{D-(h')^2}{k},
$$

with the signed residual tracked separately. Alongside this small state, two consecutive convergent coefficients satisfy

$$
(w,v)\longmapsto(v,w+qv).
$$

The corresponding numerator is recovered from the current reduced state instead of carrying a second row of arbitrary-precision matrix entries. The solver stops at the small residuals $\pm1,\pm2,\pm4$ and applies exact norm-one finishing identities. Its fast machine-word path also uses **continued-fraction midpoint symmetry** to stop earlier when possible; see [4] for midpoint criteria and performance comparisons. These stopping rules are part of the Pell specialization, not a general algorithm for representing arbitrary integers by binary quadratic forms.

### Arithmetic optimizations

The reduced-form parameters use checked 64-bit arithmetic for $D<10^{15}$, with a GMP fallback for larger inputs. The growing convergents use GMP (`mpz_t`/`mpz_class`). Several sparse transformations are accumulated into a checked machine-word **block matrix** before touching the growing integers; fused GMP multiply-adds and swaps minimize allocations and copies. Integer-root and divisibility tests are exact. This is a sequential block-accumulation scheme, not generic Strassen matrix multiplication.

### Optional square-factor reduction and modular lifting

For $D=d s^2$, the solver can first find the fundamental unit $\varepsilon=u+v\sqrt d$, then determine the **least** $n\ge1$ such that the coefficient of $\sqrt d$ in $\varepsilon^n$ is divisible by $s$. It computes this index from prime-power congruences and their least common multiple. An exact compact solution is then $(d,s,u,v,n)$; huge coordinates need not be constructed merely to report that solution.

A small, dependency-free prime table and early exact perfect-square tests detect square factors automatically for 64-bit inputs. The decomposition is complete for $D\le10^{12}$; at larger sizes detection is deliberately bounded and may return a nonmaximal $s$, but always checks $s^2\mid D$. **Full factorization is not required**: the direct solver works for nonsquarefree inputs too. PARI's `core2`/`core2partial` and FLINT's factorization routines are possible comparison implementations, not dependencies.

## Build and usage

Requires a C++17 compiler and GMP (`libgmp-dev` on Debian/Ubuntu); Python 3 is used by the command-line tests. The fast path assumes a 64-bit `unsigned long` (LP64) and supports GCC/Clang `unsigned __int128`.

```sh
make
make test
./conway_fast 61
./conway_fast 410286423278424 --no-print
./conway_fast 410286423278424 --no-print --expand
```

By default, the solver prints solutions of at most **1,000 decimal digits per coordinate** and summarizes larger ones. When square-factor lifting yields a provably oversized result, auto mode and `--no-print` avoid expanding it at all. For direct solves without a compact representation, the current code still constructs the unit.

| Option | Effect |
| --- | --- |
| `--print` / `--no-print` | Force or suppress coordinate display |
| `--print-limit N` | Change the automatic digit threshold |
| `--expand` | Construct full coordinates without necessarily printing them |
| `--verify` | Expand and verify $X^2-DY^2=1$ |
| `--direct` | Skip square-factor detection and lifting |
| `--square-factor S` | Supply a checked square factor, $S^2\mid D$ |
| `--no-discover-square-factor` | Alias for `--direct` |
| `--discover-square-factor` | Explicitly enable bounded detection (64-bit $D$ only) |
| `--repeat N` | Time $N$ computations and report average microseconds |

Reported computation time **excludes** output formatting and decimal printing; `--verify` performs its norm check after the timed interval. There is no cattle-specific flag: Archimedes' cattle problem is a benchmark for the general Pell solver, not a separate operating mode.

## Tests and status

`make test` checks square-factor extraction against an independent reference for all integers from 2 through 100,000, exercises structured and larger inputs, compares direct and lifted solutions on small parameters, and tests printing/compact modes. GitHub Actions builds and tests with GCC and Clang.

**Status:** an experimental, performance-oriented implementation. Exact norm verification establishes that an expanded pair solves the equation, but does **not by itself** prove that it is fundamental. Minimality relies on the reduced-form stopping criteria and the least-exponent modular-lifting argument; a full proof-level audit remains desirable. Benchmark results are environment-dependent and should not be interpreted as a demonstrated speed record against PARI or other libraries.

## References

1. **J. H. Conway**, assisted by **F. Y. C. Fung**, *The Sensual (Quadratic) Form*, Mathematical Association of America, 1997. The topograph and its river provide a geometric interpretation of indefinite-form reduction. [Publisher](https://doi.org/10.5948/UPO9781614440253).
2. **Will Jagy**, “[Binary quadratic forms over the rational integers and class numbers of quadratic fields](https://math.blogoverflow.com/2014/08/23/binary-quadratic-forms-over-the-rational-integers-and-class-numbers-of-quadratic-%EF%AC%81elds/),” *Stack Exchange Mathematics Blog*, August 23, 2014. Explicit reduced-form neighbor steps, cycles, and their connection with Pell equations. See also [Jagy's answer on general binary-form representation](https://math.stackexchange.com/questions/739752/how-to-solve-binary-form-ax2bxycy2-m-for-integer-and-rational-x-y), a broader problem than the one implemented here.
3. **D. A. Buell**, *Binary Quadratic Forms: Classical Theory and Modern Computations*, Springer, 1989. Indefinite-form reduction, cycles, and automorphisms. [DOI](https://doi.org/10.1007/978-1-4612-4542-1).
4. **K. Matthews, J. Robertson, and J. White**, “Midpoint criteria for solving Pell's equation using the nearest square continued fraction,” *Mathematics of Computation* **79** (2010), 485–499. Midpoint techniques and comparisons with regular and nearest-integer continued fractions. [DOI](https://doi.org/10.1090/S0025-5718-09-02286-8).

This software is released under the MIT License.