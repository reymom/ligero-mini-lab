# ligero-mini-lab 🧪

A mini Rust lab that implements a **baby version of Ligero's three tests** — proximity, multiplication, and linear — for a tiny arithmetic circuit, and uses them to see **soundness amplification** in action.

> This repo is meant as a companion to the ZK Hack S3M5 whiteboard session on **Ligero** (Muthu Venkitasubramaniam).  
> The goal is not to build a production SNARK, but to make the _matrix U + three tests_ picture concrete in a few hundred lines of Rust.

---

## 1. What this lab does

We work over a small prime field $F_p$ (currently `p = 97`) and a tiny Reed–Solomon (RS) code, then:

- Encode the circuit
  $$
  z = x · y + x
  $$
  as `m` independent **multiplication gates** over $F_p$.
- Arrange the wires into a **tableau** `U` whose rows are RS codewords:
  - for each gate $i$, we have rows $(x_i, y_i, t_i, z_i)$ with
    - $t_i = x_i · y_i$,
    - $z*i$ _meant_ to satisfy $(z_i = t_i + x_i)$.
- Implement three Ligero-style tests on the rows of `U`:
  1. **Proximity test** $V_p$: checks “rows look like RS codewords”.
  2. **Multiplication test** $V_m$: checks all $t_i = x_i · y_i$.
  3. **Linear test** $V\_\\ell$: checks all $z_i = t_i + x_i$.
- Craft an explicit **cheating strategy** that:
  - keeps all rows in the RS code (passes proximity),
  - preserves linear constraints (passes linear),
  - but breaks multiplication on two gates in a correlated way.
- Run a Monte-Carlo experiment to see that:
  - a single random challenge `r` catches the cheat with probability ≈ `1/p`,
  - repeating with fresh `r` drives cheating probability down like `≈ (1/p)^k`.

This reproduces, in toy form, the **soundness amplification** story that appears in Ligero and more broadly in IOP-style protocols.

---

## 2. High-level picture

### 2.1 Field and code

- Field: a tiny prime field $F_p$ with `p = 97` for easy debugging.
- Code: a Reed–Solomon code defined by a domain
  ${\alpha*0, …, \alpha*{n-1}} \subset F_p$, here `α_j = j`.
- Codewords: vectors of the form
  $$
  (f(\alpha_0), …, f(\alpha_{n-1})) \in F_p^n,
  $$
  for polynomials `f(X)` with `deg(f) < d`.

The RS machinery in `rs_code.rs` provides:

- `eval_poly(coeffs, domain) -> Vec<F>` – evaluate a polynomial on the domain.
- `interpolate_poly(values, domain, deg_bound) -> Option<Vec<F>>` – recover the unique polynomial of degree `< deg_bound` if it exists.
- `is_low_degree(codeword, domain, deg_bound) -> bool` – convenience wrapper:
  - interpolate,
  - re-evaluate,
  - check equality.

Ligero uses this structure to ask: _“Is this vector close to the RS code?”_ — in this lab we go with the strongest version and ask: _“Is it exactly an RS codeword of degree `< d`?”_.

### 2.2 Tableau U and rows

We use a **tableau** `U` with `4m` rows, grouped by gate:

- For gate `i` (0-based):
  - row `4i + 0`: `x_i` (constant RS codeword),
  - row `4i + 1`: `y_i` (constant RS codeword),
  - row `4i + 2`: `t_i = x_i*y_i`,
  - row `4i + 3`: `z_i` (intended: `z_i = t_i + x_i`).

Each row is the evaluation of a **low-degree polynomial** (in the lab, just degree 0) on the RS domain.

---

## 3. Code structure

```text
src/
  field.rs      # Tiny prime field F_p (p = 97), + basic arithmetic, pow, inv via EEA
  rs_code.rs    # Mini Reed–Solomon utilities: domain, eval, interpolate, low-degree check
  ligero_mini.rs
                # Core lab:
                #   - Witness / WitnessVec
                #   - LigeroParams
                #   - Row / Tableau
                #   - Proximity, multiplication, linear tests
  lib.rs        # Re-exports modules
  main.rs       # Experiments: single-gate smoke test + multi-gate soundness experiment
```

### 3.1 `field.rs`

- Implements `F` with:
  - `add`, `sub`, `mul`, `neg`,
  - `pow(self, exp: u64)` via square-and-multiply,
  - `inv(self)` using the **Extended Euclidean Algorithm**.
- This is intentionally tiny and hand-checkable; it can later be swapped for a serious field (e.g. Goldilocks, BabyBear via arkworks) to compare performance / ergonomics.

### 3.2 `rs_code.rs`

- `Domain::new_n(n)` – domain `[0, 1, 2, ..., n-1]` in `F_p`.
- `eval_poly` – evaluate `f(X)` at all points.
- `interpolate_poly` – small Vandermonde-style interpolation with degree bound check.
- `is_low_degree` – the core primitive for the **proximity test**.

### 3.3 `ligero_mini.rs`

Core types:

- `Witness` – scalar witness for a single gate (x, y, z) with `z = x*y + x`.
- `WitnessVec` – vector witness for `m` gates:
  - `xs: Vec<F>`, `ys: Vec<F>`, `zs: Vec<F>` with `zs[i] = xs[i]*ys[i] + xs[i]` for honest witnesses.
- `LigeroParams` – tiny parameter bundle (m, n, deg_bound).
- `Row` – one RS codeword (a vector in `F^n`).
- `Tableau` – the matrix `U`, containing `4m` rows, each `Row`.

Builders:

- `Tableau::from_witness_single_gate` – build a 1-gate tableau `(x, y, t, z)` with constant polynomials.
- `Tableau::from_witness_vec` – multi-gate version for `m > 1`.

Tests on rows:

- `proximity_test(r: F) -> bool`
  - Computes `V_p = Σ r^i u_i` over rows `u_i` and returns `is_low_degree(V_p)`.
  - Only checks **code structure**; honest and cheating provers both pass if they keep rows as valid RS codewords.
- `multiplication_test(r: F) -> bool`
  - For each gate `i`:
    - takes rows `x_i, y_i, t_i`,
    - forms `diff_i[j] = x_i[j] * y_i[j] - t_i[j]`,
    - sets `V_m[j] = Σ r^i * diff_i[j]`.
  - Requires `V_m` to be identically zero.
- `linear_test(r: F) -> bool`
  - For each gate `i`:
    - takes rows `x_i, t_i, z_i`,
    - forms `diff_i[j] = z_i[j] - t_i[j] - x_i[j]`,
    - sets `V_ℓ[j] = Σ r^i * diff_i[j]`.
  - Requires `V_ℓ` to be identically zero.

These three tests capture the same separation as in the Ligero whiteboard:

- proximity: “Are you in the code?”
- multiplication: “Are the products correct?”
- linear: “Do the linear constraints on wires hold?”

---

## 4. Building and running

### 4.1 Prerequisites

Clone and test

```bash
git clone https://github.com/reymom/ligero-mini-lab.git
cd ligero-mini-lab

# Run unit tests (field, RS, basic Ligero tests)
cargo test
```

You should see something like:

```text
running 4 tests
test field::tests::basic_field_arithmetic ... ok
test ligero_mini::cheating_prover ... ok
test rs_code::tests::low_degree_check ... ok
test rs_code::tests::rs_roundtrip ... ok

test result: ok. 4 passed; 0 failed; ...
```

### 4.2 Run the experiments

`main.rs` currently runs the single-gate smoke test and the multi-gate soundness experiment. To run it:

```bash
cargo run
```

Expected shape of the multi-gate output (values may vary):

```text
Honest (single r): prox=true, mult=true, lin=true
Cheat  (single r): prox=true, mult=false, lin=true

Estimating acceptance probabilities over 5000 trials:
  rounds = 1: honest_accept ≈ 1.0000, cheat_accept ≈ 0.008600
  rounds = 2: honest_accept ≈ 1.0000, cheat_accept ≈ 0.000000
  rounds = 3: honest_accept ≈ 1.0000, cheat_accept ≈ 0.000000
```

Interpretation:

- Honest prover passes (≈1.0) regardless of the number of rounds.
- Cheating prover passes about `≈ 1/p` in one round (`p = 97` here), and dies off as you repeat with fresh randomness — matching the standard “soundness ≈ (1/p)^k” story.

**Optional: single-gate smoke test**

- Builds a 1-gate tableau for `x = 2, y = 3, z = 2*3 + 2 = 8`.
- Run all three tests once and print their booleans.

This is handy to connect the math with a concrete example you can compute by hand.

## 5. How honest / cheating behaviour is wired

### 5.1 Honest prover

- For each gate `i`:

  - sample `x_i, y_i ∈ F_p`,
  - set `t_i = x_i * y_i`,
  - set `z_i = t_i + x_i`.

- Build `U` with rows `x_i, y_i, t_i, z_i` as **constant RS codewords**:

  ```rust
  let xs: Vec<F> = (0..m).map(|_| F::from(rng.gen_range(1..F::P))).collect();
  let ys: Vec<F> = (0..m).map(|_| F::from(rng.gen_range(1..F::P))).collect();

  let w_vec = WitnessVec::from_xy(xs, ys);
  let honest_tab = Tableau::from_witness_vec(&w_vec, params.clone());
  ```

- Proximity, multiplication and linear tests all pass for every `r`.

### 5.2 Correlated cheating prover

To see soundness, we construct a correlated cheat hitting two gates:

- Start from the honest tableau.
- Choose some non-zero `δ ∈ F_p`.
- Modify gate 0 and gate 1 as follows (for all columns `j`):

  ```rust
  // Gate 0
  t0'[j] = t0[j] + δ
  z0'[j] = z0[j] + δ

  // Gate 1
  t1'[j] = t1[j] - δ
  z1'[j] = z1[j] - δ
  ```

Properties:

- All rows are still constant RS codewords → proximity passes.
- For both gates, `z_i' = t_i' + x_i` still holds → linear test passes.
- Products are broken:
  - gate 0: `x0*y0 - t0' = -δ`,
  - gate 1: `x1*y1 - t1' = +δ`.

So for each column `j`:

[\\V_m[j] = (-\\delta)r^0 + (+\\delta)r^1 = \\delta (r - 1), \\] which is zero **iff** `r = 1`:

- Cheater passes one multiplication test with probability exactly `1/p` over random `r`.
- Repeating the protocol `k` times with fresh `r` multiplies probabilities: `≈ (1/p)^k`. This is the core phenomenon we observe in the Monte-Carlo experiment.

## 6. How this relates to real Ligero

This repo is deliberately tiny and does **not** implement a full Ligero SNARK.

What it does capture from the whiteboard:

- The idea of a **tableau** `U` whose rows live in a linear code (RS).
- Three tests defined as **random linear combinations** of those rows:
  - `V_p` – code structure / proximity.
  - `V_m` – multiplication constraints.
  - `V_ℓ` – linear constraints.
  - Soundness amplification via **fresh randomness** `r` across multiple rounds.

Things it does **not** include (yet): - Full **MPC-in-the-head** construction and how `U` arises from views of a secure computation. - Packed secret sharing schemes beyond the simplest `m`-gate layout. - The exact `A w = b` linear system and `\u005chat r = (1, r, r², …) A` construction. - Efficient low-degree / proximity tests with slack (we require exact RS codewords). - Zero-knowledge, commitments, or any pre-processing / public parameters.

Think of this as an **explainer lab**: enough to make the whiteboard's `U` + tests feel concrete, without committing to a full protocol.

## 7. Possible extensions

If you want to extend this lab, some natural next steps:

- Replace the toy field `F_p` with:
  - a larger prime field, - or a dedicated SNARK-friendly field from a library (e.g. arkworks).
  - Increase the RS domain size and play with **real proximity** (allow a few errors instead of requiring exact codewords).
  - Add more constraints:
    - multi-gate linear relations (`A w = b`),
    - small circuits composed of several `z = x*y + x` blocks.
  - Connect to a real MPC-in-the-head scheme:
    - generate `U` from simulated MPC views,
    - use the same three tests as part of an end-to-end argument.

## 8. References

- ZK Hack — **S3M5: The Ligero Proof System**, Muthu Venkitasubramaniam <https://zkhack.dev/whiteboard/s3m5/>
- Ames, Hazay, Ishai, Venkitasubramaniam — _Ligero: Lightweight Sublinear Arguments Without a Trusted Setup_, ePrint 2022/1608.
- Ishai, Kushilevitz, Ostrovsky, Sahai — _Zero-Knowledge from Secure Multiparty Computation_, STOC 2007.
