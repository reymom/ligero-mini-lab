use rand::Rng;

use ligero_mini_lab::field::F;
use ligero_mini_lab::ligero_mini::{LigeroParams, Tableau, Witness, WitnessVec};

/// Single-gate smoke test: z = x*y + x with x=2, y=3, z=8.
/// Uses the m=1 convenience constructor and a fixed challenge r.
fn single_gate_smoke() {
    // Tiny smoke test: build an honest tableau for x=2, y=3.
    let w = Witness::new(2, 3); // z = 2*3 + 2 = 8
    let params = LigeroParams {
        m: 1,
        n: 5,         // 5 "virtual parties" / positions
        deg_bound: 1, // constant polynomials for now
    };

    let tab = Tableau::from_witness_single_gate(&w, params);
    println!("== Single-gate smoke test ==\n");
    println!("Tableau rows: {:?}", tab.rows);

    // Fiat–Shamir challenge: r = hash(public inputs) (fixed for demo).
    let r = F::from(42);

    println!("Proximity test (V_p low-deg?): {}", tab.proximity_test(r));
    println!(
        "Multiplication test (V_m = 0?): {}",
        tab.multiplication_test(r)
    );
    println!("Linear test (constraint holds?): {}", tab.linear_test(r));
}

/// Multi-gate experiment: m>1, correlated cheating strategy that
/// passes proximity + linear tests, and passes multiplication only
/// with probability ≈ 1/p per round.
fn multigate_soundness_experiment() {
    let mut rng = rand::thread_rng();

    // ----- 1) Build an honest multi-gate tableau -----
    let m: usize = 4;
    let params = LigeroParams {
        m,
        n: 7,         // RS codeword length (tiny domain)
        deg_bound: 1, // rows are constant polynomials: degree 0 < 1
    };

    // Random x,y per gate; avoid trivial zeros just for nicer numbers.
    let xs: Vec<F> = (0..m).map(|_| F::from(rng.gen_range(1..F::P))).collect();
    let ys: Vec<F> = (0..m).map(|_| F::from(rng.gen_range(1..F::P))).collect();

    let w_vec = WitnessVec::from_xy(xs, ys);
    let honest_tab = Tableau::from_witness_vec(&w_vec, params.clone());

    // ----- 2) Build a correlated cheating tableau -----
    // Clone honest tableau, then corrupt product wires of gate 0 and gate 1
    // in opposite directions, while keeping z = t + x so linear test still holds.
    assert!(m >= 2, "need at least two gates for the correlated cheat");
    let mut cheat_tab = honest_tab.clone();
    let n = cheat_tab.domain.len();

    let delta = F::from(1); // non-zero perturbation

    let idx_t0 = 4 * 0 + 2;
    let idx_z0 = 4 * 0 + 3;
    let idx_t1 = 4 * 1 + 2;
    let idx_z1 = 4 * 1 + 3;

    for j in 0..n {
        cheat_tab.rows[idx_t0].values[j] = cheat_tab.rows[idx_t0].values[j] + delta;
        cheat_tab.rows[idx_z0].values[j] = cheat_tab.rows[idx_z0].values[j] + delta;

        cheat_tab.rows[idx_t1].values[j] = cheat_tab.rows[idx_t1].values[j] - delta;
        cheat_tab.rows[idx_z1].values[j] = cheat_tab.rows[idx_z1].values[j] - delta;
    }

    // Sanity check for a single r
    let r_test = F::from(42);
    println!("\n== Multi-gate soundness experiment ==\n");
    println!(
        "Honest (single r): prox={}, mult={}, lin={}",
        honest_tab.proximity_test(r_test),
        honest_tab.multiplication_test(r_test),
        honest_tab.linear_test(r_test),
    );
    println!(
        "Cheat  (single r): prox={}, mult={}, lin={}",
        cheat_tab.proximity_test(r_test),
        cheat_tab.multiplication_test(r_test),
        cheat_tab.linear_test(r_test),
    );

    // ----- 3) Monte-Carlo over rounds and random r -----
    let trials = 5_000;
    let rounds_list = [1usize, 2, 3];

    println!("\nEstimating acceptance probabilities over {trials} trials:");
    for rounds in rounds_list {
        let mut honest_accept = 0usize;
        let mut cheat_accept = 0usize;

        for _ in 0..trials {
            let mut honest_ok = true;
            let mut cheat_ok = true;

            for _ in 0..rounds {
                // Sample a fresh random challenge r in F each round.
                let r = F::from(rng.gen_range(0..F::P));

                let honest_round = honest_tab.proximity_test(r)
                    && honest_tab.multiplication_test(r)
                    && honest_tab.linear_test(r);

                let cheat_round = cheat_tab.proximity_test(r)
                    && cheat_tab.multiplication_test(r)
                    && cheat_tab.linear_test(r);

                honest_ok &= honest_round;
                cheat_ok &= cheat_round;
            }

            if honest_ok {
                honest_accept += 1;
            }
            if cheat_ok {
                cheat_accept += 1;
            }
        }

        let honest_rate = honest_accept as f64 / trials as f64;
        let cheat_rate = cheat_accept as f64 / trials as f64;

        println!(
            "  rounds = {rounds}: honest_accept ≈ {:.4}, cheat_accept ≈ {:.6}",
            honest_rate, cheat_rate
        );
    }
}

fn main() {
    single_gate_smoke();
    println!("\n----------------------------------------\n");
    multigate_soundness_experiment();
}
