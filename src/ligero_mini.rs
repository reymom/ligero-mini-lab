use crate::field::F;
use crate::rs_code::{Domain, eval_poly, is_low_degree};

/// Tiny witness for statement: z = x * y + x in F.
///
/// In a real Ligero instance, x, y, z would be *vectors* (length m), but we start
/// with one coordinate and later extend to multiple.
#[derive(Clone, Debug)]
pub struct Witness {
    pub x: F,
    pub y: F,
    pub z: F,
}

impl Witness {
    pub fn new(x: u32, y: u32) -> Self {
        let x_f = F::from(x);
        let y_f = F::from(y);
        let z_f = x_f * y_f + x_f;
        Self {
            x: x_f,
            y: y_f,
            z: z_f,
        }
    }
}

/// Vector witness for m gates: for gate i,
///   z_i is meant to satisfy  z_i = x_i * y_i + x_i.
#[derive(Clone, Debug)]
pub struct WitnessVec {
    pub xs: Vec<F>,
    pub ys: Vec<F>,
    pub zs: Vec<F>,
}

impl WitnessVec {
    /// Build an honest witness from x,y coordinates; z_i := x_i*y_i + x_i.
    pub fn from_xy(xs: Vec<F>, ys: Vec<F>) -> Self {
        assert_eq!(xs.len(), ys.len());
        let zs = xs.iter().zip(&ys).map(|(&x, &y)| x * y + x).collect();
        Self { xs, ys, zs }
    }

    pub fn len(&self) -> usize {
        self.xs.len()
    }
}

/// Parameters for the baby Ligero instance.
#[derive(Clone, Debug)]
pub struct LigeroParams {
    /// Number of multiplication gates (m in the whiteboard).
    pub m: usize,
    /// RS codeword length (number of positions / "virtual parties").
    pub n: usize,
    /// Maximum degree allowed for each row polynomial.
    pub deg_bound: usize,
}

/// A single row of the Ligero tableau U: evaluations of a low-degree polynomial at domain points.
#[derive(Clone, Debug)]
pub struct Row {
    pub values: Vec<F>,
}

/// A tiny Ligero tableau: 3m rows, each a low-degree codeword in F^n.
///
/// Rows are grouped as (u_1, u_2, u_3) for gate 1, ..., (u_{3m-2}, u_{3m-1}, u_{3m}) for gate m.
#[derive(Clone, Debug)]
pub struct Tableau {
    pub params: LigeroParams,
    pub domain: Domain,
    pub rows: Vec<Row>, // length = 3m (for now: just one multiplication gate -> 3 rows)
}

impl Tableau {
    /// Build an honest tableau for a single gate z = x * y + x.
    /// We introduce an intermediate product wire t = x * y so that:
    ///  - mult. test enforces t = x * y,
    ///  - linear test enforces z = t + x.
    pub fn from_witness_single_gate(w: &Witness, params: LigeroParams) -> Self {
        assert_eq!(params.m, 1, "from_witness_single_gate assumes m = 1");
        let domain = Domain::new_n(params.n);

        let x = w.x;
        let y = w.y;
        let t = x * y;
        let z = w.z; // should be t + x for an honest witness

        // Constant polynomials for each wire.
        let fx = vec![x];
        let fy = vec![y];
        let ft = vec![t];
        let fz = vec![z];

        let row_x = Row {
            values: eval_poly(&fx, &domain),
        };
        let row_y = Row {
            values: eval_poly(&fy, &domain),
        };
        let row_t = Row {
            values: eval_poly(&ft, &domain),
        };
        let row_z = Row {
            values: eval_poly(&fz, &domain),
        };

        Tableau {
            params,
            domain,
            rows: vec![row_x, row_y, row_t, row_z],
        }
    }

    /// Build an honest tableau for m gates using a vector witness.
    /// For each gate i:
    ///   t_i = x_i * y_i
    ///   z_i is taken from w.zs[i] (expected to be x_i*y_i + x_i).
    pub fn from_witness_vec(w: &WitnessVec, params: LigeroParams) -> Self {
        assert_eq!(params.m, w.len(), "params.m must match witness length");
        let domain = Domain::new_n(params.n);
        let mut rows = Vec::with_capacity(4 * params.m);

        for i in 0..params.m {
            let x = w.xs[i];
            let y = w.ys[i];
            let t = x * y;
            let z = w.zs[i];

            let fx = [x];
            let fy = [y];
            let ft = [t];
            let fz = [z];

            let row_x = Row {
                values: eval_poly(&fx, &domain),
            }; // constant
            let row_y = Row {
                values: eval_poly(&fy, &domain),
            };
            let row_t = Row {
                values: eval_poly(&ft, &domain),
            };
            let row_z = Row {
                values: eval_poly(&fz, &domain),
            };

            rows.push(row_x);
            rows.push(row_y);
            rows.push(row_t);
            rows.push(row_z);
        }

        Tableau {
            params,
            domain,
            rows,
        }
    }

    /// V_p = (1, r, r², ..., r^{3m-1}) · U
    /// Check: is V_p a low-degree (< deg_bound) RS codeword?
    pub fn proximity_test(&self, r: F) -> bool {
        let mut v_p = vec![F::from(0); self.domain.len()];

        // V_p[j] = ∑_{i=0}^{3m-1} r^i * u_i[j]  for each domain point j
        for j in 0..self.domain.len() {
            let mut sum = F::from(0);
            let mut power = F::from(1);
            for row in &self.rows {
                sum = sum + row.values[j] * power;
                power = power * r;
            }
            v_p[j] = sum;
        }

        is_low_degree(&v_p, &self.domain, self.params.deg_bound)
    }

    /// Multiplication test: checks all gates satisfy t_i = x_i * y_i.
    ///
    /// For gate i, we use rows (4i, 4i+1, 4i+2) = (x_i, y_i, t_i).
    /// We form:
    ///   V_m[j] = Σ_{i=0}^{m-1} r^i * (x_i[j] * y_i[j] - t_i[j])
    /// and require V_m to be identically zero.
    pub fn multiplication_test(&self, r: F) -> bool {
        let m = self.params.m;
        let n = self.domain.len();
        let mut v_m = vec![F::from(0); n];

        for j in 0..n {
            let mut sum = F::from(0);
            let mut power = F::from(1); // r^0, r^1, ..., r^{m-1}

            for gate in 0..m {
                let idx_x = 4 * gate;
                let idx_y = 4 * gate + 1;
                let idx_t = 4 * gate + 2;

                let xj = self.rows[idx_x].values[j];
                let yj = self.rows[idx_y].values[j];
                let tj = self.rows[idx_t].values[j];

                let diff = xj * yj - tj;
                sum = sum + power * diff;
                power = power * r;
            }

            v_m[j] = sum;
        }

        // In Ligero, V_m should decode to the zero polynomial; here we enforce the
        // strongest possible condition: every coordinate is zero.
        v_m.iter().all(|&v| v == F::from(0))
    }

    /// Linear test: checks all gates satisfy the linear relation z_i = t_i + x_i.
    ///
    /// For gate i, we use rows (4i, 4i+2, 4i+3) = (x_i, t_i, z_i) and form:
    ///   V_l[j] = Σ_{i=0}^{m-1} r^i * (z_i[j] - t_i[j] - x_i[j])
    /// We require V_l to be identically zero.
    pub fn linear_test(&self, r: F) -> bool {
        let m = self.params.m;
        let n = self.domain.len();
        let mut v_l = vec![F::from(0); n];

        for j in 0..n {
            let mut sum = F::from(0);
            let mut power = F::from(1); // r^0, r^1, ..., r^{m-1}

            for gate in 0..m {
                let idx_x = 4 * gate;
                let idx_t = 4 * gate + 2;
                let idx_z = 4 * gate + 3;

                let xj = self.rows[idx_x].values[j];
                let tj = self.rows[idx_t].values[j];
                let zj = self.rows[idx_z].values[j];

                // Linear constraint: z - t - x = 0
                let diff = zj - tj - xj;
                sum = sum + power * diff;
                power = power * r;
            }

            v_l[j] = sum;
        }

        v_l.iter().all(|&v| v == F::from(0))
    }
}

#[test]
fn cheating_prover() {
    let honest_w = Witness::new(2, 3); // x=2, y=3, z=8 correct
    let cheat_w = Witness {
        z: F::from(9),
        ..honest_w
    }; // same x,y, wrong z

    let params = LigeroParams {
        m: 1,
        n: 7,
        deg_bound: 1,
    };
    let honest_tab = Tableau::from_witness_single_gate(&honest_w, params.clone());
    let cheat_tab = Tableau::from_witness_single_gate(&cheat_w, params);

    let r = F::from(42);

    // Honest prover passes everything.
    assert!(honest_tab.multiplication_test(r));
    assert!(honest_tab.linear_test(r));

    // Cheater still passes multiplication (t = x*y), but fails the linear relation z = t + x.
    assert!(cheat_tab.multiplication_test(r));
    assert!(!cheat_tab.linear_test(r));
}
