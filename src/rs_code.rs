/// Minimal Reed–Solomon-style utilities over F.
/// Domain points are fixed once for the entire lab.
use crate::field::F;

/// RS domain used for all codewords: α_0, ..., α_{n-1}.
#[derive(Clone, Debug)]
pub struct Domain {
    pub points: Vec<F>,
}

impl Domain {
    /// Create a domain with points {0,1,...,n-1} in F.
    pub fn new_n(n: usize) -> Self {
        let points = (0..n).map(|i| F::from(i as u32)).collect();
        Self { points }
    }

    pub fn len(&self) -> usize {
        self.points.len()
    }
}

/// Evaluate polynomial with coefficients coeffs[0..=deg] on all domain points.
/// coeffs[j] is the coefficient of X^j.
pub fn eval_poly(coeffs: &[F], domain: &Domain) -> Vec<F> {
    domain
        .points
        .iter()
        .map(|&x| {
            let mut acc = F::from(0);
            let mut power = F::from(1);
            for &c in coeffs {
                acc = acc + c * power;
                power = power * x;
            }
            acc
        })
        .collect()
}

/// Interpolate polynomial from evaluations on given domain, assuming degree < deg_bound.
/// This uses a straightforward O(n^3) Gaussian elimination — fine for tiny n.
pub fn interpolate_poly(values: &[F], domain: &Domain, deg_bound: usize) -> Option<Vec<F>> {
    let n = domain.len();
    if values.len() != n {
        return None;
    }
    if deg_bound + 1 > n {
        // Under-determined, multiple polynomials; in the lab we avoid this case.
        return None;
    }

    // Build Vandermonde matrix V[i][j] = α_i^j for i in [0..n), j in [0..deg_bound]
    let mut mat = vec![vec![F::from(0); deg_bound + 1]; n];
    for i in 0..n {
        let x = domain.points[i];
        let mut power = F::from(1);
        for j in 0..=deg_bound {
            mat[i][j] = power;
            power = power * x;
        }
    }

    // Right-hand side is the 'values' vector.
    let mut rhs = values.to_vec();

    // Solve V * coeffs = values via Gaussian elimination over F.
    let m = deg_bound + 1;
    let mut row = 0;
    for col in 0..m {
        // Find pivot.
        let mut pivot = None;
        for r in row..n {
            if mat[r][col] != F::from(0) {
                pivot = Some(r);
                break;
            }
        }
        if pivot.is_none() {
            // Column is all zeros; skip. In practice, tiny domains should avoid this.
            continue;
        }
        let pivot = pivot.unwrap();
        mat.swap(row, pivot);
        rhs.swap(row, pivot);

        // Normalize pivot row.
        let inv = mat[row][col].inv();
        for c in col..m {
            mat[row][c] = mat[row][c] * inv;
        }
        rhs[row] = rhs[row] * inv;

        // Eliminate below.
        for r in (row + 1)..n {
            if mat[r][col] != F::from(0) {
                let factor = mat[r][col];
                for c in col..m {
                    mat[r][c] = mat[r][c] - factor * mat[row][c];
                }
                rhs[r] = rhs[r] - factor * rhs[row];
            }
        }
        row += 1;
        if row == m {
            break;
        }
    }

    // Back substitution to solve for coeffs[0..m).
    let mut coeffs = vec![F::from(0); m];
    for i in (0..m).rev() {
        // Find leading 1 in row i, if any.
        let mut lead_col = None;
        for c in 0..m {
            if mat[i][c] != F::from(0) {
                lead_col = Some(c);
                break;
            }
        }
        if let Some(col) = lead_col {
            let mut acc = rhs[i];
            for c in (col + 1)..m {
                acc = acc - mat[i][c] * coeffs[c];
            }
            coeffs[col] = acc;
        }
    }

    // Quick consistency check: re-evaluate and compare.
    let check_vals = eval_poly(&coeffs, domain);
    if check_vals == values {
        Some(coeffs)
    } else {
        None
    }
}

/// Check if a vector is a valid RS codeword for some polynomial of degree < deg_bound.
pub fn is_low_degree(values: &[F], domain: &Domain, deg_bound: usize) -> bool {
    interpolate_poly(values, domain, deg_bound).is_some()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::field::F;

    #[test]
    fn rs_roundtrip() {
        let domain = Domain::new_n(5);
        // f(X) = 3 + 2X + X^2 over F
        let coeffs = vec![F::from(3), F::from(2), F::from(1)];
        let evals = eval_poly(&coeffs, &domain);
        let recovered = interpolate_poly(&evals, &domain, 2).unwrap();
        assert_eq!(coeffs, recovered);
    }

    #[test]
    fn low_degree_check() {
        let domain = Domain::new_n(5);
        let coeffs = vec![F::from(1), F::from(0)]; // constant 1
        let evals = eval_poly(&coeffs, &domain);
        assert!(is_low_degree(&evals, &domain, 1));
    }
}
