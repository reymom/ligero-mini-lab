/// Small prime field F_p.
/// We keep p tiny so we can debug by hand while still exercising modular arithmetic.
#[derive(Copy, Clone, Debug, PartialEq, Eq, Hash)]
pub struct F {
    /// Always kept reduced modulo P.
    value: u32,
}

impl F {
    /// Small prime modulus; something convenient for hand debugging.
    pub const P: u32 = 97;

    /// Create a new field element reduced modulo P.
    pub fn new(x: u32) -> Self {
        Self { value: x % Self::P }
    }

    /// Return the canonical representative in [0, P).
    pub fn as_u32(self) -> u32 {
        self.value
    }

    /// Add two field elements.
    pub fn add(self, other: F) -> F {
        let mut s = self.value as u64 + other.value as u64;
        s %= Self::P as u64;
        F::new(s as u32)
    }

    /// Subtract two field elements.
    pub fn sub(self, other: F) -> F {
        let p = Self::P as i64;
        let mut s = self.value as i64 - other.value as i64;
        s %= p;
        if s < 0 {
            s += p;
        }
        F::new(s as u32)
    }

    /// Multiply two field elements.
    pub fn mul(self, other: F) -> F {
        let mut p = self.value as u64 * other.value as u64;
        p %= Self::P as u64;
        F::new(p as u32)
    }

    /// Compute multiplicative inverse using extended Euclidean algorithm.
    /// Panics if called on zero.
    pub fn inv(self) -> F {
        assert!(self.value != 0, "cannot invert zero in field");
        // Extended Euclidean algorithm for a * x + p * y = gcd(a, p) = 1
        let mut t: i64 = 0;
        let mut new_t: i64 = 1;
        let mut r: i64 = Self::P as i64;
        let mut new_r: i64 = self.value as i64;

        while new_r != 0 {
            let q = r / new_r;
            let tmp_t = t - q * new_t;
            t = new_t;
            new_t = tmp_t;

            let tmp_r = r - q * new_r;
            r = new_r;
            new_r = tmp_r;
        }

        if r != 1 {
            panic!("element not invertible, gcd != 1");
        }
        if t < 0 {
            t += Self::P as i64;
        }
        F::new(t as u32)
    }

    /// Exponentiation by square-and-multiply: self^exp.
    pub fn pow(self, mut exp: u64) -> F {
        let mut base = self;
        let mut acc = F::new(1);
        while exp > 0 {
            if exp & 1 == 1 {
                acc = acc.mul(base);
            }
            base = base.mul(base);
            exp >>= 1;
        }
        acc
    }
}

use std::ops::{Add, Div, Mul, Neg, Sub};

impl Add for F {
    type Output = F;
    fn add(self, rhs: F) -> F {
        self.add(rhs)
    }
}

impl Sub for F {
    type Output = F;
    fn sub(self, rhs: F) -> F {
        self.sub(rhs)
    }
}

impl Mul for F {
    type Output = F;
    fn mul(self, rhs: F) -> F {
        self.mul(rhs)
    }
}

impl Div for F {
    type Output = F;
    fn div(self, rhs: F) -> F {
        self.mul(rhs.inv())
    }
}

impl Neg for F {
    type Output = F;
    fn neg(self) -> F {
        F::new(0) - self
    }
}

impl From<u32> for F {
    fn from(x: u32) -> F {
        F::new(x)
    }
}

#[cfg(test)]
mod tests {
    use super::F;

    #[test]
    fn basic_field_arithmetic() {
        let a = F::new(5);
        let b = F::new(96); // -1 mod 97
        assert_eq!((a + b).as_u32(), 4); // 5 + (-1) = 4
        assert_eq!((a * b).as_u32(), (5 * 96) % F::P);
        let inv = F::new(5).inv();
        assert_eq!((F::new(5) * inv).as_u32(), 1);
    }
}
