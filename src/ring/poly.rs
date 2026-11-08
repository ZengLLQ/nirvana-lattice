//! Polynomial representation and arithmetic over cyclotomic quotient ring R_q = Z_q[X] / (X^N + 1).

use std::fmt;
use std::ops::{Add, AddAssign, Index, IndexMut, Neg, Sub, SubAssign};
use rand::Rng;
use serde::{Deserialize, Serialize};

use crate::ring::ntt::NttContext;
use crate::ring::reduction::{
    center, freeze, mod_add, mod_mul, mod_neg, mod_sub, DEFAULT_N, DEFAULT_Q,
};

/// A polynomial in the cyclotomic quotient ring R_q = Z_q[X] / (X^N + 1).
#[derive(Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct Poly {
    pub coeffs: Vec<i32>,
    pub q: i32,
}

impl Poly {
    /// Creates a zero polynomial of degree `n` modulo `q`.
    pub fn zero(n: usize, q: i32) -> Self {
        Self {
            coeffs: vec![0; n],
            q,
        }
    }

    /// Creates a zero polynomial with default Dilithium parameters (N=256, q=8380417).
    pub fn default_zero() -> Self {
        Self::zero(DEFAULT_N, DEFAULT_Q)
    }

    /// Creates a polynomial from a slice of coefficients, canonicalizing each modulo `q`.
    pub fn from_coeffs(coeffs: &[i32], q: i32) -> Self {
        let canonical: Vec<i32> = coeffs.iter().map(|&c| freeze(c as i64, q)).collect();
        Self {
            coeffs: canonical,
            q,
        }
    }

    /// Returns the degree bound N (number of coefficients).
    pub fn len(&self) -> usize {
        self.coeffs.len()
    }

    /// Returns whether the polynomial has zero length.
    pub fn is_empty(&self) -> bool {
        self.coeffs.is_empty()
    }

    /// Checks if all coefficients are zero modulo q.
    pub fn is_zero(&self) -> bool {
        self.coeffs.iter().all(|&c| c == 0)
    }

    /// Infinity norm ||p||_infinity: maximum absolute value of centered coefficients in [-(q-1)/2, (q-1)/2].
    pub fn infinity_norm(&self) -> i32 {
        self.coeffs
            .iter()
            .map(|&c| center(c, self.q).abs())
            .max()
            .unwrap_or(0)
    }

    /// Squared Euclidean norm ||p||_2^2: sum of squared centered coefficients.
    pub fn l2_norm_sq(&self) -> i64 {
        self.coeffs
            .iter()
            .map(|&c| {
                let v = center(c, self.q) as i64;
                v * v
            })
            .sum()
    }

    /// Euclidean L2 norm ||p||_2.
    pub fn l2_norm(&self) -> f64 {
        (self.l2_norm_sq() as f64).sqrt()
    }

    /// Returns centered coefficients in [-(q-1)/2, (q-1)/2].
    pub fn centered_coeffs(&self) -> Vec<i32> {
        self.coeffs.iter().map(|&c| center(c, self.q)).collect()
    }

    /// Multiplies the polynomial by an integer scalar modulo `q`.
    pub fn scalar_mul(&self, scalar: i32) -> Self {
        let s = freeze(scalar as i64, self.q);
        let coeffs: Vec<i32> = self
            .coeffs
            .iter()
            .map(|&c| mod_mul(c, s, self.q))
            .collect();
        Self {
            coeffs,
            q: self.q,
        }
    }

    /// Negacyclic schoolbook multiplication in Z_q[X] / (X^N + 1).
    /// Used as ground-truth oracle for testing and verifying NTT multiplication.
    pub fn schoolbook_mul(&self, other: &Self) -> Self {
        assert_eq!(self.len(), other.len(), "Polynomial lengths must match");
        assert_eq!(self.q, other.q, "Moduli must match");
        let n = self.len();
        let q = self.q;
        let mut res = vec![0i64; n];

        for i in 0..n {
            for j in 0..n {
                let prod = (self.coeffs[i] as i64) * (other.coeffs[j] as i64);
                if i + j < n {
                    res[i + j] += prod;
                } else {
                    // X^N == -1 mod (X^N + 1)
                    res[i + j - n] -= prod;
                }
            }
        }

        let canonical: Vec<i32> = res.iter().map(|&v| freeze(v, q)).collect();
        Self {
            coeffs: canonical,
            q,
        }
    }

    /// Converts to NTT representation using the given NTT context.
    pub fn to_ntt(&self, ctx: &NttContext) -> NttPoly {
        assert_eq!(self.len(), ctx.n);
        assert_eq!(self.q, ctx.q);
        let mut evals = self.coeffs.clone();
        ctx.forward(&mut evals);
        NttPoly {
            evals,
            q: self.q,
        }
    }

    /// Multiplies two polynomials using Negacyclic NTT in O(N log N).
    pub fn ntt_mul(&self, other: &Self, ctx: &NttContext) -> Self {
        assert_eq!(self.len(), other.len());
        assert_eq!(self.q, other.q);
        let a_ntt = self.to_ntt(ctx);
        let b_ntt = other.to_ntt(ctx);
        let c_ntt = a_ntt.pointwise_mul(&b_ntt);
        c_ntt.to_poly(ctx)
    }

    /// Uniform random polynomial sampling with coefficients in [0, q - 1].
    pub fn random_uniform<R: Rng>(n: usize, q: i32, rng: &mut R) -> Self {
        let coeffs: Vec<i32> = (0..n).map(|_| rng.gen_range(0..q)).collect();
        Self { coeffs, q }
    }

    /// Random polynomial with small ternary coefficients in {-1, 0, 1}.
    pub fn random_ternary<R: Rng>(n: usize, q: i32, rng: &mut R) -> Self {
        let coeffs: Vec<i32> = (0..n)
            .map(|_| {
                let r: i32 = rng.gen_range(-1..=1);
                freeze(r as i64, q)
            })
            .collect();
        Self { coeffs, q }
    }

    /// Random sparse challenge polynomial with exactly `tau` non-zero coefficients in {-1, 1}.
    /// Used for Fiat-Shamir challenge generation in lattice zero-knowledge arguments.
    pub fn random_sparse_challenge<R: Rng>(n: usize, q: i32, tau: usize, rng: &mut R) -> Self {
        assert!(tau <= n, "Weight tau cannot exceed polynomial degree N");
        let mut coeffs = vec![0i32; n];
        let mut indices: Vec<usize> = (0..n).collect();

        // Partial Fisher-Yates shuffle to pick tau unique positions
        for i in 0..tau {
            let j = rng.gen_range(i..n);
            indices.swap(i, j);
            let sign = if rng.gen_bool(0.5) { 1 } else { -1 };
            coeffs[indices[i]] = freeze(sign as i64, q);
        }

        Self { coeffs, q }
    }
}

impl Index<usize> for Poly {
    type Output = i32;
    fn index(&self, index: usize) -> &Self::Output {
        &self.coeffs[index]
    }
}

impl IndexMut<usize> for Poly {
    fn index_mut(&mut self, index: usize) -> &mut Self::Output {
        &mut self.coeffs[index]
    }
}

impl Add for &Poly {
    type Output = Poly;
    fn add(self, rhs: Self) -> Self::Output {
        assert_eq!(self.len(), rhs.len(), "Polynomial dimensions must match");
        assert_eq!(self.q, rhs.q, "Moduli must match");
        let coeffs: Vec<i32> = self
            .coeffs
            .iter()
            .zip(&rhs.coeffs)
            .map(|(&a, &b)| mod_add(a, b, self.q))
            .collect();
        Poly {
            coeffs,
            q: self.q,
        }
    }
}

impl Add for Poly {
    type Output = Poly;
    fn add(self, rhs: Self) -> Self::Output {
        &self + &rhs
    }
}

impl AddAssign<&Poly> for Poly {
    fn add_assign(&mut self, rhs: &Poly) {
        assert_eq!(self.len(), rhs.len());
        assert_eq!(self.q, rhs.q);
        for i in 0..self.coeffs.len() {
            self.coeffs[i] = mod_add(self.coeffs[i], rhs.coeffs[i], self.q);
        }
    }
}

impl Sub for &Poly {
    type Output = Poly;
    fn sub(self, rhs: Self) -> Self::Output {
        assert_eq!(self.len(), rhs.len(), "Polynomial dimensions must match");
        assert_eq!(self.q, rhs.q, "Moduli must match");
        let coeffs: Vec<i32> = self
            .coeffs
            .iter()
            .zip(&rhs.coeffs)
            .map(|(&a, &b)| mod_sub(a, b, self.q))
            .collect();
        Poly {
            coeffs,
            q: self.q,
        }
    }
}

impl Sub for Poly {
    type Output = Poly;
    fn sub(self, rhs: Self) -> Self::Output {
        &self - &rhs
    }
}

impl SubAssign<&Poly> for Poly {
    fn sub_assign(&mut self, rhs: &Poly) {
        assert_eq!(self.len(), rhs.len());
        assert_eq!(self.q, rhs.q);
        for i in 0..self.coeffs.len() {
            self.coeffs[i] = mod_sub(self.coeffs[i], rhs.coeffs[i], self.q);
        }
    }
}

impl Neg for &Poly {
    type Output = Poly;
    fn neg(self) -> Self::Output {
        let coeffs: Vec<i32> = self.coeffs.iter().map(|&c| mod_neg(c, self.q)).collect();
        Poly {
            coeffs,
            q: self.q,
        }
    }
}

impl Neg for Poly {
    type Output = Poly;
    fn neg(self) -> Self::Output {
        -&self
    }
}

impl fmt::Debug for Poly {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "Poly(N={}, q={}, norm_inf={})",
            self.len(),
            self.q,
            self.infinity_norm()
        )
    }
}

/// A polynomial in the NTT (evaluation) domain.
#[derive(Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct NttPoly {
    pub evals: Vec<i32>,
    pub q: i32,
}

impl NttPoly {
    /// Zero polynomial in NTT domain.
    pub fn zero(n: usize, q: i32) -> Self {
        Self {
            evals: vec![0; n],
            q,
        }
    }

    /// Pointwise multiplication in NTT domain: C_i = A_i * B_i mod q.
    pub fn pointwise_mul(&self, other: &Self) -> Self {
        assert_eq!(self.evals.len(), other.evals.len());
        assert_eq!(self.q, other.q);
        let evals: Vec<i32> = self
            .evals
            .iter()
            .zip(&other.evals)
            .map(|(&a, &b)| mod_mul(a, b, self.q))
            .collect();
        Self { evals, q: self.q }
    }

    /// Pointwise addition in NTT domain.
    pub fn add(&self, other: &Self) -> Self {
        assert_eq!(self.evals.len(), other.evals.len());
        assert_eq!(self.q, other.q);
        let evals: Vec<i32> = self
            .evals
            .iter()
            .zip(&other.evals)
            .map(|(&a, &b)| mod_add(a, b, self.q))
            .collect();
        Self { evals, q: self.q }
    }

    /// Pointwise subtraction in NTT domain.
    pub fn sub(&self, other: &Self) -> Self {
        assert_eq!(self.evals.len(), other.evals.len());
        assert_eq!(self.q, other.q);
        let evals: Vec<i32> = self
            .evals
            .iter()
            .zip(&other.evals)
            .map(|(&a, &b)| mod_sub(a, b, self.q))
            .collect();
        Self { evals, q: self.q }
    }

    /// Converts NTT polynomial back to coefficient domain using INTT.
    pub fn to_poly(&self, ctx: &NttContext) -> Poly {
        assert_eq!(self.evals.len(), ctx.n);
        assert_eq!(self.q, ctx.q);
        let mut coeffs = self.evals.clone();
        ctx.inverse(&mut coeffs);
        Poly {
            coeffs,
            q: self.q,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use rand::SeedableRng;
    use rand_chacha::ChaCha20Rng;

    #[test]
    fn test_poly_neg_and_zero() {
        let mut rng = ChaCha20Rng::seed_from_u64(999);
        let p = Poly::random_uniform(DEFAULT_N, DEFAULT_Q, &mut rng);
        let neg_p = -&p;
        let sum = &p + &neg_p;
        assert!(sum.is_zero());
    }

    #[test]
    fn test_poly_scalar_distributivity() {
        let mut rng = ChaCha20Rng::seed_from_u64(1010);
        let p = Poly::random_uniform(DEFAULT_N, DEFAULT_Q, &mut rng);
        let s1 = 17;
        let s2 = 42;
        let left = p.scalar_mul(s1 + s2);
        let right = &p.scalar_mul(s1) + &p.scalar_mul(s2);
        assert_eq!(left, right);
    }

    #[test]
    fn test_sparse_challenge_properties() {
        let mut rng = ChaCha20Rng::seed_from_u64(2020);
        let tau = 40;
        let p = Poly::random_sparse_challenge(DEFAULT_N, DEFAULT_Q, tau, &mut rng);
        let non_zero = p.coeffs.iter().filter(|&&c| c != 0).count();
        assert_eq!(non_zero, tau);
        for &c in &p.coeffs {
            let centered = center(c, DEFAULT_Q);
            assert!(centered == -1 || centered == 0 || centered == 1);
        }
    }
}
