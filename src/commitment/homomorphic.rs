//! Lattice-based homomorphic commitment scheme over R_q = Z_q[X] / (X^N + 1).
//!
//! Implements Com(m; r) = A1 * r + A2 * m (mod q).
//! Supports additive and scalar polynomial homomorphism:
//! - Com(m1; r1) + Com(m2; r2) = Com(m1 + m2; r1 + r2)
//! - c * Com(m; r) = Com(c * m; c * r)

use std::ops::{Add, Sub};
use rand::Rng;
use serde::{Deserialize, Serialize};

use crate::lwe::matrix::{PolyMatrix, PolyVector};
use crate::lwe::sampler::LatticeSampler;
use crate::ring::ntt::NttContext;
use crate::ring::poly::Poly;
use crate::ring::reduction::{DEFAULT_N, DEFAULT_Q};

/// Public setup parameters for lattice homomorphic commitment scheme.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct CommitmentParams {
    pub k: usize,      // Output commitment dimension (rows)
    pub l_r: usize,    // Randomness dimension (cols of A1)
    pub l_m: usize,    // Message dimension (cols of A2)
    pub n: usize,      // Polynomial degree N = 256
    pub q: i32,        // Modulus q
    pub matrix_a1: PolyMatrix, // Binding matrix for randomness: k x l_r
    pub matrix_a2: PolyMatrix, // Projection matrix for message: k x l_m
}

impl CommitmentParams {
    /// Generates public commitment parameters uniformly at random.
    pub fn setup<R: Rng>(
        k: usize,
        l_r: usize,
        l_m: usize,
        n: usize,
        q: i32,
        rng: &mut R,
    ) -> Self {
        let matrix_a1 = LatticeSampler::sample_uniform_matrix(k, l_r, n, q, rng);
        let matrix_a2 = LatticeSampler::sample_uniform_matrix(k, l_m, n, q, rng);
        Self {
            k,
            l_r,
            l_m,
            n,
            q,
            matrix_a1,
            matrix_a2,
        }
    }

    /// Default parameters for Module-LWE zk-SNARK commitments (k=3, l_r=2, l_m=2, N=256, q=8380417).
    pub fn default_params<R: Rng>(rng: &mut R) -> Self {
        Self::setup(3, 2, 2, DEFAULT_N, DEFAULT_Q, rng)
    }

    /// Computes homomorphic commitment: C = A1 * r + A2 * m (mod q).
    pub fn commit(
        &self,
        message: &PolyVector,
        randomness: &PolyVector,
        ctx: &NttContext,
    ) -> Commitment {
        assert_eq!(message.dim(), self.l_m, "Message dimension must match l_m");
        assert_eq!(randomness.dim(), self.l_r, "Randomness dimension must match l_r");

        let a1_r = self.matrix_a1.mul_vector(randomness, ctx);
        let a2_m = self.matrix_a2.mul_vector(message, ctx);
        let c_vec = &a1_r + &a2_m;

        Commitment {
            vector: c_vec,
            q: self.q,
        }
    }

    /// Verifies opening of a commitment Com(m; r).
    pub fn verify_open(
        &self,
        commitment: &Commitment,
        message: &PolyVector,
        randomness: &PolyVector,
        beta_r: i32,
        ctx: &NttContext,
    ) -> bool {
        // 1. Norm check on randomness: ||r||_infinity <= beta_r
        if randomness.infinity_norm() > beta_r {
            return false;
        }

        // 2. Exact algebraic verification: C == A1 * r + A2 * m (mod q)
        let expected = self.commit(message, randomness, ctx);
        expected == *commitment
    }
}

/// A commitment value in R_q^k.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct Commitment {
    pub vector: PolyVector,
    pub q: i32,
}

impl Commitment {
    /// Dimension of the commitment vector.
    pub fn dim(&self) -> usize {
        self.vector.dim()
    }

    /// Infinity norm of commitment vector.
    pub fn infinity_norm(&self) -> i32 {
        self.vector.infinity_norm()
    }

    /// Scalar polynomial multiplication on commitment: c * Com(m; r)
    pub fn poly_scalar_mul(&self, scalar_poly: &Poly, ctx: &NttContext) -> Self {
        Self {
            vector: self.vector.poly_scalar_mul(scalar_poly, ctx),
            q: self.q,
        }
    }
}

impl Add for &Commitment {
    type Output = Commitment;
    fn add(self, rhs: Self) -> Self::Output {
        assert_eq!(self.q, rhs.q);
        Commitment {
            vector: &self.vector + &rhs.vector,
            q: self.q,
        }
    }
}

impl Add for Commitment {
    type Output = Commitment;
    fn add(self, rhs: Self) -> Self::Output {
        &self + &rhs
    }
}

impl Sub for &Commitment {
    type Output = Commitment;
    fn sub(self, rhs: Self) -> Self::Output {
        assert_eq!(self.q, rhs.q);
        Commitment {
            vector: &self.vector - &rhs.vector,
            q: self.q,
        }
    }
}

impl Sub for Commitment {
    type Output = Commitment;
    fn sub(self, rhs: Self) -> Self::Output {
        &self - &rhs
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use rand::SeedableRng;
    use rand_chacha::ChaCha20Rng;

    #[test]
    fn test_verify_open_success_and_failures() {
        let ctx = NttContext::default_params();
        let mut rng = ChaCha20Rng::seed_from_u64(9876);
        let params = CommitmentParams::default_params(&mut rng);

        let m = LatticeSampler::sample_gaussian_vector(params.l_m, DEFAULT_N, DEFAULT_Q, 2.0, &mut rng);
        let r = LatticeSampler::sample_gaussian_vector(params.l_r, DEFAULT_N, DEFAULT_Q, 2.0, &mut rng);
        let c = params.commit(&m, &r, &ctx);

        // Valid opening with bound 50
        assert!(params.verify_open(&c, &m, &r, 50, &ctx));

        // Invalid: exceeded bound
        assert!(!params.verify_open(&c, &m, &r, 1, &ctx));

        // Invalid: modified message
        let mut altered_m = m.clone();
        altered_m.polys[0].coeffs[0] = (altered_m.polys[0].coeffs[0] + 1) % DEFAULT_Q;
        assert!(!params.verify_open(&c, &altered_m, &r, 50, &ctx));
    }
}
