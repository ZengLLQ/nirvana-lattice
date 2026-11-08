//! Module-LWE key generation and M-SIS hardness verification.

use rand::Rng;
use serde::{Deserialize, Serialize};

use crate::lwe::matrix::{PolyMatrix, PolyVector};
use crate::lwe::sampler::LatticeSampler;
use crate::ring::ntt::NttContext;
use crate::ring::reduction::{DEFAULT_N, DEFAULT_Q};

/// Configuration parameters for Module-LWE scheme.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct MlwParams {
    pub k: usize,        // Number of rows in matrix A (dimension of b and e)
    pub l: usize,        // Number of columns in matrix A (dimension of secret s)
    pub n: usize,        // Polynomial degree (e.g. 256)
    pub q: i32,          // Modulus (e.g. 8380417)
    pub sigma: f64,      // Gaussian noise standard deviation for secret and error
    pub beta_s: i32,     // Bound on secret key norm ||s||_infinity
}

impl MlwParams {
    /// Recommended standard Module-LWE parameters (k=3, l=2, N=256, q=8380417, sigma=3.0).
    pub fn standard() -> Self {
        Self {
            k: 3,
            l: 2,
            n: DEFAULT_N,
            q: DEFAULT_Q,
            sigma: 3.0,
            beta_s: 15,
        }
    }

    /// Compact Ring-LWE parameter set (k=1, l=1, N=256, q=8380417).
    pub fn ring_lwe() -> Self {
        Self {
            k: 1,
            l: 1,
            n: DEFAULT_N,
            q: DEFAULT_Q,
            sigma: 3.0,
            beta_s: 15,
        }
    }
}

/// A Module-LWE keypair comprising public parameters and secret witness.
#[derive(Clone, Serialize, Deserialize)]
pub struct MlwKeypair {
    pub public_matrix_a: PolyMatrix,
    pub public_key_b: PolyVector,
    pub secret_key_s: PolyVector,
    pub error_e: PolyVector,
    pub params: MlwParams,
}

impl MlwKeypair {
    /// Generates a new Module-LWE keypair: b = A * s + e (mod q).
    pub fn generate<R: Rng>(params: &MlwParams, ctx: &NttContext, rng: &mut R) -> Self {
        assert_eq!(params.n, ctx.n);
        assert_eq!(params.q, ctx.q);

        // Sample public matrix A uniformly in R_q^(k x l)
        let public_matrix_a = LatticeSampler::sample_uniform_matrix(
            params.k,
            params.l,
            params.n,
            params.q,
            rng,
        );

        // Sample secret s from discrete Gaussian in R_q^l
        let secret_key_s = LatticeSampler::sample_gaussian_vector(
            params.l,
            params.n,
            params.q,
            params.sigma,
            rng,
        );

        // Sample noise error e from discrete Gaussian in R_q^k
        let error_e = LatticeSampler::sample_gaussian_vector(
            params.k,
            params.n,
            params.q,
            params.sigma,
            rng,
        );

        // Compute b = A * s + e mod q
        let as_prod = public_matrix_a.mul_vector(&secret_key_s, ctx);
        let public_key_b = &as_prod + &error_e;

        Self {
            public_matrix_a,
            public_key_b,
            secret_key_s,
            error_e,
            params: *params,
        }
    }

    /// Checks whether the keypair relationship b = A * s + e (mod q) strictly holds.
    pub fn verify_relation(&self, ctx: &NttContext) -> bool {
        let as_prod = self.public_matrix_a.mul_vector(&self.secret_key_s, ctx);
        let computed_b = &as_prod + &self.error_e;
        computed_b == self.public_key_b
    }
}

/// Module Short Integer Solution (M-SIS) problem hardness verifier.
/// Verifies if candidate vector `v` is a non-trivial solution to A * v = 0 (mod q)
/// with bounded infinity norm ||v||_infinity <= beta_sis.
pub fn verify_msis(
    matrix_a: &PolyMatrix,
    candidate_v: &PolyVector,
    beta_sis: i32,
    ctx: &NttContext,
) -> bool {
    // 1. Non-zero condition
    if candidate_v.is_zero() {
        return false;
    }

    // 2. Norm bound condition
    if candidate_v.infinity_norm() > beta_sis {
        return false;
    }

    // 3. A * v == 0 (mod q)
    let av = matrix_a.mul_vector(candidate_v, ctx);
    av.is_zero()
}
