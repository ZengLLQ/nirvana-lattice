//! Lattice zk-SNARK prover implementation with Lyubashevsky rejection sampling.

use rand::Rng;

use crate::lwe::matrix::{PolyMatrix, PolyVector};
use crate::lwe::sampler::LatticeSampler;
use crate::ring::ntt::NttContext;
use crate::snark::fiat_shamir::FiatShamirOracle;
use crate::snark::proof::LatticeProof;

/// Configuration parameters for the lattice zk-SNARK protocol.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SnarkParams {
    pub gamma_1: i32,    // Masking bound for secret witness s
    pub gamma_2: i32,    // Masking bound for error e
    pub tau: usize,      // Number of non-zero coefficients in challenge polynomial c
    pub beta_s: i32,     // Bound on secret key norm ||s||_infinity
    pub beta_e: i32,     // Bound on error norm ||e||_infinity
    pub max_attempts: usize, // Max rejection sampling attempts
}

impl SnarkParams {
    /// Standard secure parameters tuned for fast zero-knowledge argument generation.
    /// Compatible with CRYSTALS-Dilithium gamma_1 = (q - 1) / 16 = 523776.
    pub fn standard() -> Self {
        Self {
            gamma_1: 523776,
            gamma_2: 523776,
            tau: 32,
            beta_s: 6,
            beta_e: 6,
            max_attempts: 128,
        }
    }

    /// Fast interactive/test parameters with high acceptance probability (>70% per attempt).
    pub fn fast() -> Self {
        Self {
            gamma_1: 523776,
            gamma_2: 523776,
            tau: 20,
            beta_s: 5,
            beta_e: 5,
            max_attempts: 64,
        }
    }

    /// Response norm threshold for secret witness z_s.
    pub fn threshold_s(&self) -> i32 {
        self.gamma_1 - (self.beta_s * self.tau as i32)
    }

    /// Response norm threshold for error vector z_e.
    pub fn threshold_e(&self) -> i32 {
        self.gamma_2 - (self.beta_e * self.tau as i32)
    }
}

/// Lattice zk-SNARK Prover engine.
pub struct LatticeProver;

impl LatticeProver {
    /// Constructs a non-interactive zero-knowledge argument of knowledge proving:
    /// Knowledge of secret `s` and error `e` such that `b = A * s + e (mod q)`.
    pub fn prove<R: Rng>(
        matrix_a: &PolyMatrix,
        vector_b: &PolyVector,
        secret_s: &PolyVector,
        error_e: &PolyVector,
        params: &SnarkParams,
        ctx: &NttContext,
        rng: &mut R,
    ) -> Result<LatticeProof, String> {
        let n = ctx.n;
        let q = ctx.q;
        let k = matrix_a.rows;
        let l = matrix_a.cols;

        assert_eq!(secret_s.dim(), l, "Secret dimension must match matrix columns");
        assert_eq!(error_e.dim(), k, "Error dimension must match matrix rows");
        assert_eq!(vector_b.dim(), k, "Public vector dimension must match matrix rows");

        let thresh_s = params.threshold_s();
        let thresh_e = params.threshold_e();

        assert!(
            thresh_s > 0 && thresh_e > 0,
            "Rejection bounds must be positive (increase gamma or decrease tau * beta)"
        );

        for attempt in 1..=params.max_attempts {
            // Step 1: Sample uniform masking vectors y_s in [-gamma_1, gamma_1]^l, y_e in [-gamma_2, gamma_2]^k
            let y_s = LatticeSampler::sample_mask_vector(l, n, q, params.gamma_1, rng);
            let y_e = LatticeSampler::sample_mask_vector(k, n, q, params.gamma_2, rng);

            // Step 2: Compute commitment w = A * y_s + y_e (mod q)
            let a_ys = matrix_a.mul_vector(&y_s, ctx);
            let w = &a_ys + &y_e;

            // Step 3: Compute challenge polynomial c = H(A, b, w) via Fiat-Shamir
            let c = FiatShamirOracle::compute_challenge(matrix_a, vector_b, &w, params.tau);

            // Step 4: Compute candidate responses: z_s = y_s + c * s, z_e = y_e + c * e
            let cs = secret_s.poly_scalar_mul(&c, ctx);
            let ce = error_e.poly_scalar_mul(&c, ctx);

            let z_s = &y_s + &cs;
            let z_e = &y_e + &ce;

            // Step 5: Lyubashevsky's rejection sampling check
            if z_s.infinity_norm() <= thresh_s && z_e.infinity_norm() <= thresh_e {
                return Ok(LatticeProof {
                    commitment_w: w,
                    challenge_c: c,
                    response_z_s: z_s,
                    response_z_e: z_e,
                    attempts: attempt,
                });
            }
        }

        Err(format!(
            "Rejection sampling failed to converge within {} iterations",
            params.max_attempts
        ))
    }
}
