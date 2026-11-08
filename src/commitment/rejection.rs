//! Bimodal and uniform rejection sampling for lattice zero-knowledge proofs.
//!
//! Based on Vadim Lyubashevsky's technique (2009/2012) to decouple the response
//! distribution from the secret witness, ensuring statistical zero-knowledge.

use rand::Rng;
use serde::{Deserialize, Serialize};

use crate::lwe::matrix::PolyVector;

/// Diagnostic statistics from a rejection sampling run.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct RejectionStats {
    pub total_attempts: usize,
    pub accepted: bool,
    pub final_norm: i32,
    pub bound_threshold: i32,
    pub acceptance_ratio: f64,
}

/// Lyubashevsky Rejection Sampler engine.
pub struct RejectionSampler;

impl RejectionSampler {
    /// Evaluates the uniform rejection sampling condition for response `z = y + c * s`.
    ///
    /// Accepts if ||z||_infinity <= gamma - beta.
    /// When this condition holds, the distribution of `z` is strictly uniform
    /// in [-(gamma - beta), gamma - beta] and completely independent of secret `s`.
    pub fn check_uniform_rejection(z: &PolyVector, gamma: i32, beta: i32) -> bool {
        assert!(gamma > beta, "Masking bound gamma must be greater than secret bound beta");
        let threshold = gamma - beta;
        z.infinity_norm() <= threshold
    }

    /// Bimodal Gaussian rejection sampling acceptance check.
    ///
    /// Computes the bimodal acceptance probability:
    /// P_acc = 1 / (M * cosh(<z, cs> / sigma^2)) * exp(-||cs||^2 / (2 * sigma^2))
    /// where M is the rejection sampling slack parameter (typically M ~ 3).
    pub fn check_bimodal_gaussian<R: Rng>(
        z: &PolyVector,
        cs: &PolyVector,
        sigma: f64,
        m_slack: f64,
        rng: &mut R,
    ) -> bool {
        assert!(sigma > 0.0, "Sigma must be strictly positive");
        assert!(m_slack >= 1.0, "Slack parameter M must be >= 1.0");

        // Compute inner product <z, cs> across centered coefficients
        let mut inner_product: f64 = 0.0;
        let mut cs_norm_sq: f64 = 0.0;

        for (pz, pcs) in z.polys.iter().zip(&cs.polys) {
            let cz = pz.centered_coeffs();
            let ccs = pcs.centered_coeffs();
            for (&vz, &vcs) in cz.iter().zip(&ccs) {
                inner_product += (vz as f64) * (vcs as f64);
                cs_norm_sq += (vcs as f64) * (vcs as f64);
            }
        }

        let sigma_sq = sigma * sigma;
        let exponent = -cs_norm_sq / (2.0 * sigma_sq);
        let cosh_term = (inner_product / sigma_sq).cosh();

        if cosh_term <= 0.0 || cosh_term.is_nan() {
            return false;
        }

        let p_acc = (1.0 / (m_slack * cosh_term)) * exponent.exp();

        // Sample uniform u in [0, 1]
        let u: f64 = rng.gen();
        u < p_acc
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::ring::poly::Poly;
    use crate::ring::reduction::{DEFAULT_N, DEFAULT_Q};
    use rand::SeedableRng;
    use rand_chacha::ChaCha20Rng;

    #[test]
    fn test_bimodal_rejection_sampling() {
        let mut rng = ChaCha20Rng::seed_from_u64(4321);
        let z = PolyVector::new(vec![Poly::zero(DEFAULT_N, DEFAULT_Q)]);
        let cs = PolyVector::new(vec![Poly::zero(DEFAULT_N, DEFAULT_Q)]);

        // For z=0, cs=0: inner_prod=0, norm_sq=0 => p_acc = 1 / (M * cosh(0)) = 1 / M
        let mut accepts = 0;
        let trials = 1000;
        let m_slack = 2.0;
        for _ in 0..trials {
            if RejectionSampler::check_bimodal_gaussian(&z, &cs, 10.0, m_slack, &mut rng) {
                accepts += 1;
            }
        }
        let empirical_p = accepts as f64 / trials as f64;
        assert!((empirical_p - 0.5).abs() < 0.1, "Expected ~0.5 acceptance, got {}", empirical_p);
    }
}
