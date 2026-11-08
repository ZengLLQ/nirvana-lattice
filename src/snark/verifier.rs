//! Lattice zk-SNARK verifier implementation.

use crate::lwe::matrix::{PolyMatrix, PolyVector};
use crate::ring::ntt::NttContext;
use crate::snark::fiat_shamir::FiatShamirOracle;
use crate::snark::proof::{LatticeProof, VerificationResult};
use crate::snark::prover::SnarkParams;

/// Lattice zk-SNARK Verifier engine.
pub struct LatticeVerifier;

impl LatticeVerifier {
    /// Verifies the validity of a non-interactive zero-knowledge argument.
    ///
    /// Checks:
    /// 1. Response norm bounds ||z_s||_infinity <= threshold_s and ||z_e||_infinity <= threshold_e.
    /// 2. Algebraic relation: w' = A * z_s + z_e - c * b (mod q).
    /// 3. Fiat-Shamir challenge consistency: c == H(A, b, w').
    pub fn verify(
        matrix_a: &PolyMatrix,
        vector_b: &PolyVector,
        proof: &LatticeProof,
        params: &SnarkParams,
        ctx: &NttContext,
    ) -> VerificationResult {
        let k = matrix_a.rows;
        let l = matrix_a.cols;

        // 1. Dimension validation
        if proof.response_z_s.dim() != l {
            return VerificationResult::DimensionMismatch {
                expected: l,
                observed: proof.response_z_s.dim(),
            };
        }
        if proof.response_z_e.dim() != k {
            return VerificationResult::DimensionMismatch {
                expected: k,
                observed: proof.response_z_e.dim(),
            };
        }
        if proof.commitment_w.dim() != k {
            return VerificationResult::DimensionMismatch {
                expected: k,
                observed: proof.commitment_w.dim(),
            };
        }

        // 2. Response norm bounds verification (Crucial for zero-knowledge and SIS hardness!)
        let thresh_s = params.threshold_s();
        let norm_s = proof.response_z_s.infinity_norm();
        if norm_s > thresh_s {
            return VerificationResult::InvalidNorm {
                variable: "z_s".to_string(),
                observed: norm_s,
                threshold: thresh_s,
            };
        }

        let thresh_e = params.threshold_e();
        let norm_e = proof.response_z_e.infinity_norm();
        if norm_e > thresh_e {
            return VerificationResult::InvalidNorm {
                variable: "z_e".to_string(),
                observed: norm_e,
                threshold: thresh_e,
            };
        }

        // 3. Reconstruct commitment: w' = A * z_s + z_e - c * b (mod q)
        let a_zs = matrix_a.mul_vector(&proof.response_z_s, ctx);
        let cb = vector_b.poly_scalar_mul(&proof.challenge_c, ctx);
        let reconstructed_w = &(&a_zs + &proof.response_z_e) - &cb;

        // Verify that reconstructed w' matches commitment w
        if reconstructed_w != proof.commitment_w {
            return VerificationResult::ChallengeMismatch;
        }

        // 4. Fiat-Shamir challenge consistency check
        let expected_c =
            FiatShamirOracle::compute_challenge(matrix_a, vector_b, &reconstructed_w, params.tau);
        if expected_c != proof.challenge_c {
            return VerificationResult::ChallengeMismatch;
        }

        VerificationResult::Valid
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::lwe::{MlwKeypair, MlwParams};
    use crate::snark::prover::LatticeProver;
    use rand::SeedableRng;
    use rand_chacha::ChaCha20Rng;

    #[test]
    fn test_dimension_mismatch_detection() {
        let ctx = NttContext::default_params();
        let mut rng = ChaCha20Rng::seed_from_u64(6543);
        let params_standard = MlwParams::standard();
        let snark_params = SnarkParams::fast();

        let keypair = MlwKeypair::generate(&params_standard, &ctx, &mut rng);
        let proof = LatticeProver::prove(
            &keypair.public_matrix_a,
            &keypair.public_key_b,
            &keypair.secret_key_s,
            &keypair.error_e,
            &snark_params,
            &ctx,
            &mut rng,
        )
        .expect("Proof");

        // Create a mismatched matrix of different dimension (rows = 4 instead of 3)
        let mismatched_params = MlwParams {
            k: 4,
            l: 2,
            ..params_standard
        };
        let mismatched_keypair = MlwKeypair::generate(&mismatched_params, &ctx, &mut rng);

        let res = LatticeVerifier::verify(
            &mismatched_keypair.public_matrix_a,
            &mismatched_keypair.public_key_b,
            &proof,
            &snark_params,
            &ctx,
        );
        assert!(matches!(res, VerificationResult::DimensionMismatch { .. }));
    }
}
