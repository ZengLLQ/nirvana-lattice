use nirvana_lattice::lwe::{MlwKeypair, MlwParams};
use nirvana_lattice::ring::{NttContext, DEFAULT_Q};
use nirvana_lattice::snark::{
    FiatShamirOracle, LatticeProof, LatticeProver, LatticeVerifier, SnarkParams, VerificationResult,
};
use rand::SeedableRng;
use rand_chacha::ChaCha20Rng;
use std::time::Instant;

#[test]
fn test_fiat_shamir_deterministic_and_unique() {
    let ctx = NttContext::default_params();
    let params = MlwParams::standard();
    let mut rng = ChaCha20Rng::seed_from_u64(1001);

    let keypair = MlwKeypair::generate(&params, &ctx, &mut rng);
    let w = keypair.public_key_b.clone();

    let c1 = FiatShamirOracle::compute_challenge(&keypair.public_matrix_a, &keypair.public_key_b, &w, 32);
    let c2 = FiatShamirOracle::compute_challenge(&keypair.public_matrix_a, &keypair.public_key_b, &w, 32);

    assert_eq!(c1, c2, "Fiat-Shamir challenge must be completely deterministic for identical inputs");

    // Count non-zero coefficients
    let non_zero_count = c1.coeffs.iter().filter(|&&c| c != 0).count();
    assert_eq!(non_zero_count, 32, "Challenge must have exactly tau=32 non-zero coefficients");
}

#[test]
fn test_zk_snark_prove_and_verify_valid() {
    let ctx = NttContext::default_params();
    let mlw_params = MlwParams::standard();
    let snark_params = SnarkParams::fast();
    let mut rng = ChaCha20Rng::seed_from_u64(2002);

    let keypair = MlwKeypair::generate(&mlw_params, &ctx, &mut rng);

    let start_prove = Instant::now();
    let proof = LatticeProver::prove(
        &keypair.public_matrix_a,
        &keypair.public_key_b,
        &keypair.secret_key_s,
        &keypair.error_e,
        &snark_params,
        &ctx,
        &mut rng,
    )
    .expect("Prover should generate valid proof");
    let prove_duration = start_prove.elapsed();

    println!("Proof generated in {:?}", prove_duration);
    println!("Proof required {} attempts", proof.attempts);

    let start_verify = Instant::now();
    let res = LatticeVerifier::verify(
        &keypair.public_matrix_a,
        &keypair.public_key_b,
        &proof,
        &snark_params,
        &ctx,
    );
    let verify_duration = start_verify.elapsed();

    println!("Verification finished in {:?}", verify_duration);
    assert_eq!(res, VerificationResult::Valid, "Honest proof must be verified as VALID");
    assert!(verify_duration.as_millis() < 50, "Verification should be sub-50ms");
}

#[test]
fn test_zk_snark_rejects_altered_response() {
    let ctx = NttContext::default_params();
    let mlw_params = MlwParams::standard();
    let snark_params = SnarkParams::fast();
    let mut rng = ChaCha20Rng::seed_from_u64(3003);

    let keypair = MlwKeypair::generate(&mlw_params, &ctx, &mut rng);

    let mut proof = LatticeProver::prove(
        &keypair.public_matrix_a,
        &keypair.public_key_b,
        &keypair.secret_key_s,
        &keypair.error_e,
        &snark_params,
        &ctx,
        &mut rng,
    )
    .expect("Proof creation");

    // Tamper with response z_s
    proof.response_z_s[0].coeffs[0] = (proof.response_z_s[0].coeffs[0] + 1) % DEFAULT_Q;

    let res = LatticeVerifier::verify(
        &keypair.public_matrix_a,
        &keypair.public_key_b,
        &proof,
        &snark_params,
        &ctx,
    );

    assert_ne!(res, VerificationResult::Valid, "Tampered response must be REJECTED");
}

#[test]
fn test_zk_snark_rejects_altered_challenge() {
    let ctx = NttContext::default_params();
    let mlw_params = MlwParams::standard();
    let snark_params = SnarkParams::fast();
    let mut rng = ChaCha20Rng::seed_from_u64(4004);

    let keypair = MlwKeypair::generate(&mlw_params, &ctx, &mut rng);

    let mut proof = LatticeProver::prove(
        &keypair.public_matrix_a,
        &keypair.public_key_b,
        &keypair.secret_key_s,
        &keypair.error_e,
        &snark_params,
        &ctx,
        &mut rng,
    )
    .expect("Proof creation");

    // Tamper with challenge c
    proof.challenge_c.coeffs[0] = (proof.challenge_c.coeffs[0] + 1) % DEFAULT_Q;

    let res = LatticeVerifier::verify(
        &keypair.public_matrix_a,
        &keypair.public_key_b,
        &proof,
        &snark_params,
        &ctx,
    );

    assert_ne!(res, VerificationResult::Valid, "Tampered challenge must be REJECTED");
}

#[test]
fn test_zk_snark_rejects_norm_overflow() {
    let ctx = NttContext::default_params();
    let mlw_params = MlwParams::standard();
    let snark_params = SnarkParams::fast();
    let mut rng = ChaCha20Rng::seed_from_u64(5005);

    let keypair = MlwKeypair::generate(&mlw_params, &ctx, &mut rng);

    let mut proof = LatticeProver::prove(
        &keypair.public_matrix_a,
        &keypair.public_key_b,
        &keypair.secret_key_s,
        &keypair.error_e,
        &snark_params,
        &ctx,
        &mut rng,
    )
    .expect("Proof creation");

    // Force norm overflow
    let huge_coeff = snark_params.threshold_s() + 500;
    proof.response_z_s[0].coeffs[0] = huge_coeff;

    let res = LatticeVerifier::verify(
        &keypair.public_matrix_a,
        &keypair.public_key_b,
        &proof,
        &snark_params,
        &ctx,
    );

    match res {
        VerificationResult::InvalidNorm { variable, observed, .. } => {
            assert_eq!(variable, "z_s");
            assert!(observed >= huge_coeff);
        }
        other => panic!("Expected InvalidNorm, got {:?}", other),
    }
}

#[test]
fn test_proof_serialization_roundtrip() {
    let ctx = NttContext::default_params();
    let mlw_params = MlwParams::standard();
    let snark_params = SnarkParams::fast();
    let mut rng = ChaCha20Rng::seed_from_u64(6006);

    let keypair = MlwKeypair::generate(&mlw_params, &ctx, &mut rng);

    let proof = LatticeProver::prove(
        &keypair.public_matrix_a,
        &keypair.public_key_b,
        &keypair.secret_key_s,
        &keypair.error_e,
        &snark_params,
        &ctx,
        &mut rng,
    )
    .expect("Proof creation");

    // JSON roundtrip
    let json_str = proof.to_json().expect("to_json");
    let decoded_json = LatticeProof::from_json(&json_str).expect("from_json");
    assert_eq!(proof, decoded_json);

    // Binary roundtrip
    let bytes = proof.to_bytes().expect("to_bytes");
    let decoded_bytes = LatticeProof::from_bytes(&bytes).expect("from_bytes");
    assert_eq!(proof, decoded_bytes);

    // Verify decoded proof
    let res = LatticeVerifier::verify(
        &keypair.public_matrix_a,
        &keypair.public_key_b,
        &decoded_bytes,
        &snark_params,
        &ctx,
    );
    assert_eq!(res, VerificationResult::Valid);
}
