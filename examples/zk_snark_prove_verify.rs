//! Example: Lattice-based zk-SNARK Prover and Verifier with Rejection Sampling.
//!
//! Run with: `cargo run --example zk_snark_prove_verify`

use std::time::Instant;
use nirvana_lattice::lwe::{MlwKeypair, MlwParams};
use nirvana_lattice::ring::{NttContext, DEFAULT_Q};
use nirvana_lattice::snark::{
    LatticeProver, LatticeVerifier, SnarkParams, VerificationResult,
};
use nirvana_lattice::visualizer::AsciiVisualizer;
use rand::thread_rng;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    println!("=========================================================================");
    println!("        NIRVANA-LATTICE: POST-QUANTUM ZK-SNARK PROVE & VERIFY DEMO       ");
    println!("=========================================================================");

    let ctx = NttContext::default_params();
    let mut rng = thread_rng();

    // 1. Setup instance
    println!("\n[1] Setting up Module-LWE Instance (k=3, l=2, N=256, q=8380417)...");
    let mlw_params = MlwParams::standard();
    let keypair = MlwKeypair::generate(&mlw_params, &ctx, &mut rng);
    println!("Public key b = A * s + e established.");

    // 2. Prover step
    println!("\n[2] Constructing Zero-Knowledge Succinct Argument of Knowledge...");
    let snark_params = SnarkParams::fast();
    let start_prove = Instant::now();
    let proof = LatticeProver::prove(
        &keypair.public_matrix_a,
        &keypair.public_key_b,
        &keypair.secret_key_s,
        &keypair.error_e,
        &snark_params,
        &ctx,
        &mut rng,
    )?;
    let prove_elapsed = start_prove.elapsed();

    println!("Proof generated successfully in {:.3} ms", prove_elapsed.as_secs_f64() * 1000.0);
    println!("Rejection sampling attempts:   {}", proof.attempts);
    println!("Proof JSON serialized size:    {} bytes", proof.to_json()?.len());
    println!("Proof binary serialized size:  {} bytes", proof.size_bytes());

    // 3. Verifier step
    println!("\n[3] Verifying zk-SNARK Argument with Verifier...");
    let start_verify = Instant::now();
    let result = LatticeVerifier::verify(
        &keypair.public_matrix_a,
        &keypair.public_key_b,
        &proof,
        &snark_params,
        &ctx,
    );
    let verify_elapsed = start_verify.elapsed();
    let verify_ms = verify_elapsed.as_secs_f64() * 1000.0;

    // 4. ASCII Certificate
    println!("\n[4] Rendering Verification Certificate Card:\n");
    let card = AsciiVisualizer::render_verification_card(&proof, &result, &snark_params, verify_ms);
    println!("{}", card);

    assert_eq!(result, VerificationResult::Valid, "Honest proof must be accepted");

    // 5. Negative Test (Tamper Resistance)
    println!("[5] Negative Security Test: Tampering with response vector...");
    let mut forged_proof = proof.clone();
    forged_proof.response_z_s[0].coeffs[5] = (forged_proof.response_z_s[0].coeffs[5] + 100) % DEFAULT_Q;

    let forged_result = LatticeVerifier::verify(
        &keypair.public_matrix_a,
        &keypair.public_key_b,
        &forged_proof,
        &snark_params,
        &ctx,
    );
    println!("Tampered proof verification result: {}", forged_result);
    assert_ne!(forged_result, VerificationResult::Valid, "Forged proof must be rejected!");
    println!("Soundness test passed: Forged arguments are strictly rejected!\n");

    Ok(())
}
