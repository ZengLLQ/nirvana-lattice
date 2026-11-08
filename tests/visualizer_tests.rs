use nirvana_lattice::lwe::{LatticeSampler, MlwKeypair, MlwParams};
use nirvana_lattice::ring::{NttContext, DEFAULT_N, DEFAULT_Q};
use nirvana_lattice::snark::{LatticeProver, LatticeVerifier, SnarkParams};
use nirvana_lattice::visualizer::AsciiVisualizer;
use rand::SeedableRng;
use rand_chacha::ChaCha20Rng;

#[test]
fn test_all_visualizations_are_strictly_7bit_ascii() {
    let ctx = NttContext::default_params();
    let mlw_params = MlwParams::standard();
    let snark_params = SnarkParams::fast();
    let mut rng = ChaCha20Rng::seed_from_u64(777);

    // 1. Butterfly diagram
    let bfly = AsciiVisualizer::render_ntt_butterfly_diagram();
    assert!(
        bfly.bytes().all(|b| b <= 127),
        "Butterfly diagram must be 100% 7-bit ASCII"
    );

    // 2. Gaussian histogram
    let samples: Vec<i32> = (0..500)
        .map(|_| LatticeSampler::sample_gaussian_int(3.0, &mut rng))
        .collect();
    let hist = AsciiVisualizer::render_gaussian_histogram(&samples, 10, 30);
    assert!(
        hist.bytes().all(|b| b <= 127),
        "Gaussian histogram must be 100% 7-bit ASCII"
    );

    // 3. Verification card
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
    .expect("Proof");
    let result = LatticeVerifier::verify(
        &keypair.public_matrix_a,
        &keypair.public_key_b,
        &proof,
        &snark_params,
        &ctx,
    );
    let card = AsciiVisualizer::render_verification_card(&proof, &result, &snark_params, 1.25);
    assert!(
        card.bytes().all(|b| b <= 127),
        "Verification card must be 100% 7-bit ASCII"
    );

    // 4. Polynomial spectrum
    let poly = LatticeSampler::sample_gaussian_poly(DEFAULT_N, DEFAULT_Q, 4.0, &mut rng);
    let spec = AsciiVisualizer::render_poly_spectrum(&poly.coeffs, DEFAULT_Q, 40);
    assert!(
        spec.bytes().all(|b| b <= 127),
        "Poly spectrum must be 100% 7-bit ASCII"
    );
}

#[test]
fn test_verification_card_contents() {
    let ctx = NttContext::default_params();
    let mlw_params = MlwParams::standard();
    let snark_params = SnarkParams::fast();
    let mut rng = ChaCha20Rng::seed_from_u64(888);

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
    .expect("Proof");
    let result = LatticeVerifier::verify(
        &keypair.public_matrix_a,
        &keypair.public_key_b,
        &proof,
        &snark_params,
        &ctx,
    );

    let card = AsciiVisualizer::render_verification_card(&proof, &result, &snark_params, 2.10);
    assert!(card.contains("[PASSED] VALID ARGUMENT OF KNOWLEDGE"));
    assert!(card.contains("NIRVANA-LATTICE POST-QUANTUM ZK-SNARK CERTIFICATE"));
    assert!(card.contains("Response Norm"));
}
