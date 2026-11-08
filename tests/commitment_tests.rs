use nirvana_lattice::commitment::{CommitmentParams, RejectionSampler};
use nirvana_lattice::lwe::LatticeSampler;
use nirvana_lattice::ring::{NttContext, Poly, DEFAULT_N, DEFAULT_Q};
use rand::SeedableRng;
use rand_chacha::ChaCha20Rng;

#[test]
fn test_homomorphic_commitment_additive() {
    let ctx = NttContext::default_params();
    let mut rng = ChaCha20Rng::seed_from_u64(111);

    let k = 3;
    let l_r = 2;
    let l_m = 2;
    let params = CommitmentParams::setup(k, l_r, l_m, DEFAULT_N, DEFAULT_Q, &mut rng);

    let m1 = LatticeSampler::sample_gaussian_vector(l_m, DEFAULT_N, DEFAULT_Q, 2.0, &mut rng);
    let r1 = LatticeSampler::sample_gaussian_vector(l_r, DEFAULT_N, DEFAULT_Q, 2.0, &mut rng);

    let m2 = LatticeSampler::sample_gaussian_vector(l_m, DEFAULT_N, DEFAULT_Q, 2.0, &mut rng);
    let r2 = LatticeSampler::sample_gaussian_vector(l_r, DEFAULT_N, DEFAULT_Q, 2.0, &mut rng);

    let c1 = params.commit(&m1, &r1, &ctx);
    let c2 = params.commit(&m2, &r2, &ctx);

    // Sum of commitments
    let c_sum = &c1 + &c2;

    // Commitment to sum of messages and randomness
    let m_sum = &m1 + &m2;
    let r_sum = &r1 + &r2;
    let c_expected = params.commit(&m_sum, &r_sum, &ctx);

    assert_eq!(
        c_sum, c_expected,
        "Homomorphic property Com(m1; r1) + Com(m2; r2) == Com(m1 + m2; r1 + r2) must hold"
    );
}

#[test]
fn test_homomorphic_commitment_subtractive() {
    let ctx = NttContext::default_params();
    let mut rng = ChaCha20Rng::seed_from_u64(222);

    let k = 3;
    let l_r = 2;
    let l_m = 2;
    let params = CommitmentParams::setup(k, l_r, l_m, DEFAULT_N, DEFAULT_Q, &mut rng);

    let m1 = LatticeSampler::sample_uniform_vector(l_m, DEFAULT_N, DEFAULT_Q, &mut rng);
    let r1 = LatticeSampler::sample_uniform_vector(l_r, DEFAULT_N, DEFAULT_Q, &mut rng);

    let m2 = LatticeSampler::sample_uniform_vector(l_m, DEFAULT_N, DEFAULT_Q, &mut rng);
    let r2 = LatticeSampler::sample_uniform_vector(l_r, DEFAULT_N, DEFAULT_Q, &mut rng);

    let c1 = params.commit(&m1, &r1, &ctx);
    let c2 = params.commit(&m2, &r2, &ctx);

    let c_diff = &c1 - &c2;
    let m_diff = &m1 - &m2;
    let r_diff = &r1 - &r2;
    let c_expected = params.commit(&m_diff, &r_diff, &ctx);

    assert_eq!(
        c_diff, c_expected,
        "Homomorphic property Com(m1; r1) - Com(m2; r2) == Com(m1 - m2; r1 - r2) must hold"
    );
}

#[test]
fn test_homomorphic_commitment_scalar_poly() {
    let ctx = NttContext::default_params();
    let mut rng = ChaCha20Rng::seed_from_u64(333);

    let k = 3;
    let l_r = 2;
    let l_m = 2;
    let params = CommitmentParams::setup(k, l_r, l_m, DEFAULT_N, DEFAULT_Q, &mut rng);

    let m = LatticeSampler::sample_gaussian_vector(l_m, DEFAULT_N, DEFAULT_Q, 2.0, &mut rng);
    let r = LatticeSampler::sample_gaussian_vector(l_r, DEFAULT_N, DEFAULT_Q, 2.0, &mut rng);
    let c = params.commit(&m, &r, &ctx);

    let challenge = Poly::random_sparse_challenge(DEFAULT_N, DEFAULT_Q, 60, &mut rng);

    let scaled_c = c.poly_scalar_mul(&challenge, &ctx);
    let scaled_m = m.poly_scalar_mul(&challenge, &ctx);
    let scaled_r = r.poly_scalar_mul(&challenge, &ctx);
    let expected_c = params.commit(&scaled_m, &scaled_r, &ctx);

    assert_eq!(
        scaled_c, expected_c,
        "Scalar homomorphic property c * Com(m; r) == Com(c * m; c * r) must hold"
    );
}

#[test]
fn test_rejection_sampling_uniform_bound() {
    let mut rng = ChaCha20Rng::seed_from_u64(444);
    let gamma = 1000;
    let beta = 100;

    // In-bound vector
    let z_in = LatticeSampler::sample_mask_vector(2, DEFAULT_N, DEFAULT_Q, 800, &mut rng);
    assert!(
        RejectionSampler::check_uniform_rejection(&z_in, gamma, beta),
        "Vector with norm <= 800 <= 900 must be accepted"
    );

    // Out-of-bound vector
    let z_out = LatticeSampler::sample_mask_vector(2, DEFAULT_N, DEFAULT_Q, 950, &mut rng);
    // Almost certainly some coefficient has absolute value > 900
    assert!(
        !RejectionSampler::check_uniform_rejection(&z_out, gamma, beta)
            || z_out.infinity_norm() <= 900
    );
}
