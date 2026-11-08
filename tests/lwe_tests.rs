use nirvana_lattice::lwe::{
    verify_msis, LatticeSampler, MlwKeypair, MlwParams, PolyVector,
};
use nirvana_lattice::ring::{NttContext, Poly, DEFAULT_N, DEFAULT_Q};
use rand::SeedableRng;
use rand_chacha::ChaCha20Rng;

#[test]
fn test_poly_vector_ops() {
    let ctx = NttContext::default_params();
    let mut rng = ChaCha20Rng::seed_from_u64(100);

    let v1 = LatticeSampler::sample_uniform_vector(3, DEFAULT_N, DEFAULT_Q, &mut rng);
    let v2 = LatticeSampler::sample_uniform_vector(3, DEFAULT_N, DEFAULT_Q, &mut rng);

    let sum = &v1 + &v2;
    let diff = &sum - &v2;
    assert_eq!(diff, v1, "Vector (v1 + v2) - v2 must equal v1");

    // Scalar polynomial multiplication
    let scalar_poly = Poly::random_uniform(DEFAULT_N, DEFAULT_Q, &mut rng);
    let scaled = v1.poly_scalar_mul(&scalar_poly, &ctx);
    assert_eq!(scaled.dim(), 3);
}

#[test]
fn test_poly_matrix_vector_mul() {
    let ctx = NttContext::default_params();
    let mut rng = ChaCha20Rng::seed_from_u64(200);

    let rows = 3;
    let cols = 2;
    let mat = LatticeSampler::sample_uniform_matrix(rows, cols, DEFAULT_N, DEFAULT_Q, &mut rng);
    let vec = LatticeSampler::sample_uniform_vector(cols, DEFAULT_N, DEFAULT_Q, &mut rng);

    let res = mat.mul_vector(&vec, &ctx);
    assert_eq!(res.dim(), rows);

    // Distributivity: A * (v1 + v2) == A * v1 + A * v2
    let vec2 = LatticeSampler::sample_uniform_vector(cols, DEFAULT_N, DEFAULT_Q, &mut rng);
    let sum_vec = &vec + &vec2;
    let left = mat.mul_vector(&sum_vec, &ctx);
    let right = &mat.mul_vector(&vec, &ctx) + &mat.mul_vector(&vec2, &ctx);
    assert_eq!(left, right, "A * (v1 + v2) must equal A * v1 + A * v2");
}

#[test]
fn test_mlw_keypair_generation_and_relation() {
    let ctx = NttContext::default_params();
    let params = MlwParams::standard();
    let mut rng = ChaCha20Rng::seed_from_u64(300);

    let keypair = MlwKeypair::generate(&params, &ctx, &mut rng);
    assert!(
        keypair.verify_relation(&ctx),
        "b = A * s + e must hold strictly for generated Module-LWE keypair"
    );

    // Verify small norm of secret and error
    assert!(
        keypair.secret_key_s.infinity_norm() < 50,
        "Secret key norm must be small"
    );
    assert!(
        keypair.error_e.infinity_norm() < 50,
        "Error vector norm must be small"
    );
}

#[test]
fn test_gaussian_sampler_properties() {
    let mut rng = ChaCha20Rng::seed_from_u64(400);
    let sigma = 3.0;
    let n = 2048;
    let samples: Vec<i32> = (0..n)
        .map(|_| LatticeSampler::sample_gaussian_int(sigma, &mut rng))
        .collect();

    let mean: f64 = samples.iter().map(|&x| x as f64).sum::<f64>() / (n as f64);
    let variance: f64 = samples
        .iter()
        .map(|&x| {
            let diff = x as f64 - mean;
            diff * diff
        })
        .sum::<f64>()
        / (n as f64);

    assert!(mean.abs() < 0.2, "Sample mean should be close to 0: got {}", mean);
    assert!(
        (variance - sigma * sigma).abs() < 1.0,
        "Sample variance should be close to sigma^2 ({}), got {}",
        sigma * sigma,
        variance
    );
}

#[test]
fn test_msis_verifier() {
    let ctx = NttContext::default_params();
    let rows = 2;
    let cols = 3;
    let mut rng = ChaCha20Rng::seed_from_u64(500);

    let mat = LatticeSampler::sample_uniform_matrix(rows, cols, DEFAULT_N, DEFAULT_Q, &mut rng);
    let zero_v = PolyVector::zero(cols, DEFAULT_N, DEFAULT_Q);

    // Zero vector is not a non-trivial solution
    assert!(!verify_msis(&mat, &zero_v, 100, &ctx));

    // Non-zero random vector won't satisfy A * v == 0 with high probability
    let rand_v = LatticeSampler::sample_gaussian_vector(cols, DEFAULT_N, DEFAULT_Q, 2.0, &mut rng);
    assert!(!verify_msis(&mat, &rand_v, 100, &ctx));
}
