use nirvana_lattice::ring::{
    barrett_reduce, center, freeze, mod_add, mod_inv, mod_mul, mod_pow, mod_sub,
    BarrettReducer, MontgomeryReducer, NttContext, Poly, DEFAULT_N, DEFAULT_Q,
};
use rand::SeedableRng;
use rand_chacha::ChaCha20Rng;

#[test]
fn test_modular_arithmetic_basic() {
    let q = DEFAULT_Q;
    assert_eq!(mod_add(10, 20, q), 30);
    assert_eq!(mod_sub(20, 30, q), q - 10);
    assert_eq!(mod_mul(1000, 2000, q), 2000000 % q);
    assert_eq!(freeze(-5, q), q - 5);
    assert_eq!(freeze((q as i64) * 2 + 7, q), 7);
}

#[test]
fn test_center_balanced() {
    let q = 17;
    assert_eq!(center(0, q), 0);
    assert_eq!(center(8, q), 8);
    assert_eq!(center(9, q), -8);
    assert_eq!(center(16, q), -1);
}

#[test]
fn test_mod_inv_and_pow() {
    let q = DEFAULT_Q;
    for a in [2, 3, 5, 17, 256, 123456, q - 2] {
        let inv = mod_inv(a, q).expect("inverse should exist");
        assert_eq!(mod_mul(a, inv, q), 1);
    }
    assert_eq!(mod_pow(3, 0, q), 1);
    assert_eq!(mod_pow(3, 1, q), 3);
    assert_eq!(mod_pow(2, 10, q), 1024);
}

#[test]
fn test_barrett_and_montgomery_reducers() {
    let q = DEFAULT_Q;
    let b_reducer = BarrettReducer::new(q);
    let m_reducer = MontgomeryReducer::new(q);

    for val in [0i64, 1, 42, 1000000, (q as i64) + 50, (q as i64) * 5 + 999] {
        let expected = freeze(val, q);
        assert_eq!(b_reducer.reduce(val), expected);
        assert_eq!(barrett_reduce(val, q), expected);

        // Montgomery domain test
        let val_i32 = freeze(val, q);
        let m_form = m_reducer.to_montgomery(val_i32);
        let back = m_reducer.from_montgomery(m_form);
        assert_eq!(back, val_i32);
    }
}

#[test]
fn test_ntt_roundtrip_identity() {
    let ctx = NttContext::default_params();
    let mut rng = ChaCha20Rng::seed_from_u64(42);

    for _ in 0..10 {
        let p = Poly::random_uniform(DEFAULT_N, DEFAULT_Q, &mut rng);
        let mut coeffs = p.coeffs.clone();

        ctx.forward(&mut coeffs);
        // After forward NTT, coefficients are in evaluation domain
        ctx.inverse(&mut coeffs);

        assert_eq!(coeffs, p.coeffs, "INTT(NTT(p)) must equal p");
    }
}

#[test]
fn test_ntt_multiplication_vs_schoolbook() {
    let ctx = NttContext::default_params();
    let mut rng = ChaCha20Rng::seed_from_u64(1337);

    for _ in 0..5 {
        let p1 = Poly::random_uniform(DEFAULT_N, DEFAULT_Q, &mut rng);
        let p2 = Poly::random_uniform(DEFAULT_N, DEFAULT_Q, &mut rng);

        let ntt_res = p1.ntt_mul(&p2, &ctx);
        let schoolbook_res = p1.schoolbook_mul(&p2);

        assert_eq!(
            ntt_res.coeffs, schoolbook_res.coeffs,
            "Negacyclic NTT multiplication must match schoolbook polynomial multiplication"
        );
    }
}
