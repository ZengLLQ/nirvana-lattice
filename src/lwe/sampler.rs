//! Noise and distribution sampling for Ring-LWE and Module-LWE.
//!
//! Provides Discrete Gaussian sampling, Centered Binomial Distribution (CBD),
//! bounded uniform masking sampling, and uniform matrix/vector sampling.

use rand::Rng;
use rand_distr::{Distribution, Normal};

use crate::lwe::matrix::{PolyMatrix, PolyVector};
use crate::ring::poly::Poly;
use crate::ring::reduction::freeze;

/// Sampler for discrete distributions in lattice cryptography.
pub struct LatticeSampler;

impl LatticeSampler {
    /// Samples a single integer from discrete Gaussian distribution D_sigma centered at 0.
    pub fn sample_gaussian_int<R: Rng>(sigma: f64, rng: &mut R) -> i32 {
        let normal = Normal::new(0.0, sigma).expect("Valid Gaussian standard deviation");
        let val: f64 = normal.sample(rng);
        val.round() as i32
    }

    /// Samples a polynomial whose coefficients follow a discrete Gaussian distribution D_sigma.
    pub fn sample_gaussian_poly<R: Rng>(n: usize, q: i32, sigma: f64, rng: &mut R) -> Poly {
        let normal = Normal::new(0.0, sigma).expect("Valid Gaussian standard deviation");
        let coeffs: Vec<i32> = (0..n)
            .map(|_| {
                let s = normal.sample(rng).round() as i64;
                freeze(s, q)
            })
            .collect();
        Poly { coeffs, q }
    }

    /// Samples a vector of polynomials from discrete Gaussian distribution D_sigma.
    pub fn sample_gaussian_vector<R: Rng>(
        dim: usize,
        n: usize,
        q: i32,
        sigma: f64,
        rng: &mut R,
    ) -> PolyVector {
        let polys = (0..dim)
            .map(|_| Self::sample_gaussian_poly(n, q, sigma, rng))
            .collect();
        PolyVector { polys, q }
    }

    /// Samples a polynomial from Centered Binomial Distribution CBD_eta.
    /// Sum of eta Bernoulli trials: sum_{i=1}^eta a_i - sum_{i=1}^eta b_i where a_i, b_i in {0, 1}.
    pub fn sample_cbd_poly<R: Rng>(n: usize, q: i32, eta: usize, rng: &mut R) -> Poly {
        let coeffs: Vec<i32> = (0..n)
            .map(|_| {
                let mut a = 0i32;
                let mut b = 0i32;
                for _ in 0..eta {
                    if rng.gen_bool(0.5) {
                        a += 1;
                    }
                    if rng.gen_bool(0.5) {
                        b += 1;
                    }
                }
                let diff = (a - b) as i64;
                freeze(diff, q)
            })
            .collect();
        Poly { coeffs, q }
    }

    /// Samples a vector of polynomials from Centered Binomial Distribution.
    pub fn sample_cbd_vector<R: Rng>(
        dim: usize,
        n: usize,
        q: i32,
        eta: usize,
        rng: &mut R,
    ) -> PolyVector {
        let polys = (0..dim)
            .map(|_| Self::sample_cbd_poly(n, q, eta, rng))
            .collect();
        PolyVector { polys, q }
    }

    /// Uniformly samples a polynomial with coefficients in [0, q - 1].
    pub fn sample_uniform_poly<R: Rng>(n: usize, q: i32, rng: &mut R) -> Poly {
        Poly::random_uniform(n, q, rng)
    }

    /// Uniformly samples a vector of polynomials in R_q^dim.
    pub fn sample_uniform_vector<R: Rng>(dim: usize, n: usize, q: i32, rng: &mut R) -> PolyVector {
        let polys = (0..dim)
            .map(|_| Poly::random_uniform(n, q, rng))
            .collect();
        PolyVector { polys, q }
    }

    /// Uniformly samples a matrix of polynomials in R_q^(rows x cols).
    pub fn sample_uniform_matrix<R: Rng>(
        rows: usize,
        cols: usize,
        n: usize,
        q: i32,
        rng: &mut R,
    ) -> PolyMatrix {
        let data = (0..rows)
            .map(|_| {
                (0..cols)
                    .map(|_| Poly::random_uniform(n, q, rng))
                    .collect()
            })
            .collect();
        PolyMatrix { rows, cols, data, q }
    }

    /// Samples a bounded uniform masking polynomial with coefficients in [-gamma, gamma].
    /// Used in Lyubashevsky's rejection sampling framework for zero-knowledge arguments.
    pub fn sample_mask_poly<R: Rng>(n: usize, q: i32, gamma: i32, rng: &mut R) -> Poly {
        let coeffs: Vec<i32> = (0..n)
            .map(|_| {
                let v = rng.gen_range(-gamma..=gamma) as i64;
                freeze(v, q)
            })
            .collect();
        Poly { coeffs, q }
    }

    /// Samples a vector of bounded uniform masking polynomials.
    pub fn sample_mask_vector<R: Rng>(
        dim: usize,
        n: usize,
        q: i32,
        gamma: i32,
        rng: &mut R,
    ) -> PolyVector {
        let polys = (0..dim)
            .map(|_| Self::sample_mask_poly(n, q, gamma, rng))
            .collect();
        PolyVector { polys, q }
    }

    /// Samples a small ternary vector with coefficients in {-1, 0, 1}.
    pub fn sample_ternary_vector<R: Rng>(dim: usize, n: usize, q: i32, rng: &mut R) -> PolyVector {
        let polys = (0..dim)
            .map(|_| Poly::random_ternary(n, q, rng))
            .collect();
        PolyVector { polys, q }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::ring::reduction::{center, DEFAULT_N, DEFAULT_Q};
    use rand::SeedableRng;
    use rand_chacha::ChaCha20Rng;

    #[test]
    fn test_cbd_sampler_bounds() {
        let mut rng = ChaCha20Rng::seed_from_u64(1234);
        let eta = 2;
        let p = LatticeSampler::sample_cbd_poly(DEFAULT_N, DEFAULT_Q, eta, &mut rng);
        for &c in &p.coeffs {
            let centered = center(c, DEFAULT_Q);
            assert!(centered.abs() <= eta as i32);
        }
    }

    #[test]
    fn test_mask_sampler_bounds() {
        let mut rng = ChaCha20Rng::seed_from_u64(5678);
        let gamma = 5000;
        let v = LatticeSampler::sample_mask_vector(3, DEFAULT_N, DEFAULT_Q, gamma, &mut rng);
        assert!(v.infinity_norm() <= gamma);
    }

    #[test]
    fn test_ternary_vector_sampler() {
        let mut rng = ChaCha20Rng::seed_from_u64(9012);
        let v = LatticeSampler::sample_ternary_vector(2, DEFAULT_N, DEFAULT_Q, &mut rng);
        assert!(v.infinity_norm() <= 1);
    }
}
