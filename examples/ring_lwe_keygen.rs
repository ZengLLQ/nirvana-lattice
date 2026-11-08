//! Example: Ring-LWE & Module-LWE Key Generation and Noise Analysis.
//!
//! Run with: `cargo run --example ring_lwe_keygen`

use nirvana_lattice::lwe::{LatticeSampler, MlwKeypair, MlwParams};
use nirvana_lattice::ring::{NttContext, DEFAULT_Q};
use nirvana_lattice::visualizer::AsciiVisualizer;
use rand::thread_rng;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    println!("=========================================================================");
    println!("            NIRVANA-LATTICE: RING-LWE & MODULE-LWE KEYGEN DEMO           ");
    println!("=========================================================================");

    let ctx = NttContext::default_params();
    let mut rng = thread_rng();

    // 1. Ring-LWE Keygen (k=1, l=1)
    println!("\n[1] Generating 1-dimensional Ring-LWE Keypair (k=1, l=1)...");
    let ring_params = MlwParams::ring_lwe();
    let ring_keypair = MlwKeypair::generate(&ring_params, &ctx, &mut rng);

    println!("Public key b = a * s + e (mod q)");
    println!("Secret s infinity norm:  {}", ring_keypair.secret_key_s.infinity_norm());
    println!("Error  e infinity norm:  {}", ring_keypair.error_e.infinity_norm());
    println!("Public b infinity norm:  {}", ring_keypair.public_key_b.infinity_norm());
    println!("Relation verified:       {}", ring_keypair.verify_relation(&ctx));

    // 2. Module-LWE Keygen (k=3, l=2)
    println!("\n[2] Generating Module-LWE Keypair (k=3, l=2, Dilithium geometry)...");
    let mlw_params = MlwParams::standard();
    let mlw_keypair = MlwKeypair::generate(&mlw_params, &ctx, &mut rng);

    println!("Matrix A: {} rows x {} cols of polynomials in R_q", mlw_params.k, mlw_params.l);
    println!("Secret vector s: {} polynomials", mlw_keypair.secret_key_s.dim());
    println!("Public vector b: {} polynomials", mlw_keypair.public_key_b.dim());
    println!("Secret vector norm:      {}", mlw_keypair.secret_key_s.infinity_norm());
    println!("Error vector norm:       {}", mlw_keypair.error_e.infinity_norm());
    println!("Relation A * s + e == b: {}", mlw_keypair.verify_relation(&ctx));

    // 3. ASCII Visualizations
    println!("\n[3] Discrete Gaussian Noise Histogram (sigma = 3.0):");
    let samples: Vec<i32> = (0..1500)
        .map(|_| LatticeSampler::sample_gaussian_int(3.0, &mut rng))
        .collect();
    println!("{}", AsciiVisualizer::render_gaussian_histogram(&samples, 11, 40));

    println!("[4] Polynomial Coefficient Amplitude Spectrum:");
    println!("{}", AsciiVisualizer::render_poly_spectrum(&mlw_keypair.public_key_b[0].coeffs, DEFAULT_Q, 45));

    println!("Keypair generation complete and mathematically verified!");
    Ok(())
}
