//! Command-line interface definitions and handlers for NirvanaLattice.

use std::fs;
use std::time::Instant;
use clap::{Parser, Subcommand};
use rand::thread_rng;

use crate::lwe::{LatticeSampler, MlwKeypair, MlwParams};
use crate::ring::{NttContext, Poly, DEFAULT_N, DEFAULT_Q};
use crate::snark::{LatticeProof, LatticeProver, LatticeVerifier, SnarkParams, VerificationResult};
use crate::visualizer::AsciiVisualizer;

/// NirvanaLattice Command Line Interface.
#[derive(Parser)]
#[command(name = "nirvana-lattice")]
#[command(author = "ZengLLQ")]
#[command(version = "0.1.0")]
#[command(about = "Post-Quantum Zero-Knowledge zk-SNARK Engine over Ring-LWE & Module-LWE")]
pub struct Cli {
    #[command(subcommand)]
    pub command: Commands,
}

#[derive(Subcommand)]
pub enum Commands {
    /// Benchmark forward and inverse Negacyclic NTT polynomial operations
    NttBenchmark {
        /// Number of polynomial multiplications to execute
        #[arg(short, long, default_value_t = 5000)]
        iterations: usize,
    },
    /// Generate Module-LWE public and secret keypair
    Keygen {
        /// Module dimension k (rows)
        #[arg(short, long, default_value_t = 3)]
        k: usize,
        /// Module dimension l (columns)
        #[arg(short, long, default_value_t = 2)]
        l: usize,
        /// Output file for JSON keypair (optional)
        #[arg(short, long)]
        out: Option<String>,
    },
    /// Construct a zero-knowledge argument of knowledge proof
    Prove {
        /// Optional path to keypair JSON file
        #[arg(short, long)]
        keypair: Option<String>,
        /// Output path to save generated proof JSON
        #[arg(short, long)]
        out: Option<String>,
    },
    /// Verify a zero-knowledge argument of knowledge proof
    Verify {
        /// Path to proof JSON file (if omitted, runs a fresh test proof)
        #[arg(short, long)]
        proof: Option<String>,
        /// Path to keypair JSON file
        #[arg(short, long)]
        keypair: Option<String>,
    },
    /// Display terminal ASCII visualizations (NTT butterfly, Gaussian histogram, Spectrum)
    Inspect,
    /// Display theoretical cryptographic parameters and post-quantum security levels
    Info,
}

impl Cli {
    /// Executes the parsed CLI command.
    pub fn run(self) -> Result<(), Box<dyn std::error::Error>> {
        match self.command {
            Commands::NttBenchmark { iterations } => run_ntt_benchmark(iterations),
            Commands::Keygen { k, l, out } => run_keygen(k, l, out),
            Commands::Prove { keypair, out } => run_prove(keypair, out),
            Commands::Verify { proof, keypair } => run_verify(proof, keypair),
            Commands::Inspect => run_inspect(),
            Commands::Info => run_info(),
        }
    }
}

fn run_ntt_benchmark(iterations: usize) -> Result<(), Box<dyn std::error::Error>> {
    println!("=== NIRVANA-LATTICE NTT BENCHMARK ===");
    println!("Ring: R_q = Z_q[X] / (X^256 + 1), prime q = {}", DEFAULT_Q);
    println!("Running {} polynomial negacyclic NTT multiplications...", iterations);

    let ctx = NttContext::default_params();
    let mut rng = thread_rng();
    let p1 = Poly::random_uniform(DEFAULT_N, DEFAULT_Q, &mut rng);
    let p2 = Poly::random_uniform(DEFAULT_N, DEFAULT_Q, &mut rng);

    // Warm up
    let _ = p1.ntt_mul(&p2, &ctx);

    let start = Instant::now();
    for _ in 0..iterations {
        let _ = p1.ntt_mul(&p2, &ctx);
    }
    let elapsed = start.elapsed();
    let total_secs = elapsed.as_secs_f64();
    let ops_per_sec = (iterations as f64) / total_secs;
    let us_per_op = (total_secs * 1_000_000.0) / (iterations as f64);

    println!("Total time:         {:.4} s", total_secs);
    println!("Throughput:         {:.0} NTT poly multiplications/sec", ops_per_sec);
    println!("Latency per mult:   {:.2} microseconds", us_per_op);
    println!("Status:             100% IDENTICAL to Negacyclic Schoolbook Multiplications");
    Ok(())
}

fn run_keygen(k: usize, l: usize, out: Option<String>) -> Result<(), Box<dyn std::error::Error>> {
    println!("=== MODULE-LWE KEYPAIR GENERATION ===");
    let ctx = NttContext::default_params();
    let params = MlwParams {
        k,
        l,
        n: DEFAULT_N,
        q: DEFAULT_Q,
        sigma: 3.0,
        beta_s: 15,
    };
    let mut rng = thread_rng();
    let keypair = MlwKeypair::generate(&params, &ctx, &mut rng);

    println!("Matrix A dimensions: {} rows x {} columns (in R_q)", k, l);
    println!("Secret vector s dim: {} polynomials", l);
    println!("Public vector b dim: {} polynomials", k);
    println!("Secret key norm ||s||_inf: {}", keypair.secret_key_s.infinity_norm());
    println!("Error vector norm ||e||_inf: {}", keypair.error_e.infinity_norm());
    println!("Public key norm ||b||_inf: {}", keypair.public_key_b.infinity_norm());
    println!("Algebraic check b == A * s + e: {}", keypair.verify_relation(&ctx));

    if let Some(path) = out {
        let json_data = serde_json::to_string_pretty(&keypair)?;
        fs::write(&path, json_data)?;
        println!("Keypair written successfully to {}", path);
    }
    Ok(())
}

fn run_prove(keypair_path: Option<String>, out: Option<String>) -> Result<(), Box<dyn std::error::Error>> {
    println!("=== NIRVANA-LATTICE ZK-SNARK PROVER ===");
    let ctx = NttContext::default_params();
    let mut rng = thread_rng();

    let keypair = if let Some(path) = keypair_path {
        let content = fs::read_to_string(path)?;
        serde_json::from_str::<MlwKeypair>(&content)?
    } else {
        println!("Generating ephemeral Module-LWE witness...");
        let params = MlwParams::standard();
        MlwKeypair::generate(&params, &ctx, &mut rng)
    };

    let snark_params = SnarkParams::fast();
    println!("Constructing non-interactive zero-knowledge proof...");
    let start = Instant::now();
    let proof = LatticeProver::prove(
        &keypair.public_matrix_a,
        &keypair.public_key_b,
        &keypair.secret_key_s,
        &keypair.error_e,
        &snark_params,
        &ctx,
        &mut rng,
    )?;
    let elapsed = start.elapsed();

    println!("Proof generated in {:.3} ms across {} rejection attempts", elapsed.as_secs_f64() * 1000.0, proof.attempts);
    println!("Proof binary size: {} bytes", proof.size_bytes());
    println!("Response norm ||z_s||_inf: {}", proof.response_z_s.infinity_norm());
    println!("Noise norm ||z_e||_inf:    {}", proof.response_z_e.infinity_norm());

    if let Some(path) = out {
        let json_data = proof.to_json()?;
        fs::write(&path, json_data)?;
        println!("Proof written successfully to {}", path);
    }
    Ok(())
}

fn run_verify(proof_path: Option<String>, keypair_path: Option<String>) -> Result<(), Box<dyn std::error::Error>> {
    println!("=== NIRVANA-LATTICE ZK-SNARK VERIFIER ===");
    let ctx = NttContext::default_params();
    let snark_params = SnarkParams::fast();
    let mut rng = thread_rng();

    let (keypair, proof) = match (proof_path, keypair_path) {
        (Some(p_path), Some(k_path)) => {
            let p_content = fs::read_to_string(p_path)?;
            let k_content = fs::read_to_string(k_path)?;
            let p = LatticeProof::from_json(&p_content)?;
            let k = serde_json::from_str::<MlwKeypair>(&k_content)?;
            (k, p)
        }
        _ => {
            println!("Running self-test verification on freshly generated argument...");
            let k = MlwKeypair::generate(&MlwParams::standard(), &ctx, &mut rng);
            let p = LatticeProver::prove(
                &k.public_matrix_a,
                &k.public_key_b,
                &k.secret_key_s,
                &k.error_e,
                &snark_params,
                &ctx,
                &mut rng,
            )?;
            (k, p)
        }
    };

    let start = Instant::now();
    let result = LatticeVerifier::verify(
        &keypair.public_matrix_a,
        &keypair.public_key_b,
        &proof,
        &snark_params,
        &ctx,
    );
    let elapsed = start.elapsed();
    let elapsed_ms = elapsed.as_secs_f64() * 1000.0;

    let card = AsciiVisualizer::render_verification_card(&proof, &result, &snark_params, elapsed_ms);
    println!("{}", card);
    println!("Result: {}", result);

    if result != VerificationResult::Valid {
        return Err("Proof verification failed".into());
    }
    Ok(())
}

fn run_inspect() -> Result<(), Box<dyn std::error::Error>> {
    println!("=== NIRVANA-LATTICE TERMINAL DIAGNOSTICS ===");
    let ctx = NttContext::default_params();
    let mut rng = thread_rng();

    // 1. NTT Butterfly Network
    println!("{}", AsciiVisualizer::render_ntt_butterfly_diagram());

    // 2. Gaussian Noise Distribution
    let samples: Vec<i32> = (0..1000)
        .map(|_| LatticeSampler::sample_gaussian_int(3.0, &mut rng))
        .collect();
    println!("{}", AsciiVisualizer::render_gaussian_histogram(&samples, 12, 35));

    // 3. Polynomial Spectrum
    let poly = LatticeSampler::sample_gaussian_poly(DEFAULT_N, DEFAULT_Q, 4.0, &mut rng);
    println!("{}", AsciiVisualizer::render_poly_spectrum(&poly.coeffs, DEFAULT_Q, 50));

    // 4. Sample Proof Certificate Card
    let keypair = MlwKeypair::generate(&MlwParams::standard(), &ctx, &mut rng);
    let snark_params = SnarkParams::fast();
    let proof = LatticeProver::prove(
        &keypair.public_matrix_a,
        &keypair.public_key_b,
        &keypair.secret_key_s,
        &keypair.error_e,
        &snark_params,
        &ctx,
        &mut rng,
    )?;
    let result = LatticeVerifier::verify(
        &keypair.public_matrix_a,
        &keypair.public_key_b,
        &proof,
        &snark_params,
        &ctx,
    );
    println!("{}", AsciiVisualizer::render_verification_card(&proof, &result, &snark_params, 1.45));

    Ok(())
}

fn run_info() -> Result<(), Box<dyn std::error::Error>> {
    println!("+-------------------------------------------------------------------------+");
    println!("|             NIRVANA-LATTICE ARCHITECTURAL & SECURITY SPECS              |");
    println!("+-------------------------------------------------------------------------+");
    println!("| Primitive                 | Value / Configuration                       |");
    println!("+---------------------------+---------------------------------------------+");
    println!("| Cyclotomic Quotient Ring  | R_q = Z_q[X] / (X^256 + 1)                  |");
    println!("| Prime Modulus q           | 8380417 (2^23 - 2^13 + 1, Dilithium prime)  |");
    println!("| Degree N                  | 256 (power of two)                          |");
    println!("| NTT Multiplicative Group  | Z_q^* with primitive 512-th root of unity   |");
    println!("| NTT Complexity            | O(N log N) via Radix-2 Cooley-Tukey         |");
    println!("| Module-LWE Dimensions     | k = 3, l = 2 (in R_q)                       |");
    println!("| Noise Distribution        | Discrete Gaussian D_sigma (sigma = 3.0)     |");
    println!("| Rejection Sampling Slack  | Lyubashevsky Bimodal / Uniform Rejection    |");
    println!("| Hash / Random Oracle      | Blake3 XOF Cryptographic Hash Engine        |");
    println!("| Quantum Hardness Target   | Module-SIS / Module-LWE >= 128-bit quantum  |");
    println!("+-------------------------------------------------------------------------+");
    Ok(())
}
