//! Pure 7-bit ASCII Visualizer and Diagnostics for NirvanaLattice.
//!
//! Renders NTT butterfly networks, Gaussian noise histograms, polynomial spectra,
//! and cryptographic verification certificate cards.
//! STRICT INVARIANT: 100% pure 7-bit ASCII (0..=127). Strictly no Unicode box drawings.

use std::fmt::Write;

use crate::snark::proof::{LatticeProof, VerificationResult};
use crate::snark::prover::SnarkParams;

/// ASCII visualizer utility functions.
pub struct AsciiVisualizer;

impl AsciiVisualizer {
    /// Renders a Cooley-Tukey Radix-2 NTT Butterfly Network flow diagram.
    pub fn render_ntt_butterfly_diagram() -> String {
        let mut s = String::new();
        s.push_str("+-------------------------------------------------------------------------+\n");
        s.push_str("|            NEGACYCLIC NTT COOLEY-TUKEY BUTTERFLY NETWORK (RADIX-2)      |\n");
        s.push_str("+-------------------------------------------------------------------------+\n");
        s.push_str("| Stage 0 (Inputs)    | Stage 1 (N/4)      | Stage 2 (N/2)      | NTT Out |\n");
        s.push_str("+---------------------+--------------------+--------------------+---------+\n");
        s.push_str("| a[0] --->(+)--------+--->(+)-------------+--->(+)-------------+-> A[0]  |\n");
        s.push_str("|           ^         |     ^              |     ^              |         |\n");
        s.push_str("|           | (W^0)   |     | (W^0)        |     | (W^0)        |         |\n");
        s.push_str("| a[1] ---(-)+--------+     |              |     |              |-> A[1]  |\n");
        s.push_str("|                     |     |              |     |              |         |\n");
        s.push_str("| a[2] --->(+)--------+----+-)-------------+     |              |-> A[2]  |\n");
        s.push_str("|           ^         |     | (W^2)        |     |              |         |\n");
        s.push_str("|           | (W^1)   |     |              |     |              |         |\n");
        s.push_str("| a[3] ---(-)+--------+-----+--------------+----+-)-------------+-> A[3]  |\n");
        s.push_str("|                                          |     | (W^1)        |         |\n");
        s.push_str("| a[4] --->(+)--------+--->(+)-------------+     |              |-> A[4]  |\n");
        s.push_str("|           ^         |     ^              |     |              |         |\n");
        s.push_str("|           | (W^0)   |     | (W^0)        |     |              |         |\n");
        s.push_str("| a[5] ---(-)+--------+     |              |     |              |-> A[5]  |\n");
        s.push_str("|                     |     |              |     |              |         |\n");
        s.push_str("| a[6] --->(+)--------+----+-)-------------+     |              |-> A[6]  |\n");
        s.push_str("|           ^         |     | (W^2)        |     |              |         |\n");
        s.push_str("|           | (W^1)   |     |              |     |              |         |\n");
        s.push_str("| a[7] ---(-)+--------+-----+--------------+-----+--------------+-> A[7]  |\n");
        s.push_str("+-------------------------------------------------------------------------+\n");
        s.push_str("| Twiddle Factors: W^k = psi^(2k) mod q, where psi^256 = -1 (mod q)       |\n");
        s.push_str("| Butterfly equations: u' = u + W*v (mod q),  v' = u - W*v (mod q)        |\n");
        s.push_str("+-------------------------------------------------------------------------+\n");
        s
    }

    /// Renders a horizontal ASCII histogram of discrete Gaussian noise distribution.
    pub fn render_gaussian_histogram(samples: &[i32], bins_count: usize, max_bar_width: usize) -> String {
        if samples.is_empty() {
            return String::from("[No samples to display]\n");
        }

        let mut s = String::new();
        let min_val = *samples.iter().min().unwrap();
        let max_val = *samples.iter().max().unwrap();

        let span = (max_val - min_val).max(1) as f64;
        let step = span / (bins_count as f64);

        let mut counts = vec![0usize; bins_count];
        for &val in samples {
            let mut idx = ((val - min_val) as f64 / step).floor() as usize;
            if idx >= bins_count {
                idx = bins_count - 1;
            }
            counts[idx] += 1;
        }

        let max_count = *counts.iter().max().unwrap_or(&1);

        let mean: f64 = samples.iter().map(|&x| x as f64).sum::<f64>() / (samples.len() as f64);
        let variance: f64 = samples
            .iter()
            .map(|&x| {
                let diff = x as f64 - mean;
                diff * diff
            })
            .sum::<f64>()
            / (samples.len() as f64);
        let std_dev = variance.sqrt();

        s.push_str("+-------------------------------------------------------------------------+\n");
        s.push_str("|             LATTICE DISCRETE GAUSSIAN NOISE DISTRIBUTION                |\n");
        s.push_str("+-------------------------------------------------------------------------+\n");
        let _ = writeln!(
            s,
            "| Sample count: {:<8}  Mean: {:>+6.3}  StdDev: {:>6.3}  Min: {:>4}  Max: {:>4} |",
            samples.len(),
            mean,
            std_dev,
            min_val,
            max_val
        );
        s.push_str("+-------------------------------------------------------------------------+\n");
        s.push_str("|   Range       | Count | Histogram Spectrum                              |\n");
        s.push_str("+---------------+-------+-------------------------------------------------+\n");

        for (i, &count) in counts.iter().enumerate() {
            let r_start = min_val as f64 + (i as f64) * step;
            let r_end = r_start + step;
            let bar_len = if max_count > 0 {
                (count * max_bar_width) / max_count
            } else {
                0
            };
            let bar: String = "#".repeat(bar_len);

            let _ = writeln!(
                s,
                "| [{:>4.0}, {:>4.0}) | {:>5} | {:<47} |",
                r_start, r_end, count, bar
            );
        }

        s.push_str("+---------------+-------+-------------------------------------------------+\n");
        s
    }

    /// Renders an ASCII verification certificate card.
    pub fn render_verification_card(
        proof: &LatticeProof,
        result: &VerificationResult,
        params: &SnarkParams,
        elapsed_ms: f64,
    ) -> String {
        let mut s = String::new();
        let status_str = if result.is_valid() {
            "[PASSED] VALID ARGUMENT OF KNOWLEDGE"
        } else {
            "[FAILED] REJECTED BY VERIFIER"
        };

        let w_norm = proof.commitment_w.infinity_norm();
        let zs_norm = proof.response_z_s.infinity_norm();
        let ze_norm = proof.response_z_e.infinity_norm();
        let thresh_s = params.threshold_s();
        let thresh_e = params.threshold_e();

        s.push_str("+-------------------------------------------------------------------------+\n");
        s.push_str("|          NIRVANA-LATTICE POST-QUANTUM ZK-SNARK CERTIFICATE              |\n");
        s.push_str("+-------------------------------------------------------------------------+\n");
        let _ = writeln!(s, "| Protocol:  Lattice-Based zk-SNARK over Module-LWE (Ring R_q)            |");
        let _ = writeln!(s, "| Ring:      Z_q[X] / (X^256 + 1) with Dilithium prime q = 8380417        |");
        let _ = writeln!(s, "| Transform: Fast Negacyclic NTT in O(N log N) with Twiddle Factors       |");
        s.push_str("+-------------------------------------------------------------------------+\n");
        let _ = writeln!(s, "| Verification Status:   {:<48} |", status_str);
        let _ = writeln!(s, "| Verification Time:     {:<44.2} ms |", elapsed_ms);
        let _ = writeln!(s, "| Rejection Attempts:    {:<48} |", proof.attempts);
        s.push_str("+-------------------------------------------------------------------------+\n");
        let _ = writeln!(s, "| Commitment Norm ||w||_inf:      {:<39} |", w_norm);
        let _ = writeln!(
            s,
            "| Response Norm   ||z_s||_inf:    {:<12} (Threshold: {:<12}) |",
            zs_norm, thresh_s
        );
        let _ = writeln!(
            s,
            "| Noise Norm      ||z_e||_inf:    {:<12} (Threshold: {:<12}) |",
            ze_norm, thresh_e
        );
        let _ = writeln!(s, "| Challenge Weight tau:           {:<39} |", params.tau);
        let _ = writeln!(s, "| Proof Binary Size:              {:<41} bytes |", proof.size_bytes());
        s.push_str("+-------------------------------------------------------------------------+\n");
        let _ = writeln!(s, "| Quantum Hardness: Module-SIS / Module-LWE (NIST Security Level 3/5)     |");
        s.push_str("+-------------------------------------------------------------------------+\n");
        s
    }

    /// Renders an ASCII spectrum bar of polynomial coefficient amplitudes.
    pub fn render_poly_spectrum(coeffs: &[i32], q: i32, width: usize) -> String {
        let mut s = String::new();
        s.push_str("+-------------------------------------------------------------------------+\n");
        s.push_str("|                     POLYNOMIAL COEFFICIENT SPECTRUM                     |\n");
        s.push_str("+-------------------------------------------------------------------------+\n");

        let chunk_size = (coeffs.len() / width).max(1);
        let half_q = q / 2;

        let mut row_chars = String::new();
        for chunk in coeffs.chunks(chunk_size) {
            let max_c = chunk
                .iter()
                .map(|&c| {
                    let centered = if c > half_q { c - q } else { c };
                    centered.abs()
                })
                .max()
                .unwrap_or(0);

            let ch = if max_c == 0 {
                '.'
            } else if max_c < 10 {
                '-'
            } else if max_c < 100 {
                '+'
            } else if max_c < 1000 {
                '*'
            } else if max_c < 10000 {
                '='
            } else {
                '#'
            };
            row_chars.push(ch);
        }

        let _ = writeln!(s, "| Spectrum: [{:<60}] |", row_chars);
        s.push_str("| Legend: . (0)  - (<10)  + (<100)  * (<1k)  = (<10k)  # (>=10k)          |\n");
        s.push_str("+-------------------------------------------------------------------------+\n");
        s
    }
}
