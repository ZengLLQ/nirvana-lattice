//! Fiat-Shamir challenge generation using Blake3 cryptographic hashing.
//!
//! Transforms the interactive 3-pass sigma-protocol into a non-interactive
//! succinct argument of knowledge (zk-SNARK) in the Random Oracle Model.

use blake3::Hasher;

use crate::lwe::matrix::{PolyMatrix, PolyVector};
use crate::ring::poly::Poly;
use crate::ring::reduction::freeze;

/// Fiat-Shamir random oracle challenge generator.
pub struct FiatShamirOracle;

impl FiatShamirOracle {
    /// Computes a deterministic sparse ternary challenge polynomial c in R_q
    /// from public matrix A, public key b, and prover commitment w.
    ///
    /// Challenge polynomial has degree N and exactly `tau` coefficients in {-1, 1},
    /// with all other coefficients equal to 0.
    pub fn compute_challenge(
        matrix_a: &PolyMatrix,
        vector_b: &PolyVector,
        commitment_w: &PolyVector,
        tau: usize,
    ) -> Poly {
        let n = commitment_w.poly_len();
        let q = commitment_w.q;
        assert!(tau <= n, "Challenge weight tau cannot exceed polynomial degree");

        let mut hasher = Hasher::new();
        hasher.update(b"NIRVANA_LATTICE_FIAT_SHAMIR_V1");

        // Hash public matrix A
        for row in &matrix_a.data {
            for poly in row {
                for &c in &poly.coeffs {
                    hasher.update(&c.to_le_bytes());
                }
            }
        }

        // Hash public vector b
        for poly in &vector_b.polys {
            for &c in &poly.coeffs {
                hasher.update(&c.to_le_bytes());
            }
        }

        // Hash prover commitment w
        for poly in &commitment_w.polys {
            for &c in &poly.coeffs {
                hasher.update(&c.to_le_bytes());
            }
        }

        let mut output_reader = hasher.finalize_xof();
        let mut coeffs = vec![0i32; n];
        let mut positions: Vec<usize> = Vec::with_capacity(tau);

        // Derive deterministic pseudo-random unique positions using XOF stream
        let mut byte_buf = [0u8; 4];
        while positions.len() < tau {
            output_reader.fill(&mut byte_buf);
            let raw_idx = u32::from_le_bytes(byte_buf) as usize;
            let idx = raw_idx % n;
            if !positions.contains(&idx) {
                positions.push(idx);
            }
        }

        // Assign +/- 1 signs to the selected positions
        let mut sign_byte = [0u8; 1];
        for pos in positions {
            output_reader.fill(&mut sign_byte);
            let sign = if (sign_byte[0] & 1) == 1 { 1 } else { -1 };
            coeffs[pos] = freeze(sign as i64, q);
        }

        Poly { coeffs, q }
    }

    /// Hashes arbitrary serialized data to a 32-byte commitment digest.
    pub fn hash_digest(data: &[u8]) -> [u8; 32] {
        *blake3::hash(data).as_bytes()
    }
}
