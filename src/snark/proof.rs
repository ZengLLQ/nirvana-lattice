//! Compact proof data structures and serialization for lattice zk-SNARKs.

use std::fmt;
use serde::{Deserialize, Serialize};

use crate::lwe::matrix::PolyVector;
use crate::ring::poly::Poly;

/// Verification outcome of a lattice zk-SNARK proof.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub enum VerificationResult {
    Valid,
    InvalidNorm {
        variable: String,
        observed: i32,
        threshold: i32,
    },
    ChallengeMismatch,
    DimensionMismatch {
        expected: usize,
        observed: usize,
    },
}

impl VerificationResult {
    pub fn is_valid(&self) -> bool {
        matches!(self, VerificationResult::Valid)
    }
}

impl fmt::Display for VerificationResult {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            VerificationResult::Valid => write!(f, "VALID (Proof accepted)"),
            VerificationResult::InvalidNorm {
                variable,
                observed,
                threshold,
            } => write!(
                f,
                "REJECTED: Norm check failed for {} (observed ||{}||_inf = {}, bound = {})",
                variable, variable, observed, threshold
            ),
            VerificationResult::ChallengeMismatch => {
                write!(f, "REJECTED: Fiat-Shamir challenge polynomial mismatch")
            }
            VerificationResult::DimensionMismatch {
                expected,
                observed,
            } => write!(
                f,
                "REJECTED: Dimension mismatch (expected {}, got {})",
                expected, observed
            ),
        }
    }
}

/// Zero-Knowledge Succinct Argument of Knowledge (zk-SNARK) proof structure.
///
/// Proves knowledge of secret witness `s` and small error `e` such that `b = A * s + e (mod q)`.
#[derive(Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct LatticeProof {
    pub commitment_w: PolyVector, // Prover commitment w = A * y_s + y_e (mod q)
    pub challenge_c: Poly,        // Fiat-Shamir challenge polynomial
    pub response_z_s: PolyVector, // Masked secret response z_s = y_s + c * s
    pub response_z_e: PolyVector, // Masked error response z_e = y_e + c * e
    pub attempts: usize,          // Rejection sampling attempts required to satisfy norm bound
}

impl LatticeProof {
    /// Serializes the proof to a compact JSON string.
    pub fn to_json(&self) -> Result<String, serde_json::Error> {
        serde_json::to_string(self)
    }

    /// Deserializes proof from a JSON string.
    pub fn from_json(json_str: &str) -> Result<Self, serde_json::Error> {
        serde_json::from_str(json_str)
    }

    /// Serializes proof to bytes.
    pub fn to_bytes(&self) -> Result<Vec<u8>, serde_json::Error> {
        serde_json::to_vec(self)
    }

    /// Deserializes proof from bytes.
    pub fn from_bytes(bytes: &[u8]) -> Result<Self, serde_json::Error> {
        serde_json::from_slice(bytes)
    }

    /// Approximate proof size in bytes.
    pub fn size_bytes(&self) -> usize {
        self.to_bytes().map(|b| b.len()).unwrap_or(0)
    }

    /// Maximum infinity norm across all response polynomials.
    pub fn max_response_norm(&self) -> i32 {
        self.response_z_s
            .infinity_norm()
            .max(self.response_z_e.infinity_norm())
    }
}

impl fmt::Debug for LatticeProof {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "LatticeProof(w_dim={}, z_s_dim={}, z_e_dim={}, max_norm={}, attempts={})",
            self.commitment_w.dim(),
            self.response_z_s.dim(),
            self.response_z_e.dim(),
            self.max_response_norm(),
            self.attempts
        )
    }
}
