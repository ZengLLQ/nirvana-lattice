//! Lattice-based zero-knowledge succinct non-interactive arguments of knowledge (zk-SNARK).
//!
//! Provides Fiat-Shamir transformation, non-interactive prover, verifier, and proof structures.

pub mod fiat_shamir;
pub mod proof;
pub mod prover;
pub mod verifier;

pub use fiat_shamir::FiatShamirOracle;
pub use proof::{LatticeProof, VerificationResult};
pub use prover::{LatticeProver, SnarkParams};
pub use verifier::LatticeVerifier;
