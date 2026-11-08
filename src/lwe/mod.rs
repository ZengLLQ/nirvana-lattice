//! Module-LWE and Module-SIS cryptographic primitives.
//!
//! Provides polynomial vectors and matrices, discrete Gaussian samplers,
//! key generation, and M-SIS hardness checking.

pub mod keygen;
pub mod matrix;
pub mod sampler;

pub use keygen::{verify_msis, MlwKeypair, MlwParams};
pub use matrix::{PolyMatrix, PolyVector};
pub use sampler::LatticeSampler;
