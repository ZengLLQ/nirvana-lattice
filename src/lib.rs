//! NirvanaLattice: Lattice-based zero-knowledge arguments of knowledge (zk-SNARK)
//! over Ring-LWE & Module-LWE with negacyclic NTT polynomial arithmetic in Rust.

pub mod cli;
pub mod commitment;
pub mod lwe;
pub mod ring;
pub mod snark;
pub mod visualizer;

pub fn version() -> &'static str {
    env!("CARGO_PKG_VERSION")
}
