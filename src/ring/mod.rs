//! Cyclotomic quotient ring R_q = Z_q[X] / (X^N + 1) arithmetic.
//!
//! Includes modular reduction, negacyclic NTT, and polynomial operations.

pub mod ntt;
pub mod poly;
pub mod reduction;

pub use ntt::NttContext;
pub use poly::{NttPoly, Poly};
pub use reduction::{
    barrett_reduce, center, freeze, mod_add, mod_inv, mod_mul, mod_neg, mod_pow, mod_sub,
    montgomery_reduce, BarrettReducer, MontgomeryReducer, DEFAULT_N, DEFAULT_Q, RLWE_Q,
};
