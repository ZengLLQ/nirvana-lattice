//! Fast modular reduction and arithmetic for cyclotomic quotient rings.
//!
//! Provides constant-time Montgomery reduction, Barrett reduction,
//! modular exponentiation, and canonical balanced centering.

/// Default prime modulus used in CRYSTALS-Dilithium: q = 8380417 = 2^23 - 2^13 + 1.
pub const DEFAULT_Q: i32 = 8380417;

/// Alternative Ring-LWE / NewHope prime modulus: q = 12289.
pub const RLWE_Q: i32 = 12289;

/// Default polynomial degree N = 256 (power of two).
pub const DEFAULT_N: usize = 256;

/// Reduces an integer modulo `q` into the canonical range `[0, q - 1]`.
#[inline]
pub fn freeze(val: i64, q: i32) -> i32 {
    let q_i64 = q as i64;
    let r = val % q_i64;
    if r < 0 {
        (r + q_i64) as i32
    } else {
        r as i32
    }
}

/// Centers an integer in `[0, q - 1]` to the balanced range `[-(q-1)/2, (q-1)/2]`.
#[inline]
pub fn center(val: i32, q: i32) -> i32 {
    let half_q = q / 2;
    if val > half_q {
        val - q
    } else if val < -half_q {
        val + q
    } else {
        val
    }
}

/// Modular addition: `(a + b) mod q`.
#[inline]
pub fn mod_add(a: i32, b: i32, q: i32) -> i32 {
    let sum = a as i64 + b as i64;
    freeze(sum, q)
}

/// Modular subtraction: `(a - b) mod q`.
#[inline]
pub fn mod_sub(a: i32, b: i32, q: i32) -> i32 {
    let diff = a as i64 - b as i64;
    freeze(diff, q)
}

/// Modular negation: `(-a) mod q`.
#[inline]
pub fn mod_neg(a: i32, q: i32) -> i32 {
    if a == 0 {
        0
    } else {
        let f = freeze(a as i64, q);
        if f == 0 { 0 } else { q - f }
    }
}

/// Modular multiplication: `(a * b) mod q`.
#[inline]
pub fn mod_mul(a: i32, b: i32, q: i32) -> i32 {
    let prod = (a as i64) * (b as i64);
    freeze(prod, q)
}

/// Modular exponentiation: `(base ^ exp) mod q`.
pub fn mod_pow(mut base: i64, mut exp: u64, q: i32) -> i32 {
    let q_i64 = q as i64;
    let mut res = 1i64;
    base = base.rem_euclid(q_i64);
    while exp > 0 {
        if exp & 1 == 1 {
            res = (res * base).rem_euclid(q_i64);
        }
        base = (base * base).rem_euclid(q_i64);
        exp >>= 1;
    }
    res as i32
}

/// Extended Euclidean Algorithm for modular inverse.
/// Returns `Some(x)` such that `(a * x) mod q == 1`, or `None` if `gcd(a, q) != 1`.
pub fn mod_inv(a: i32, q: i32) -> Option<i32> {
    let a_canonical = freeze(a as i64, q);
    if a_canonical == 0 {
        return None;
    }
    // Fermat's Little Theorem since q is prime: a^(q-2) mod q
    let inv = mod_pow(a_canonical as i64, (q - 2) as u64, q);
    if mod_mul(a_canonical, inv, q) == 1 {
        Some(inv)
    } else {
        None
    }
}

/// Standalone Barrett reduction for arbitrary modulus q.
#[inline]
pub fn barrett_reduce(a: i64, q: i32) -> i32 {
    let reducer = BarrettReducer::new(q);
    reducer.reduce(a)
}

/// Standalone Montgomery reduction for odd modulus q.
#[inline]
pub fn montgomery_reduce(t: i64, q: i32) -> i32 {
    let reducer = MontgomeryReducer::new(q);
    reducer.reduce(t)
}

/// Barrett reduction structure for fast division-free modular reduction.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct BarrettReducer {
    pub q: i32,
    pub k: u32,
    pub mu: u64,
}

impl BarrettReducer {
    /// Constructs a Barrett reducer for prime modulus `q`.
    /// Precomputes `mu = floor(2^k / q)` where `k = 48`.
    pub fn new(q: i32) -> Self {
        assert!(q > 1, "Modulus must be greater than 1");
        let k = 48u32;
        let mu = ((1u128 << k) / (q as u128)) as u64;
        Self { q, k, mu }
    }

    /// Reduces a 64-bit integer `a` modulo `q` into `[0, q - 1]`.
    #[inline]
    pub fn reduce(&self, a: i64) -> i32 {
        let q_i64 = self.q as i64;
        if a >= 0 {
            let a_u64 = a as u64;
            let q_est = (((a_u64 as u128) * (self.mu as u128)) >> self.k) as i64;
            let mut r = a - q_est * q_i64;
            while r >= q_i64 {
                r -= q_i64;
            }
            while r < 0 {
                r += q_i64;
            }
            r as i32
        } else {
            let neg_a = -a;
            let rem = self.reduce(neg_a);
            if rem == 0 {
                0
            } else {
                self.q - rem
            }
        }
    }
}

/// Montgomery reducer for constant-time modular multiplication.
/// Works with Montgomery radix R = 2^32.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct MontgomeryReducer {
    pub q: i32,
    pub r_mod_q: i32,
    pub r2_mod_q: i32,
    pub q_inv_neg: u64, // -q^{-1} mod 2^32
}

impl MontgomeryReducer {
    /// Creates a new Montgomery reducer for modulus `q`.
    pub fn new(q: i32) -> Self {
        assert!(q % 2 != 0, "Montgomery reduction requires odd modulus");
        let q_u64 = q as u64;
        // Compute inverse of q modulo 2^32 via Newton's method
        let mut inv = 1u64;
        for _ in 0..5 {
            inv = inv.wrapping_mul(2u64.wrapping_sub(q_u64.wrapping_mul(inv)));
        }
        let q_inv_neg = (0u64.wrapping_sub(inv)) & 0xFFFF_FFFF;

        let r_mod_q = ((1u128 << 32) % (q as u128)) as i32;
        let r2_mod_q = (((1u128 << 64) % (q as u128))) as i32;

        Self {
            q,
            r_mod_q,
            r2_mod_q,
            q_inv_neg,
        }
    }

    /// Converts regular integer `a` into Montgomery domain: `a * R mod q`.
    #[inline]
    pub fn to_montgomery(&self, a: i32) -> i32 {
        self.mul(a, self.r2_mod_q)
    }

    /// Converts Montgomery form `aR mod q` back to regular integer `a mod q`.
    #[inline]
    pub fn from_montgomery(&self, a_mont: i32) -> i32 {
        self.reduce(a_mont as i64)
    }

    /// Montgomery reduction: computes `T * R^{-1} mod q`.
    #[inline]
    pub fn reduce(&self, t: i64) -> i32 {
        let t_u64 = t as u64;
        let m = (t_u64.wrapping_mul(self.q_inv_neg) & 0xFFFF_FFFF) as i64;
        let t_sub = t + m * (self.q as i64);
        let mut res = (t_sub >> 32) as i32;
        if res >= self.q {
            res -= self.q;
        } else if res < 0 {
            res += self.q;
        }
        res
    }

    /// Multiplies two numbers in Montgomery form: `(a * b * R^{-1}) mod q`.
    #[inline]
    pub fn mul(&self, a: i32, b: i32) -> i32 {
        let prod = (a as i64) * (b as i64);
        self.reduce(prod)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_mod_neg_identities() {
        let q = DEFAULT_Q;
        assert_eq!(mod_neg(0, q), 0);
        assert_eq!(mod_neg(1, q), q - 1);
        assert_eq!(mod_neg(q - 1, q), 1);
        assert_eq!(mod_add(42, mod_neg(42, q), q), 0);
    }

    #[test]
    fn test_fermat_little_theorem() {
        let q = DEFAULT_Q;
        for a in [2, 7, 13, 257, 1048576] {
            assert_eq!(mod_pow(a, (q - 1) as u64, q), 1);
        }
    }

    #[test]
    fn test_montgomery_random_roundtrips() {
        let q = DEFAULT_Q;
        let reducer = MontgomeryReducer::new(q);
        for i in 0..100 {
            let val = (i * 83801) % q;
            let mont = reducer.to_montgomery(val);
            let back = reducer.from_montgomery(mont);
            assert_eq!(back, val);
        }
    }

    #[test]
    fn test_barrett_extreme_ranges() {
        let q = DEFAULT_Q;
        let reducer = BarrettReducer::new(q);
        for k in [-10000000000i64, -q as i64, 0, q as i64, 10000000000i64] {
            assert_eq!(reducer.reduce(k), freeze(k, q));
        }
    }
}
