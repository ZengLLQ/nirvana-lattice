//! Negacyclic Number Theoretic Transform (NTT) over cyclotomic rings.
//!
//! Implements fast negacyclic NTT and Inverse NTT in O(N log N) time
//! over R_q = Z_q[X] / (X^N + 1) using precomputed twiddle factors
//! and radix-2 Cooley-Tukey / Gentleman-Sande butterflies.

use crate::ring::reduction::{
    mod_add, mod_inv, mod_mul, mod_pow, mod_sub, DEFAULT_N, DEFAULT_Q,
};

/// Precomputed context for Negacyclic NTT and INTT.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct NttContext {
    pub n: usize,
    pub q: i32,
    pub psi: i32,         // Primitive 2N-th root of unity modulo q
    pub psi_inv: i32,     // Inverse of psi modulo q
    pub omega: i32,       // Primitive N-th root of unity (psi^2 mod q)
    pub omega_inv: i32,   // Inverse of omega modulo q
    pub n_inv: i32,       // Inverse of N modulo q
    pub psi_powers: Vec<i32>,      // psi^i mod q for i in 0..N
    pub psi_inv_powers: Vec<i32>,  // psi^(-i) mod q for i in 0..N
    pub twiddles: Vec<i32>,        // twiddle factors for cyclic NTT
    pub inv_twiddles: Vec<i32>,    // twiddle factors for cyclic INTT
    pub bit_rev: Vec<usize>,       // Bit-reversal permutation indices
}

impl NttContext {
    /// Constructs an NTT context for polynomial degree `n` and prime modulus `q`.
    /// `n` must be a power of two, and `q` must satisfy `q = 1 (mod 2n)`.
    pub fn new(n: usize, q: i32) -> Self {
        assert!(n.is_power_of_two(), "NTT degree n must be a power of two");
        let two_n = (2 * n) as i64;
        assert_eq!(
            (q as i64 - 1) % two_n,
            0,
            "Modulus q must satisfy q = 1 (mod 2n)"
        );

        let psi = find_primitive_2n_root(n, q).expect("Failed to find primitive 2N-th root of unity");
        let psi_inv = mod_inv(psi, q).expect("psi must have modular inverse");
        let omega = mod_mul(psi, psi, q);
        let omega_inv = mod_inv(omega, q).expect("omega must have modular inverse");
        let n_inv = mod_inv(n as i32, q).expect("n must have modular inverse modulo q");

        // Precompute powers of psi and psi^{-1} for negacyclic twisting/untwisting
        let mut psi_powers = Vec::with_capacity(n);
        let mut psi_inv_powers = Vec::with_capacity(n);
        let mut cur_psi = 1;
        let mut cur_psi_inv = 1;
        for _ in 0..n {
            psi_powers.push(cur_psi);
            psi_inv_powers.push(cur_psi_inv);
            cur_psi = mod_mul(cur_psi, psi, q);
            cur_psi_inv = mod_mul(cur_psi_inv, psi_inv, q);
        }

        // Bit-reversal permutation table
        let log_n = n.trailing_zeros() as usize;
        let mut bit_rev = Vec::with_capacity(n);
        for i in 0..n {
            bit_rev.push(reverse_bits(i, log_n));
        }

        // Precompute twiddle factors for standard Cooley-Tukey
        let mut twiddles = vec![0i32; n];
        let mut inv_twiddles = vec![0i32; n];
        for i in 0..n {
            twiddles[i] = mod_pow(omega as i64, i as u64, q);
            inv_twiddles[i] = mod_pow(omega_inv as i64, i as u64, q);
        }

        Self {
            n,
            q,
            psi,
            psi_inv,
            omega,
            omega_inv,
            n_inv,
            psi_powers,
            psi_inv_powers,
            twiddles,
            inv_twiddles,
            bit_rev,
        }
    }

    /// Creates default NTT context with N = 256 and Dilithium prime q = 8380417.
    pub fn default_params() -> Self {
        Self::new(DEFAULT_N, DEFAULT_Q)
    }

    /// In-place forward Negacyclic NTT on coefficients of polynomial `a`.
    /// Transforms `a` from coefficient domain to evaluation domain (NTT domain).
    pub fn forward(&self, a: &mut [i32]) {
        assert_eq!(a.len(), self.n, "Array length must match NTT degree N");

        // Step 1: Negacyclic twist: a[j] = a[j] * psi^j mod q
        for j in 0..self.n {
            a[j] = mod_mul(a[j], self.psi_powers[j], self.q);
        }

        // Step 2: Cyclic Cooley-Tukey NTT
        self.cyclic_ntt(a, false);
    }

    /// In-place inverse Negacyclic NTT (INTT) on polynomial `a`.
    /// Transforms `a` from NTT domain back to standard coefficient domain.
    pub fn inverse(&self, a: &mut [i32]) {
        assert_eq!(a.len(), self.n, "Array length must match NTT degree N");

        // Step 1: Cyclic INTT
        self.cyclic_ntt(a, true);

        // Step 2: Scale by 1/N and untwist: a[j] = a[j] * N^{-1} * psi^{-j} mod q
        for j in 0..self.n {
            let scaled = mod_mul(a[j], self.n_inv, self.q);
            a[j] = mod_mul(scaled, self.psi_inv_powers[j], self.q);
        }
    }

    /// Radix-2 Decimation-In-Frequency (DIF) / Decimation-In-Time (DIT) Cyclic NTT.
    fn cyclic_ntt(&self, a: &mut [i32], inverse: bool) {
        let n = self.n;

        // Apply bit-reversal permutation
        for i in 0..n {
            let rev = self.bit_rev[i];
            if i < rev {
                a.swap(i, rev);
            }
        }

        // Cooley-Tukey butterflies
        let mut len = 2;
        while len <= n {
            let half = len / 2;
            let step = n / len;
            for i in (0..n).step_by(len) {
                for j in 0..half {
                    let twiddle_idx = j * step;
                    let w = if inverse {
                        self.inv_twiddles[twiddle_idx]
                    } else {
                        self.twiddles[twiddle_idx]
                    };
                    let u = a[i + j];
                    let v = mod_mul(a[i + j + half], w, self.q);
                    a[i + j] = mod_add(u, v, self.q);
                    a[i + j + half] = mod_sub(u, v, self.q);
                }
            }
            len <<= 1;
        }
    }

    /// Pointwise multiplication of two polynomials in the NTT domain.
    pub fn pointwise_mul(&self, a: &[i32], b: &[i32], out: &mut [i32]) {
        assert_eq!(a.len(), self.n);
        assert_eq!(b.len(), self.n);
        assert_eq!(out.len(), self.n);
        for i in 0..self.n {
            out[i] = mod_mul(a[i], b[i], self.q);
        }
    }

    /// Polynomial multiplication in R_q = Z_q[X] / (X^N + 1) using negacyclic NTT.
    pub fn poly_mul(&self, a: &[i32], b: &[i32]) -> Vec<i32> {
        let mut a_ntt = a.to_vec();
        let mut b_ntt = b.to_vec();
        self.forward(&mut a_ntt);
        self.forward(&mut b_ntt);
        let mut c_ntt = vec![0i32; self.n];
        self.pointwise_mul(&a_ntt, &b_ntt, &mut c_ntt);
        self.inverse(&mut c_ntt);
        c_ntt
    }
}

/// Helper function to reverse bits of index `val` with `bits` width.
pub fn reverse_bits(mut val: usize, bits: usize) -> usize {
    let mut res = 0;
    for _ in 0..bits {
        res = (res << 1) | (val & 1);
        val >>= 1;
    }
    res
}

/// Finds a primitive `2n`-th root of unity modulo prime `q`.
pub fn find_primitive_2n_root(n: usize, q: i32) -> Option<i32> {
    let two_n = (2 * n) as u64;
    let q_minus_1 = (q - 1) as u64;
    if q_minus_1 % two_n != 0 {
        return None;
    }
    let exponent = q_minus_1 / two_n;

    // Search for a generator of Z_q^*
    for g in 2..q {
        let candidate = mod_pow(g as i64, exponent, q);
        // Check if candidate^(2n) == 1 and candidate^n == q - 1 (i.e. -1 mod q)
        if mod_pow(candidate as i64, two_n, q) == 1
            && mod_pow(candidate as i64, n as u64, q) == (q - 1)
        {
            // Verify order is strictly 2n
            let mut is_primitive = true;
            let mut d = 2;
            while d * d <= two_n {
                if two_n % d == 0 {
                    if mod_pow(candidate as i64, d, q) == 1
                        || mod_pow(candidate as i64, two_n / d, q) == 1
                    {
                        is_primitive = false;
                        break;
                    }
                }
                d += 1;
            }
            if is_primitive {
                return Some(candidate);
            }
        }
    }
    None
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::ring::reduction::RLWE_Q;

    #[test]
    fn test_bit_reversal_involution() {
        let bits = 8;
        for i in 0..256 {
            let rev = reverse_bits(i, bits);
            let rev_rev = reverse_bits(rev, bits);
            assert_eq!(i, rev_rev);
        }
    }

    #[test]
    fn test_primitive_roots_properties() {
        let n = DEFAULT_N;
        let q = DEFAULT_Q;
        let psi = find_primitive_2n_root(n, q).expect("Root exists");
        assert_eq!(mod_pow(psi as i64, (2 * n) as u64, q), 1);
        assert_eq!(mod_pow(psi as i64, n as u64, q), q - 1);
    }

    #[test]
    fn test_alternate_ring_lwe_modulus_ntt() {
        // RLWE_Q = 12289, 12288 % 512 == 0
        let ctx = NttContext::new(DEFAULT_N, RLWE_Q);
        let mut data = vec![0i32; DEFAULT_N];
        data[0] = 1;
        data[1] = 2;
        let original = data.clone();
        ctx.forward(&mut data);
        ctx.inverse(&mut data);
        assert_eq!(data, original);
    }
}
