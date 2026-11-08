//! Module-LWE polynomial vectors and matrices over R_q = Z_q[X] / (X^N + 1).

use std::fmt;
use std::ops::{Add, AddAssign, Index, IndexMut, Sub, SubAssign};
use serde::{Deserialize, Serialize};

use crate::ring::ntt::NttContext;
use crate::ring::poly::Poly;

/// A vector of polynomials in R_q^k.
#[derive(Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct PolyVector {
    pub polys: Vec<Poly>,
    pub q: i32,
}

impl PolyVector {
    /// Creates a zero vector of `dim` polynomials of degree `n` modulo `q`.
    pub fn zero(dim: usize, n: usize, q: i32) -> Self {
        Self {
            polys: vec![Poly::zero(n, q); dim],
            q,
        }
    }

    /// Creates a polynomial vector from a slice of polynomials.
    pub fn new(polys: Vec<Poly>) -> Self {
        assert!(!polys.is_empty(), "Vector cannot be empty");
        let q = polys[0].q;
        let n = polys[0].len();
        for p in &polys {
            assert_eq!(p.q, q, "All polynomials must share the same modulus");
            assert_eq!(p.len(), n, "All polynomials must have the same degree");
        }
        Self { polys, q }
    }

    /// Dimension (number of polynomial elements) in the vector.
    pub fn dim(&self) -> usize {
        self.polys.len()
    }

    /// Degree of polynomial elements.
    pub fn poly_len(&self) -> usize {
        if self.polys.is_empty() {
            0
        } else {
            self.polys[0].len()
        }
    }

    /// Infinity norm of the vector: max_i ||p_i||_infinity.
    pub fn infinity_norm(&self) -> i32 {
        self.polys.iter().map(|p| p.infinity_norm()).max().unwrap_or(0)
    }

    /// Squared Euclidean norm ||v||_2^2: sum of squared norms of all polynomials.
    pub fn l2_norm_sq(&self) -> i64 {
        self.polys.iter().map(|p| p.l2_norm_sq()).sum()
    }

    /// Euclidean L2 norm ||v||_2.
    pub fn l2_norm(&self) -> f64 {
        (self.l2_norm_sq() as f64).sqrt()
    }

    /// Multiplies the entire vector by a scalar polynomial in R_q using NTT.
    pub fn poly_scalar_mul(&self, scalar_poly: &Poly, ctx: &NttContext) -> Self {
        let scalar_ntt = scalar_poly.to_ntt(ctx);
        let polys: Vec<Poly> = self
            .polys
            .iter()
            .map(|p| {
                let p_ntt = p.to_ntt(ctx);
                let prod_ntt = p_ntt.pointwise_mul(&scalar_ntt);
                prod_ntt.to_poly(ctx)
            })
            .collect();
        Self { polys, q: self.q }
    }

    /// Multiplies the vector by an integer scalar.
    pub fn scalar_mul(&self, scalar: i32) -> Self {
        let polys: Vec<Poly> = self.polys.iter().map(|p| p.scalar_mul(scalar)).collect();
        Self { polys, q: self.q }
    }

    /// Checks if all polynomials in the vector are zero.
    pub fn is_zero(&self) -> bool {
        self.polys.iter().all(|p| p.is_zero())
    }
}

impl Index<usize> for PolyVector {
    type Output = Poly;
    fn index(&self, index: usize) -> &Self::Output {
        &self.polys[index]
    }
}

impl IndexMut<usize> for PolyVector {
    fn index_mut(&mut self, index: usize) -> &mut Self::Output {
        &mut self.polys[index]
    }
}

impl Add for &PolyVector {
    type Output = PolyVector;
    fn add(self, rhs: Self) -> Self::Output {
        assert_eq!(self.dim(), rhs.dim(), "Vector dimensions must match");
        assert_eq!(self.q, rhs.q, "Moduli must match");
        let polys: Vec<Poly> = self
            .polys
            .iter()
            .zip(&rhs.polys)
            .map(|(a, b)| a + b)
            .collect();
        PolyVector { polys, q: self.q }
    }
}

impl Add for PolyVector {
    type Output = PolyVector;
    fn add(self, rhs: Self) -> Self::Output {
        &self + &rhs
    }
}

impl AddAssign<&PolyVector> for PolyVector {
    fn add_assign(&mut self, rhs: &PolyVector) {
        assert_eq!(self.dim(), rhs.dim());
        for i in 0..self.dim() {
            self.polys[i] += &rhs.polys[i];
        }
    }
}

impl Sub for &PolyVector {
    type Output = PolyVector;
    fn sub(self, rhs: Self) -> Self::Output {
        assert_eq!(self.dim(), rhs.dim(), "Vector dimensions must match");
        assert_eq!(self.q, rhs.q, "Moduli must match");
        let polys: Vec<Poly> = self
            .polys
            .iter()
            .zip(&rhs.polys)
            .map(|(a, b)| a - b)
            .collect();
        PolyVector { polys, q: self.q }
    }
}

impl Sub for PolyVector {
    type Output = PolyVector;
    fn sub(self, rhs: Self) -> Self::Output {
        &self - &rhs
    }
}

impl SubAssign<&PolyVector> for PolyVector {
    fn sub_assign(&mut self, rhs: &PolyVector) {
        assert_eq!(self.dim(), rhs.dim());
        for i in 0..self.dim() {
            self.polys[i] -= &rhs.polys[i];
        }
    }
}

impl fmt::Debug for PolyVector {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "PolyVector(dim={}, poly_n={}, norm_inf={})",
            self.dim(),
            self.poly_len(),
            self.infinity_norm()
        )
    }
}

/// A matrix of polynomials in R_q^(k x l).
/// Stored in row-major order: rows = k, cols = l.
#[derive(Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct PolyMatrix {
    pub rows: usize,
    pub cols: usize,
    pub data: Vec<Vec<Poly>>, // [rows][cols]
    pub q: i32,
}

impl PolyMatrix {
    /// Creates a zero polynomial matrix of size `rows x cols`.
    pub fn zero(rows: usize, cols: usize, n: usize, q: i32) -> Self {
        let data = vec![vec![Poly::zero(n, q); cols]; rows];
        Self { rows, cols, data, q }
    }

    /// Constructs a matrix from a 2D grid of polynomials.
    pub fn new(data: Vec<Vec<Poly>>) -> Self {
        assert!(!data.is_empty(), "Matrix rows cannot be empty");
        let rows = data.len();
        let cols = data[0].len();
        assert!(cols > 0, "Matrix columns cannot be empty");
        let q = data[0][0].q;
        let n = data[0][0].len();

        for row in &data {
            assert_eq!(row.len(), cols, "All rows must have the same column count");
            for p in row {
                assert_eq!(p.q, q, "All polynomials must share the same modulus");
                assert_eq!(p.len(), n, "All polynomials must have the same degree");
            }
        }

        Self { rows, cols, data, q }
    }

    /// Matrix-vector multiplication: `A * v mod (X^N + 1, q)`
    /// Computes result in R_q^k where `A` is k x l and `v` is l x 1.
    pub fn mul_vector(&self, v: &PolyVector, ctx: &NttContext) -> PolyVector {
        assert_eq!(self.cols, v.dim(), "Matrix column count must match vector dimension");
        let n = ctx.n;
        let q = self.q;

        // Transform input vector to NTT domain once
        let v_ntt: Vec<_> = v.polys.iter().map(|p| p.to_ntt(ctx)).collect();

        let mut res_polys = Vec::with_capacity(self.rows);

        for row_idx in 0..self.rows {
            let mut acc_ntt = vec![0i32; n];
            let mut temp_prod = vec![0i32; n];

            for col_idx in 0..self.cols {
                let a_poly = &self.data[row_idx][col_idx];
                let a_ntt = a_poly.to_ntt(ctx);
                ctx.pointwise_mul(&a_ntt.evals, &v_ntt[col_idx].evals, &mut temp_prod);

                for i in 0..n {
                    acc_ntt[i] = crate::ring::reduction::mod_add(acc_ntt[i], temp_prod[i], q);
                }
            }

            ctx.inverse(&mut acc_ntt);
            res_polys.push(Poly {
                coeffs: acc_ntt,
                q,
            });
        }

        PolyVector {
            polys: res_polys,
            q,
        }
    }
}

impl fmt::Debug for PolyMatrix {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "PolyMatrix({}x{}, poly_len={})",
            self.rows,
            self.cols,
            if self.rows > 0 && self.cols > 0 { self.data[0][0].len() } else { 0 }
        )
    }
}
