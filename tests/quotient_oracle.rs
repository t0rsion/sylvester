//! Independent determinant and minimal polynomial checks for finite quotients.

use num_rational::BigRational;
use num_traits::{One, Zero};
use sylvester::{
    Budget, Felt, FiniteQuotient, GroebnerBasis, MultiplicationMatrix, PolynomialRing, PrimeField,
    Rationals, UnivariatePolynomial,
};

fn prime_quotient(
    modulus: u64,
    variables: &[&str],
    generators: &[&str],
) -> (PolynomialRing, FiniteQuotient) {
    let ring = PolynomialRing::prime_field(modulus, variables.iter().copied())
        .expect("the prime ring builds");
    let polynomials = generators
        .iter()
        .map(|text| ring.parse_polynomial(text).expect("the generator parses"))
        .collect();
    let basis = GroebnerBasis::<PrimeField>::from_polynomials(&ring, polynomials, Budget::new())
        .expect("the supplied list is a reduced Groebner basis");
    let quotient = basis
        .finite_quotient(Budget::new())
        .expect("the quotient is finite");
    (ring, quotient)
}

fn rational_quotient(
    variables: &[&str],
    generators: &[&str],
) -> (PolynomialRing<Rationals>, FiniteQuotient<Rationals>) {
    let ring =
        PolynomialRing::rationals(variables.iter().copied()).expect("the rational ring builds");
    let polynomials = generators
        .iter()
        .map(|text| ring.parse_polynomial(text).expect("the generator parses"))
        .collect();
    let basis = GroebnerBasis::<Rationals>::from_polynomials(&ring, polynomials, Budget::new())
        .expect("the supplied list is a reduced Groebner basis");
    let quotient = basis
        .finite_quotient(Budget::new())
        .expect("the quotient is finite");
    (ring, quotient)
}

fn coefficient_values(polynomial: &UnivariatePolynomial, modulus: u64) -> Vec<u64> {
    polynomial
        .coefficients()
        .iter()
        .map(|coefficient| coefficient.map(Felt::value).unwrap_or(0) % modulus)
        .collect()
}

fn matrix_value(matrix: &MultiplicationMatrix, row: usize, column: usize) -> u64 {
    matrix
        .entry(row, column)
        .map(|coefficient| coefficient.value())
        .unwrap_or(0)
}

fn add_polynomial(target: &mut [u64], source: &[u64], sign: i64, modulus: u64) {
    for (left, right) in target.iter_mut().zip(source) {
        let value = if sign < 0 {
            (modulus - *right % modulus) % modulus
        } else {
            *right % modulus
        };
        *left = (*left + value) % modulus;
    }
}

fn multiply_polynomials(left: &[u64], right: &[u64], modulus: u64) -> Vec<u64> {
    let mut product = vec![0; left.len() + right.len() - 1];
    for (left_degree, left_value) in left.iter().enumerate() {
        for (right_degree, right_value) in right.iter().enumerate() {
            product[left_degree + right_degree] =
                (product[left_degree + right_degree] + left_value * right_value) % modulus;
        }
    }
    product
}

fn permutation_sign(permutation: &[usize]) -> i64 {
    let mut inversions: usize = 0;
    for left in 0..permutation.len() {
        for right in left + 1..permutation.len() {
            if permutation[left] > permutation[right] {
                inversions += 1;
            }
        }
    }
    if inversions.is_multiple_of(2) { 1 } else { -1 }
}

fn visit_permutations<F>(permutation: &mut [usize], start: usize, visit: &mut F)
where
    F: FnMut(&[usize]),
{
    if start == permutation.len() {
        visit(permutation);
        return;
    }
    for index in start..permutation.len() {
        permutation.swap(start, index);
        visit_permutations(permutation, start + 1, visit);
        permutation.swap(start, index);
    }
}

/// Return `det(tI - A)` by the Leibniz formula.
fn determinant_polynomial(matrix: &MultiplicationMatrix, modulus: u64) -> Vec<u64> {
    let dimension = matrix.dimension();
    let mut result = vec![0; dimension + 1];
    let mut permutation: Vec<usize> = (0..dimension).collect();
    visit_permutations(&mut permutation, 0, &mut |permutation| {
        let mut term = vec![1];
        for (row, &column) in permutation.iter().enumerate() {
            let value = matrix_value(matrix, row, column);
            let factor = if row == column {
                vec![(modulus - value) % modulus, 1]
            } else {
                vec![(modulus - value) % modulus]
            };
            term = multiply_polynomials(&term, &factor, modulus);
        }
        add_polynomial(&mut result, &term, permutation_sign(permutation), modulus);
    });
    result
}

fn matrix_vector(matrix: &MultiplicationMatrix, vector: &[u64], modulus: u64) -> Vec<u64> {
    (0..matrix.dimension())
        .map(|row| {
            (0..matrix.dimension()).fold(0, |sum, column| {
                (sum + matrix_value(matrix, row, column) * vector[column]) % modulus
            })
        })
        .collect()
}

fn horner_matrix_column(
    matrix: &MultiplicationMatrix,
    coefficients: &[u64],
    column: usize,
    modulus: u64,
) -> Vec<u64> {
    let mut value = vec![0; matrix.dimension()];
    for degree in (0..coefficients.len()).rev() {
        value = matrix_vector(matrix, &value, modulus);
        value[column] = (value[column] + coefficients[degree]) % modulus;
    }
    value
}

fn inverse_mod(mut value: u64, modulus: u64) -> u64 {
    let mut exponent = modulus - 2;
    let mut result = 1;
    while exponent != 0 {
        if exponent & 1 == 1 {
            result = result * value % modulus;
        }
        value = value * value % modulus;
        exponent >>= 1;
    }
    result
}

fn rank_mod(columns: &[Vec<u64>], modulus: u64) -> usize {
    if columns.is_empty() {
        return 0;
    }
    let rows = columns[0].len();
    let count = columns.len();
    let mut work = vec![vec![0; count]; rows];
    for row in 0..rows {
        for column in 0..count {
            work[row][column] = columns[column][row] % modulus;
        }
    }

    let mut rank = 0;
    for column in 0..count {
        let Some(pivot) = (rank..rows).find(|&row| work[row][column] != 0) else {
            continue;
        };
        work.swap(rank, pivot);
        let inverse = inverse_mod(work[rank][column], modulus);
        for row in rank + 1..rows {
            if work[row][column] == 0 {
                continue;
            }
            let factor = work[row][column] * inverse % modulus;
            let mut right = column;
            while right < count {
                let subtraction = factor * work[rank][right] % modulus;
                work[row][right] = (work[row][right] + modulus - subtraction) % modulus;
                right += 1;
            }
        }
        rank += 1;
    }
    rank
}

fn krylov_columns(matrix: &MultiplicationMatrix, degree: usize, modulus: u64) -> Vec<Vec<u64>> {
    let mut current = vec![0; matrix.dimension()];
    if !current.is_empty() {
        current[0] = 1;
    }
    let mut columns = Vec::with_capacity(degree + 1);
    for _ in 0..=degree {
        columns.push(current.clone());
        current = matrix_vector(matrix, &current, modulus);
    }
    columns
}

fn evaluate_on_one(coefficients: &[u64], columns: &[Vec<u64>], modulus: u64) -> Vec<u64> {
    let mut result = vec![0; columns[0].len()];
    for (coefficient, vector) in coefficients.iter().zip(columns) {
        for (value, term) in result.iter_mut().zip(vector) {
            *value = (*value + coefficient * term) % modulus;
        }
    }
    result
}

fn assert_prime_minimality(
    matrix: &MultiplicationMatrix,
    polynomial: &UnivariatePolynomial,
    modulus: u64,
) {
    let coefficients = coefficient_values(polynomial, modulus);
    let degree = polynomial
        .degree()
        .expect("a nonzero minimal polynomial has a degree");
    assert_eq!(coefficients[degree], 1, "the minimal polynomial is monic");
    let columns = krylov_columns(matrix, degree, modulus);
    for prefix in 0..=degree {
        assert_eq!(rank_mod(&columns[..prefix], modulus), prefix);
    }
    assert_eq!(rank_mod(&columns, modulus), degree);
    assert_eq!(
        evaluate_on_one(&coefficients, &columns, modulus),
        vec![0; matrix.dimension()]
    );
    for column in 0..matrix.dimension() {
        assert_eq!(
            horner_matrix_column(matrix, &coefficients, column, modulus),
            vec![0; matrix.dimension()]
        );
    }
}

fn rational_matrix_vector(
    matrix: &MultiplicationMatrix<Rationals>,
    vector: &[BigRational],
) -> Vec<BigRational> {
    (0..matrix.dimension())
        .map(|row| {
            (0..matrix.dimension()).fold(BigRational::zero(), |sum, column| {
                let entry = matrix
                    .entry(row, column)
                    .cloned()
                    .unwrap_or_else(BigRational::zero);
                sum + entry * vector[column].clone()
            })
        })
        .collect()
}

fn rational_horner_matrix_column(
    matrix: &MultiplicationMatrix<Rationals>,
    coefficients: &[Option<BigRational>],
    column: usize,
) -> Vec<BigRational> {
    let mut value = vec![BigRational::zero(); matrix.dimension()];
    for degree in (0..coefficients.len()).rev() {
        value = rational_matrix_vector(matrix, &value);
        if let Some(coefficient) = &coefficients[degree] {
            value[column] += coefficient.clone();
        }
    }
    value
}

fn rational_rank(columns: &[Vec<BigRational>]) -> usize {
    if columns.is_empty() {
        return 0;
    }
    let rows = columns[0].len();
    let count = columns.len();
    let mut work = vec![vec![BigRational::zero(); count]; rows];
    for row in 0..rows {
        for column in 0..count {
            work[row][column] = columns[column][row].clone();
        }
    }
    let mut rank = 0;
    for column in 0..count {
        let Some(pivot) = (rank..rows).find(|&row| !work[row][column].is_zero()) else {
            continue;
        };
        work.swap(rank, pivot);
        let pivot_value = work[rank][column].clone();
        for row in rank + 1..rows {
            if work[row][column].is_zero() {
                continue;
            }
            let factor = work[row][column].clone() / pivot_value.clone();
            let mut right = column;
            while right < count {
                work[row][right] =
                    work[row][right].clone() - factor.clone() * work[rank][right].clone();
                right += 1;
            }
        }
        rank += 1;
    }
    rank
}

#[test]
fn characteristic_polynomial_matches_an_independent_oracle_over_small_fields() {
    for modulus in [2, 3, 7] {
        for degree in 1..=6 {
            for seed in 0..5 {
                let ring = PolynomialRing::prime_field(modulus, ["x"]).expect("the ring builds");
                let mut terms = vec![(1, vec![degree as u16])];
                for exponent in 0..degree {
                    let coefficient = (seed + 3 * exponent + degree) % modulus as usize;
                    if coefficient != 0 {
                        terms.push((coefficient as u64, vec![exponent as u16]));
                    }
                }
                let generator = ring.polynomial(terms).expect("the polynomial builds");
                let basis = GroebnerBasis::<PrimeField>::from_polynomials(
                    &ring,
                    vec![generator],
                    Budget::new(),
                )
                .expect("the univariate basis is reduced");
                let quotient = basis
                    .finite_quotient(Budget::new())
                    .expect("the quotient is finite");
                let x = ring.generator(0).expect("the ring has x");
                let matrix = quotient
                    .multiplication_matrix(&x, Budget::new())
                    .expect("the matrix computes");
                let expected = determinant_polynomial(&matrix, modulus);
                let mut expected_generator = vec![0; degree + 1];
                for (exponent, coefficient) in
                    expected_generator.iter_mut().enumerate().take(degree)
                {
                    *coefficient = (seed + 3 * exponent + degree) as u64 % modulus;
                }
                expected_generator[degree] = 1;
                assert_eq!(expected, expected_generator);
                let characteristic = quotient
                    .characteristic_polynomial(&x, Budget::new())
                    .expect("the characteristic polynomial computes");
                assert_eq!(coefficient_values(&characteristic, modulus), expected);
                let minimal = quotient
                    .minimal_polynomial(&x, Budget::new())
                    .expect("the minimal polynomial computes");
                assert_eq!(coefficient_values(&minimal, modulus), expected_generator);
                assert_prime_minimality(&matrix, &minimal, modulus);
            }
        }
    }
}

#[test]
fn multivariate_nilpotent_and_separable_minimal_polynomials_pass_rank_checks() {
    let (ring, nilpotent) = prime_quotient(2, &["x", "y"], &["x^2", "y^2"]);
    let x = ring.parse_polynomial("x").expect("x parses");
    let matrix = nilpotent
        .multiplication_matrix(&x, Budget::new())
        .expect("the nilpotent matrix computes");
    let minimal = nilpotent
        .minimal_polynomial(&x, Budget::new())
        .expect("the nilpotent minimal polynomial computes");
    let characteristic = nilpotent
        .characteristic_polynomial(&x, Budget::new())
        .expect("the nilpotent characteristic polynomial computes");
    assert_eq!(
        coefficient_values(&characteristic, 2),
        determinant_polynomial(&matrix, 2)
    );
    assert_eq!(
        matrix
            .entries()
            .iter()
            .map(|value| value.map(Felt::value))
            .collect::<Vec<_>>(),
        vec![
            None,
            None,
            None,
            None,
            None,
            None,
            None,
            None,
            Some(1),
            None,
            None,
            None,
            None,
            Some(1),
            None,
            None,
        ]
    );
    assert_prime_minimality(&matrix, &minimal, 2);
    assert_eq!(coefficient_values(&minimal, 2), vec![0, 0, 1]);

    let (ring, separable) = prime_quotient(7, &["x", "y"], &["x^2 - 1", "y^2 - 1"]);
    let element = ring.parse_polynomial("x + y").expect("x + y parses");
    let matrix = separable
        .multiplication_matrix(&element, Budget::new())
        .expect("the separable matrix computes");
    let minimal = separable
        .minimal_polynomial(&element, Budget::new())
        .expect("the separable minimal polynomial computes");
    let characteristic = separable
        .characteristic_polynomial(&element, Budget::new())
        .expect("the separable characteristic polynomial computes");
    assert_eq!(
        coefficient_values(&characteristic, 7),
        determinant_polynomial(&matrix, 7)
    );
    assert_eq!(
        matrix
            .entries()
            .iter()
            .map(|value| value.map(Felt::value))
            .collect::<Vec<_>>(),
        vec![
            None,
            Some(1),
            Some(1),
            None,
            Some(1),
            None,
            None,
            Some(1),
            Some(1),
            None,
            None,
            Some(1),
            None,
            Some(1),
            Some(1),
            None,
        ]
    );
    assert_prime_minimality(&matrix, &minimal, 7);
    assert_eq!(coefficient_values(&minimal, 7), vec![0, 3, 0, 1]);
}

#[test]
fn determinant_oracle_handles_unit_and_zero_variable_quotients() {
    let (ring, unit) = prime_quotient(7, &["x"], &["1"]);
    let x = ring.parse_polynomial("x").expect("x parses");
    let matrix = unit
        .multiplication_matrix(&x, Budget::new())
        .expect("the unit matrix computes");
    assert_eq!(determinant_polynomial(&matrix, 7), vec![1]);
    let characteristic = unit
        .characteristic_polynomial(&x, Budget::new())
        .expect("the unit characteristic polynomial computes");
    let minimal = unit
        .minimal_polynomial(&x, Budget::new())
        .expect("the unit minimal polynomial computes");
    assert_eq!(coefficient_values(&characteristic, 7), vec![1]);
    assert_eq!(coefficient_values(&minimal, 7), vec![1]);

    let (ring, constants) = prime_quotient(2, &[], &[]);
    let one = ring.one();
    let matrix = constants
        .multiplication_matrix(&one, Budget::new())
        .expect("the constant matrix computes");
    assert_eq!(determinant_polynomial(&matrix, 2), vec![1, 1]);
    let characteristic = constants
        .characteristic_polynomial(&one, Budget::new())
        .expect("the constant characteristic polynomial computes");
    let minimal = constants
        .minimal_polynomial(&one, Budget::new())
        .expect("the constant minimal polynomial computes");
    assert_eq!(coefficient_values(&characteristic, 2), vec![1, 1]);
    assert_eq!(coefficient_values(&minimal, 2), vec![1, 1]);
}

#[test]
fn rational_nilpotent_minimal_polynomial_has_the_first_rank_dependence() {
    let (ring, quotient) = rational_quotient(&["x", "y"], &["x^2", "y^2"]);
    let x = ring.parse_polynomial("x").expect("x parses");
    let matrix = quotient
        .multiplication_matrix(&x, Budget::new())
        .expect("the rational matrix computes");
    assert_eq!(matrix.dimension(), 4);
    let minimal = quotient
        .minimal_polynomial(&x, Budget::new())
        .expect("the rational minimal polynomial computes");
    let degree = minimal.degree().expect("the polynomial has a degree");
    assert_eq!(degree, 2);
    assert!(minimal.coefficient(degree).is_some_and(BigRational::is_one));
    assert_rational_minimality(&matrix, &minimal, degree);
}

fn assert_rational_minimality(
    matrix: &MultiplicationMatrix<Rationals>,
    minimal: &UnivariatePolynomial<Rationals>,
    degree: usize,
) {
    let mut current = vec![BigRational::zero(); matrix.dimension()];
    current[0] = BigRational::one();
    let mut columns = Vec::with_capacity(degree + 1);
    for _ in 0..=degree {
        columns.push(current.clone());
        current = rational_matrix_vector(matrix, &current);
    }
    for prefix in 0..=degree {
        assert_eq!(rational_rank(&columns[..prefix]), prefix);
    }
    assert_eq!(rational_rank(&columns), degree);
    let result = rational_linear_combination(minimal.coefficients(), &columns, matrix.dimension());
    assert!(result.iter().all(BigRational::is_zero));
    for column in 0..matrix.dimension() {
        assert!(
            rational_horner_matrix_column(matrix, minimal.coefficients(), column)
                .iter()
                .all(BigRational::is_zero)
        );
    }
}

fn rational_linear_combination(
    coefficients: &[Option<BigRational>],
    vectors: &[Vec<BigRational>],
    dimension: usize,
) -> Vec<BigRational> {
    let mut result = vec![BigRational::zero(); dimension];
    for (coefficient, vector) in coefficients.iter().zip(vectors) {
        if let Some(coefficient) = coefficient {
            for (value, term) in result.iter_mut().zip(vector) {
                *value += coefficient.clone() * term.clone();
            }
        }
    }
    result
}
