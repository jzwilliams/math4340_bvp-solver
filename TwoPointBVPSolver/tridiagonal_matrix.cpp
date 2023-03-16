#ifndef TRIDIAGONAL_MATRIX_H
#define TRIDIAGONAL_MATRIX_H
#include "tridiagonal_matrix.h"
#endif

tridiagonal_matrix::tridiagonal_matrix(int m) {
	dimension = m;
	diag.resize(m);
	upperdiag.resize(m - 1);
	lowerdiag.resize(m - 1);
	hatupperdiag.resize(m - 1);
	r.resize(m - 1);
	transformed = false;
}

// Copy constructor
tridiagonal_matrix::tridiagonal_matrix(const tridiagonal_matrix *mat) {
	dimension = mat->get_dimension();
	diag.resize(dimension);
	upperdiag.resize(dimension - 1);
	lowerdiag.resize(dimension - 1);
	hatupperdiag.resize(dimension - 1);
	r.resize(dimension - 1);
	for (int i = 0; i < dimension; i++) {
		diag[i] = mat->get_diagonal_entry(i);
	}
	for (int i = 0; i < dimension - 1; i++) {
		lowerdiag[i] = mat->get_lower_diagonal_entry(i);
		upperdiag[i] = mat->get_upper_diagonal_entry(i);
	}

	transformed = mat->is_transformed();

	if (transformed) {
		for (int i = 0; i < dimension - 1; i++) {
			r.push_back(mat->get_r_entry(i));
			hatupperdiag.push_back(mat->get_hatupperdiag_entry(i));
		}

		// transformed = true; Unnecessary
	}

}

int tridiagonal_matrix::get_dimension() const {
	return dimension;
}

void tridiagonal_matrix::set_diagonal_entry(int i, double val) {
	diag[i] = val;
}

void tridiagonal_matrix::set_upper_diagonal_entry(int i, double val) {
	upperdiag[i] = val;
}

void tridiagonal_matrix::set_lower_diagonal_entry(int i, double val) {
	lowerdiag[i] = val;
}

double tridiagonal_matrix::get_diagonal_entry(int i) const {
	return diag[i];
}

double tridiagonal_matrix::get_upper_diagonal_entry(int i) const {
	return upperdiag[i];
}

double tridiagonal_matrix::get_lower_diagonal_entry(int i) const {
	return lowerdiag[i];
}

bool tridiagonal_matrix::is_transformed() const {
	return transformed;
}

double tridiagonal_matrix::get_r_entry(int i) const {
	return r[i];
}

double tridiagonal_matrix::get_hatupperdiag_entry(int i) const {
	return hatupperdiag[i];
}

void tridiagonal_matrix::add_to_diagonal_entry(int i, double val) {
	diag[i] += val;
}

void tridiagonal_matrix::add_to_upper_diagonal_entry(int i, double val) {
	upperdiag[i] += val;
}

void tridiagonal_matrix::add_to_lower_diagonal_entry(int i, double val) {
	lowerdiag[i] += val;
}

void tridiagonal_matrix::scalarmult_diagonal_entry(int i, double val) {
	diag[i] *= val;
}

void tridiagonal_matrix::scalarmult_upper_diagonal_entry(int i, double val) {
	upperdiag[i] *= val;
}

void tridiagonal_matrix::scalarmult_lower_diagonal_entry(int i, double val) {
	lowerdiag[i] *= val;
}

void tridiagonal_matrix::transform() {
	// Must implement the transformation using r and hatupperdiag
	int m = dimension;
	hatupperdiag[0] = upperdiag[0] / diag[0];
	for (int i = 1; i <= (m - 2); i++) {
		r[i - 1] = diag[i] - (lowerdiag[i - 1] * hatupperdiag[i - 1]);
		hatupperdiag[i] = upperdiag[i] / r[i - 1];
	}
	r[m - 2] = diag[m - 1] - (lowerdiag[m - 2] * hatupperdiag[m - 2]);

	transformed = true;
}

// Given RHS vector, return sol satisfying A*sol = RHS
vector <double> tridiagonal_matrix::solve_linear_system(const vector <double> & rhs) {
	int m = dimension;
	vector <double> hatrhs(dimension,0.0);			// the modified RHS
	vector <double> sol(dimension,0.0);				// the solution

	if (!transformed) {
		transform();
	}

	// Transform RHS vector
	hatrhs[0] = rhs[0] / diag[0];
	for (int i = 1; i <= m - 2; i++) {
		hatrhs[i] = (rhs[i] - (lowerdiag[i - 1] * hatrhs[i - 1])) / r[i - 1];
	}
	hatrhs[m - 1] = (rhs[m - 1] - (lowerdiag[m - 2] * hatrhs[m - 2])) / r[m - 2];

	// Solve the system
	sol[m - 1] = hatrhs[m - 1];
	for (int i = (m - 2); i >= 0; i--) {
		sol[i] = hatrhs[i] - (hatupperdiag[i] * sol[i + 1]);
	}

	return sol;
}

// Given vector LHS, return RHS = A*LHS
vector <double> tridiagonal_matrix::tridiagonal_multiply(const vector <double> & lhs) {
	int m = dimension;
	vector <double> rhs;
	rhs.resize(m);
	rhs[0] = diag[0] * lhs[0] + upperdiag[0] * lhs[1];
	for (int i = 1; i <= m - 2; i++) {
		rhs[i] = lowerdiag[i - 1] * lhs[i - 1] + diag[i] * lhs[i] + upperdiag[i] * lhs[i + 1];
	}
	rhs[m - 1] = lowerdiag[m - 2] * lhs[m - 2] + diag[m - 1] * lhs[m - 1];

	return rhs;
}

tridiagonal_matrix::~tridiagonal_matrix() {		// Destructor can be empty
	;
}