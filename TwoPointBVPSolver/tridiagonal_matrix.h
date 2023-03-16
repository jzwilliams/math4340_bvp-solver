#include <iostream>
#include <fstream>
#include <vector>
#include <iomanip>
#include <cmath>

using namespace std;

class tridiagonal_matrix {
private:
	int dimension;
	vector <double> diag;
	vector <double> upperdiag;
	vector <double> lowerdiag;
	vector <double> hatupperdiag;			// transformed values for upper diagonal
	vector <double> r;						// vector of denominators for transformation
	bool transformed;

public:
	tridiagonal_matrix(int m);
	tridiagonal_matrix(const tridiagonal_matrix *mat);

	int get_dimension(void) const;
	void set_diagonal_entry(int i, double val);
	void set_upper_diagonal_entry(int i, double val);
	void set_lower_diagonal_entry(int i, double val);

	double get_diagonal_entry(int i) const;
	double get_upper_diagonal_entry(int i) const;
	double get_lower_diagonal_entry(int i) const;

	bool is_transformed() const;

	double get_r_entry(int i) const;
	double get_hatupperdiag_entry(int i) const;

	void add_to_diagonal_entry(int i, double val);								// diag[i] += val
	void add_to_upper_diagonal_entry(int i, double val);						// upperdiag[i] += val
	void add_to_lower_diagonal_entry(int i, double val);						// lowerdiag[i] += val

	void scalarmult_diagonal_entry(int i, double val);
	void scalarmult_upper_diagonal_entry(int i, double val);
	void scalarmult_lower_diagonal_entry(int i, double val);

	vector <double> solve_linear_system(const vector<double> & rhs);

	void transform(void);

	vector <double> tridiagonal_multiply(const vector <double> & lhs);				// matrix multiplication

	~tridiagonal_matrix();
};