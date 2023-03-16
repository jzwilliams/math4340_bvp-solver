#include <iostream>
#include <fstream>
#include <cmath>
#include <iomanip>
#include <vector>

#ifndef TRIDIAGONAL_MATRIX_H
#define TRIDIAGONAL_MATRIX_H
#include "tridiagonal_matrix.h"
#endif

#ifndef TWOPOINTBVP_H
#define TWOPOINTBVP_H
#include "twopointbvp.h"
#endif

using namespace std;

class ApproxBVP {
protected:
	int numSubInt;
	const double * stepLengths;
	const TwoPointBVP * theProblem;
	vector <double> xcoord;
	vector <double> midcoord;
	vector <double> Deltax;

	void AssembleDiffusion(tridiagonal_matrix * tmat);

	void AssembleReaction(vector <double> & W, vector <double> & R, vector <double> & prpu);
	vector <double> AssembleForce();

public:
	ApproxBVP(int N, const double * subIntervalLengths, 
		const TwoPointBVP * prob);

	// Time-derivative functions
	tridiagonal_matrix* calcDiffusion(void);
	vector<double> calcLumpedMass(void);
	void calcReaction(vector<double>&, vector<double>&, vector<double>&);
	vector<double> calcForce(double timeLevel);

	// Two-point BVP accessor
	const TwoPointBVP * getBVP(void) const;

	vector <double> getXCoord(void);

	vector <double> solve(int maxIter, double tol);

	// Basic functions
	bool left_bound_Dirichlet(void);
	bool right_bound_Dirichlet(void);

	~ApproxBVP();
};