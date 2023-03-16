#include <iostream>
#include <fstream>
#include <vector>
#include <iomanip>
#include <cmath>


#ifndef TRIDIAGONAL_MATRIX_H
#define TRIDIAGONAL_MATRIX_H
#include "tridiagonal_matrix.h"
#endif

#ifndef TWOPOINTBVP_H
#define TWOPOINTBVP_H
#include "twopointbvp.h"
#endif

#ifndef APPROX_BVP_H
#define APPROX_BVP_H
#include "approx_bvp.h"
#endif

using namespace std;

// Base class
class HyperbolicIBVPSolver {
protected:
	ApproxBVP * appr_bvp;
	tridiagonal_matrix * diffusmat;
	const TwoPointBVP * theProblem;

public:
	HyperbolicIBVPSolver(ApproxBVP * prob);

	virtual vector<double> two_step_march(vector<double> &timeinterval,
		vector<double> &Unm1, vector<double> Un) = 0;

	virtual vector<double> first_march(vector<double> &timeinterval,
		vector<double> U0) = 0;

	~HyperbolicIBVPSolver();
};

// Two-step implicit right-point rule implementation
class TwoStepImplicitRPR : public HyperbolicIBVPSolver {
public:
	TwoStepImplicitRPR(ApproxBVP * prob);

	virtual vector<double> two_step_march(vector<double> &timeinterval,
		vector<double> &Unm1, vector<double> &Un);

	virtual vector<double> first_march(vector<double> &timeinterval,
		vector<double> &U0, vector<double> & didt, vector<double> & d2idt2);
};