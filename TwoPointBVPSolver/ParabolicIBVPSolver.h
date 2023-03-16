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

// Base class (has pure virtual function!!)
class ParabolicIBVPSolver {
protected:
	ApproxBVP * appr_bvp;
	tridiagonal_matrix * diffusmat;

public:
	ParabolicIBVPSolver(ApproxBVP * prob);

	// A pure virtual function
	virtual vector<double> one_step_march(vector<double> & timeinterval,
		vector<double> & Ul) = 0;

	~ParabolicIBVPSolver();
};

// Forward Euler implementation
class ForwardEuler : public ParabolicIBVPSolver {
public:
	ForwardEuler(ApproxBVP * prob);

	virtual vector<double> one_step_march(vector<double> & timeinterval,
		vector<double> & Ul);
};

class BackwardEuler : public ParabolicIBVPSolver {
public:
	BackwardEuler(ApproxBVP * prob);

	virtual vector<double> one_step_march(vector<double> & timeinterval,
		vector<double> & Ul);
};

class ImplicitMidpoint : public ParabolicIBVPSolver{
	ImplicitMidpoint(ApproxBVP * prob);

virtual vector<double> one_step_march(vector<double> & timeinterval,
	vector<double> & Ul);
};

class ExplicitMidpoint : public ParabolicIBVPSolver {
	ExplicitMidpoint(ApproxBVP * prob);

	virtual vector<double> one_step_march(vector<double> & timeinterval,
		vector<double> & Ul);
};