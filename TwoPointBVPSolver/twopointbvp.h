#include <iostream>
#include <cmath>
#include <iomanip>
#include <vector>

// Jacob Williams
// MATH 4340
// A C++ class which defines a two-point boundary value problem

using namespace std;

class TwoPointBVP {
protected: // outside class can't access these directly
	// These are: the members of the class (properties)
	double * domain;
	bool leftBoundIsDirichlet;		// TRUE if left boundary condition is Dirichlet
	bool rightBoundIsDirichlet;		// TRUE if right boundary condition is Dirichlet
	double * leftBoundVals;			// contains gamma_a and g_a
	double * rightBoundVals;		// contains gamma_b and g_b
		// if gamma_0 == 0 => Neumann boundary; if nonzero => Robin boundary
	bool containsReaction;			// TRUE if has reaction coefficient
	bool containsExtForce;			// TRUE if has external force
	double (*diffusion)(vector<double> &);	// diffusion coefficient k(x)
	double(*partialdiffusion)(vector<double> &);	// derivative of diffusion coefficient
	double(*reaction)(vector<double> &);	// reaction coefficient r(x)
	double(*extForce)(vector<double> &);	// external force f(x)
	double(*partialreactionpartialu)(vector<double> &);	// Jacobian of external force wrt u

	//Time-dependent stuff
	double gammaa;
	double gammab;

	double(*ga)(vector<double>&);
	double(*gb)(vector<double>&);

public:
	// Constructor for regular BVP or parabolic IBVP
	TwoPointBVP(double *dom, double(*func)(vector<double> &));
		// Domain, vector for diffusion coefficient

	// Constructor for hyperbolic IBVP
	TwoPointBVP(double *dom, double(*func)(vector<double> &), double(*pdfunc)(vector<double> &));
		// Domain, vector for diffusion coefficient, vector for derivative of diffusion coefficient

	// Set the boundaries
	void set_left_bdry(bool _leftIsDirichlet, double *val);
	void set_right_bdry(bool _rightIsDirichlet, double *val);

	// Overloaded boundaries for time-dependent problems
	void set_left_bdry(bool _leftIsDirichlet, double gammaa, double(*gleft)(vector<double>&));
	void set_right_bdry(bool _rightIsDirichlet, double gammab, double(*gright)(vector<double> &par));
	vector<double> calcLeftBdry(double t) const; // Evaluates left boundary value at certain time
	vector<double> calcRightBdry(double t) const; // Likewise for right boundary

	// Set the reaction (automatically sets containsReaction to FALSE if not called)
	void set_reaction(double(*funcone) (vector<double> &),
		double(*functwo) (vector<double> &));

	// Set the external force
	void set_ext_force(double(*func)(vector<double> &));

	// Accessor functions:

	// Return domain: see above
	double *get_domain() const;

	// Return whether boundaries are Dirichlet
	bool left_bdry_is_Dirichlet() const;
	bool right_bdry_is_Dirichlet() const;


	// Return boundary values; one contains (gamma_a, g_a), other contains (gamma_b, g_b)
	double * get_left_bdry() const;
	double * get_right_bdry() const;

	// Evaluate diffusion term k(x) and k'(x)
	double eval_diffusion(vector<double> &x) const;
	double eval_partial_diffusion(vector<double> &x) const;

	// Return whether nonzero reaction is present
	bool reaction_present() const;

	// Evaluate reaction term r(x)
	vector<double> eval_reaction(vector <double> &) const;

	// Return whether nonzero external force is present
	bool extForce_present() const;

	// Evaluate external force f(x)
	double eval_ext_force(vector <double> &x) const;


	// Destructor at the end (memory management)
	~TwoPointBVP();

};