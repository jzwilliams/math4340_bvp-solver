#ifndef PARABOLIC_IBVP_SOLVER_H
#define PARABOLIC_IBVP_SOLVER_H
#include "ParabolicIBVPSolver.h"
#endif

// -------------------------------------------------------------------------------- //
//						Base class for parabolic IBVP								//
// -------------------------------------------------------------------------------- //
ParabolicIBVPSolver::ParabolicIBVPSolver(ApproxBVP * prob) {
	appr_bvp = prob;
	diffusmat = appr_bvp->calcDiffusion();	// Creates memory! See destructor
}

ParabolicIBVPSolver::~ParabolicIBVPSolver() {
	delete diffusmat;
}

// -------------------------------------------------------------------------------- //
//								Forward Euler scheme								//
// -------------------------------------------------------------------------------- //

ForwardEuler::ForwardEuler(ApproxBVP * prob) : ParabolicIBVPSolver(prob) {
	cout << "A forward Euler scheme is utilized." << endl;
}

vector<double> ForwardEuler::one_step_march(vector<double> & timeinterval,
	vector<double> & Ul) {
	
	int k = Ul.size();
	vector<double> Ur(k);

	double tl = timeinterval[0];
	double tr = timeinterval[1];
	double dt = tr - tl;


	// Correction for Dirichlet boundary conditions
	int left = 1; // Tells us whether we include left boundary in internal row
	int right = 0; // Likewise for right boundary

	vector<double> AU = diffusmat->tridiagonal_multiply(Ul);
	
	vector<double> M = appr_bvp->calcLumpedMass();

	vector<double> Ftr = appr_bvp->calcForce(tr);
	if (appr_bvp->left_bound_Dirichlet()) {
		Ur[0] = Ftr[0];
		left = 0;
	}
	if (appr_bvp->right_bound_Dirichlet()) {
		Ur[k - 1] = Ftr[k - 1];
		right = 1;
	}

	// Invert boundary of M if not Dirichlet
	/*vector<double> Minv(M.size());
	for (int i = 1 - left; i < k - right; i++) {
		Minv[i] = 1 / M[i]; // M-values > 0 so okay
	}
	*/

	// Reaction and force terms calculated at tl, Ul
	vector<double> R, prpu;
	appr_bvp->calcReaction(Ul, R, prpu);
	vector<double> Ftl = appr_bvp->calcForce(tl);

	vector<double> Gtl(AU.size());

	// Internal rows
	for (int i = 1 - left; i < k - right; i++) {
		Gtl[i] = AU[i] + R[i] - Ftl[i];
		Gtl[i] *= (dt / M[i]);

		Ur[i] = Ul[i] - Gtl[i];
	}

	return Ur;
}

// -------------------------------------------------------------------------------- //
//								Backward Euler scheme								//
// -------------------------------------------------------------------------------- //

BackwardEuler::BackwardEuler(ApproxBVP * prob) : ParabolicIBVPSolver(prob) {
	cout << "A backward Euler scheme is utilized." << endl;
}


vector<double> BackwardEuler::one_step_march(vector<double> & timeinterval,
	vector<double> & Ul) {
	
	// Attempt 2: retry

	// Newton's parameters
	double tol = 1e-12;
	int maxIter = 20;

	// Basic parameters
	double tl = timeinterval[0];
	double tr = timeinterval[1];
	double dt = tr - tl;
	int k = Ul.size();

	// The eventual output
	vector<double> Ur(k);

	// Parameters that don't change each Newton iteration:
	vector<double> M = appr_bvp->calcLumpedMass();

	vector<double> MUl(k);
	for (int i = 0; i < k; i++) MUl[i] = M[i] * Ul[i];

	tridiagonal_matrix Adt(diffusmat);
	for (int i = 0; i < k - 1; i++) {
		Adt.scalarmult_diagonal_entry(i, dt);
		Adt.scalarmult_lower_diagonal_entry(i, dt);
		Adt.scalarmult_upper_diagonal_entry(i, dt);
	}
	Adt.scalarmult_diagonal_entry(k - 1, dt);

	vector<double> Frdt = appr_bvp->calcForce(tr);
	for (int i = 0; i < k; i++) Frdt[i] *= dt;
	vector<double> Ftr = appr_bvp->calcForce(tr);	// Inefficient, but whatever

	// Parameters that do change each Newton iteration
	vector<double> r, prpu;				// scaled by dt also
	vector<double> AUr(k);				// A*dt*U0
	vector<double> MUr(k);				// M*U0
	vector<double> WUr(k);				// Store -G(tr,Ur)
	vector<double> h(k);				// The Newton step
	double err = -1000;					// The discrepancy (L_infty norm of h)

	// The initial guess
	for (int j = 0; j < k; j++) Ur[j] = Ul[j];

	// The Newton iterator
	for (int i = 0; i <= maxIter; i++) {

		// Dirichlet corrections to Ur
		if (appr_bvp->left_bound_Dirichlet()) Ur[0] = Ftr[0];
		if (appr_bvp->right_bound_Dirichlet()) Ur[k - 1] = Ftr[k - 1];

		appr_bvp->calcReaction(Ur, r, prpu);

		for (int j = 0; j < k; j++) {
			r[j] *= dt;
			prpu[j] *= dt;
			MUr[j] = M[j] * Ur[j];
		}
		tridiagonal_matrix WpUr(Adt);		// Store G'(tr,Ur)

		AUr = Adt.tridiagonal_multiply(Ur);

		// Get overall functions W and W'
		for (int j = 0; j < k; j++) {
			WUr[j] = MUr[j] - MUl[j] + AUr[j] + r[j] - Frdt[j];
			WUr[j] *= -1;
			WpUr.add_to_diagonal_entry(j, (M[j] + prpu[j]));
		}

		// Dirichlet corrections to Jacobian
		if (appr_bvp->left_bound_Dirichlet()) {
			WpUr.set_diagonal_entry(0, 1.0);
			WpUr.set_upper_diagonal_entry(0, 0.0);
		}
		if (appr_bvp->right_bound_Dirichlet()) {
			WpUr.set_diagonal_entry(k - 1, 1.0);
			WpUr.set_lower_diagonal_entry(k - 2, 0.0);
		}

		h = WpUr.solve_linear_system(WUr);

		for (int j = 0; j < k; j++) {
			Ur[j] += h[j];
			if (abs(h[j]) > err) err = h[j];
		}

		if (abs(err) < tol) {
			//cout << "Convergence after "<< i << " iterations" << endl;
			break;
		}
		if (i == maxIter) cout << "No convergence, iterations = " << maxIter << " , err = " << err << endl;
	}

	return Ur;

}

// -------------------------------------------------------------------------------- //
//						Implicit Midpoint Rule scheme								//
// -------------------------------------------------------------------------------- //

ImplicitMidpoint::ImplicitMidpoint(ApproxBVP * prob) : ParabolicIBVPSolver(prob) {
	cout << "An implicit Midpoint Rule scheme is utilized." << endl;
}

vector<double> ImplicitMidpoint::one_step_march(vector<double> & timeinterval,
	vector<double> & Ul) {
	
	vector<double> Ur;

	// IMPLEMENT!!

	return Ur;
}

// -------------------------------------------------------------------------------- //
//						Explicit Midpoint Rule scheme								//
// -------------------------------------------------------------------------------- //
ExplicitMidpoint::ExplicitMidpoint(ApproxBVP * prob) : ParabolicIBVPSolver(prob) {
	cout << "An explicit Midpoint Rule scheme is utilized." << endl;
}

vector<double> ExplicitMidpoint::one_step_march(vector<double> & timeinterval,
	vector<double> & Ul) {

	// Basic parameters
	double tl = timeinterval[0];
	double tr = timeinterval[1];
	double tm = (tr + tl) / 2;
	double dt = tr - tl;
	double dtm = tm - tl;
	int k = Ul.size();

	vector<double> Ur(k);
	vector<double> Um(k);

	// IMPLEMENT!

	return Ur;
}