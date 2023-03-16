#ifndef APPROX_BVP_H
#define APPROX_BVP_H
#include "approx_bvp.h"
#endif

ApproxBVP::ApproxBVP(int N, const double * subIntLengths,
	const TwoPointBVP * prob) {
	// Constructor: constructs stepLengths (hi matrix); xcoord (x_i); midcoord (x_(i+-1/2))

	numSubInt = N;
	stepLengths = subIntLengths; // length N
	theProblem = prob;
	// Note on boundary conditions: Assigned to leftBoundVals/rightBoundVals [gamma_a, g_a]
		// gamma_a = 0 => Neumann; else Robin
		// (Dirichlet assigned by separate boolean)

	// Calculate Deltax_i
	Deltax.resize(numSubInt + 1); // Deltax: length N+1
	Deltax[0] = 0.5 * stepLengths[0];
	for (int i = 1; i < numSubInt; i++) {
		Deltax[i] = 0.5 * (stepLengths[i] + stepLengths[i - 1]);
	}
	Deltax[numSubInt] = 0.5 * stepLengths[numSubInt - 1];

	// Initialize x_i
	xcoord.resize(numSubInt + 1); // xcoord: length numSubInt
	double * domain = prob->get_domain();
	xcoord[0] = domain[0];
	for (int i = 1; i <= numSubInt; i++) {
		xcoord[i] = xcoord[i - 1] + stepLengths[i - 1];
	}

	// Initialize x_i+-(1/2)
	midcoord.resize(numSubInt + 2); // Length N+2 because we add x0 and xN as endpoints
	midcoord[0] = domain[0];
	for (int i = 1; i <= numSubInt; i++) { // Useful range
		midcoord[i] = 0.5 * (xcoord[i] + xcoord[i - 1]);
	}
	midcoord[numSubInt + 1] = domain[1];
}

bool ApproxBVP::left_bound_Dirichlet(void) {
	return theProblem->left_bdry_is_Dirichlet();
}

bool ApproxBVP::right_bound_Dirichlet(void) {
	return theProblem->right_bdry_is_Dirichlet();
}

tridiagonal_matrix * ApproxBVP::calcDiffusion(void) {
	tridiagonal_matrix * A = new tridiagonal_matrix(numSubInt + 1);
	AssembleDiffusion(A);

	// We delete this in ParabolicIBVPSolver.

	return A;
}

vector<double> ApproxBVP::calcLumpedMass(void) {
	vector<double> D(numSubInt + 1);

	for (int i = 0; i < numSubInt + 1; i++) {
		D[i] = Deltax[i];
	}
	if (theProblem->left_bdry_is_Dirichlet()) D[0] = 0;
	if (theProblem->right_bdry_is_Dirichlet()) D[numSubInt] = 0;

	/*
	if (theProblem->left_bdry_is_Dirichlet()) {
		D[0] = 0;
	} else {
		D[0] = 0.5*stepLengths[0];
	}
	for (int i = 1; i < numSubInt; i++) {
		D[i] = 0.5*stepLengths[i - 1] + 0.5*stepLengths[i];
	}
	if (theProblem->right_bdry_is_Dirichlet()) {
		D[numSubInt] = 0;
	} else {
		D[numSubInt] = 0.5*stepLengths[numSubInt - 1];
	}
	*/

	return D;
}

void ApproxBVP::calcReaction(vector<double>& W, vector<double>& R, vector<double>& prpu) {
	AssembleReaction(W, R, prpu);
}

vector<double> ApproxBVP::calcForce(double timeLevel) {
	vector <double> force(numSubInt + 1);
	vector <double> par(2); // x, t
	par[1] = timeLevel;

	bool forcePresent = theProblem->extForce_present();
	bool leftDirichlet = theProblem->left_bdry_is_Dirichlet();
	bool rightDirichlet = theProblem->right_bdry_is_Dirichlet();

	// Get boundary values
	vector<double> leftBound = theProblem->calcLeftBdry(timeLevel);
	vector<double> rightBound = theProblem->calcRightBdry(timeLevel);

	// Default to zero
	if (!forcePresent) {
		for (int i = 0; i < numSubInt + 1; i++) force[i] = 0.0;
		return force;
	}

	// Left boundary
	if (leftDirichlet) {
		force[0] = leftBound[1];
	} else {
		force[0] = -leftBound[1];
		if (theProblem->extForce_present()) {
			par[0] = xcoord[0];
			force[0] += theProblem->eval_ext_force(par) * Deltax[0];
		}
	}

	// Right boundary
	if (rightDirichlet) {
		force[numSubInt] = rightBound[1];
	} else {
		force[numSubInt] = -rightBound[1];
		if (theProblem->extForce_present()) {
			par[0] = xcoord[numSubInt];
			force[numSubInt] += theProblem->eval_ext_force(par) * Deltax[numSubInt];
		}
	}

	// Internal rows
	if (theProblem->extForce_present()) {
		for (int i = 1; i < numSubInt; i++) {
			par[0] = xcoord[i];
			force[i] = theProblem->eval_ext_force(par) * Deltax[i];
		}
	}

	return force;
}

// BVP accessor
const TwoPointBVP * ApproxBVP::getBVP(void) const {
	return theProblem;
}

vector <double> ApproxBVP::getXCoord() {
	return xcoord;
}


// Assign the diffusion coefficients to the tridiagonal matrix
void ApproxBVP::AssembleDiffusion(tridiagonal_matrix * tmat) {
	// tmat should of dimension (numSubInt + 1). Thus:
		// diag has length numSubInt + 1
		// lowerdiag and upperdiag have length numSubInt

	// Initialize diffusion vector
	vector <double> kappa(numSubInt + 2);
	vector <double> par(1);

	// Get boundary values
	double * leftBound = theProblem->get_left_bdry();
	double * rightBound = theProblem->get_right_bdry();

	// Get diffusion from BVP, including from x_0 and x_N
	for (int i = 0; i < numSubInt + 2; i++) {
		par[0] = midcoord[i];
		kappa[i] = theProblem->eval_diffusion(par);
	}

	// Divide by step length h_i, observing that both kappa(x_0) and kappa(x_(1/2)) are divided by h_1;
		// similarly at right boundary
	kappa[0] /= stepLengths[0];
	kappa[numSubInt + 1] /= stepLengths[numSubInt - 1];
	for (int i = 1; i <= numSubInt; i++) {
		kappa[i] /= stepLengths[i - 1];
	}

	// Fill the inner lines
	for (int i = 1; i < numSubInt; i++) {
		tmat->set_lower_diagonal_entry(i - 1, -kappa[i]);			// Misses last lowerdiag entry
		tmat->set_diagonal_entry(i, kappa[i] + kappa[i + 1]);		// Misses first and last diag entries
		tmat->set_upper_diagonal_entry(i, -kappa[i + 1]);			// Misses first upperdiag entry
	}

	// Assemble left boundary
	if (theProblem->left_bdry_is_Dirichlet()) {						// Dirichlet left boundary
		tmat->set_diagonal_entry(0,1);
		tmat->set_upper_diagonal_entry(0,0);
	} else {														// Neumann left boundary
		tmat->set_diagonal_entry(0, kappa[1]);
		tmat->set_upper_diagonal_entry(0, -kappa[1]);
		tmat->add_to_diagonal_entry(0, -leftBound[0]);
	}

	// Right boundary
	if (theProblem->right_bdry_is_Dirichlet()) {					// Dirichlet right boundary
		tmat->set_diagonal_entry(numSubInt, 1);
		tmat->set_lower_diagonal_entry(numSubInt - 1, 0);
	} else {														// Neumann right boundary
		tmat->set_diagonal_entry(numSubInt, kappa[numSubInt]);
		tmat->set_lower_diagonal_entry(numSubInt - 1, -kappa[numSubInt]);
		tmat->add_to_diagonal_entry(numSubInt, -rightBound[0]);		// adjust for Robin
	}
}

// Build the reaction vector
void ApproxBVP::AssembleReaction(vector<double> & W, vector<double> & R, vector<double> & prpu) {
	// W: contains guess for u
	// R: reaction vector
	// prpu: partial-reaction-partial-u
	// We'll also need to pull xcoord from the ApproxBVP object
	 
	// If no reaction, just set R = prpu = 0
	R.resize(numSubInt + 1);
	prpu.resize(numSubInt + 1);

	if (!theProblem->reaction_present()) {
		for (int i = 0; i < numSubInt + 1; i++) {
			R[i] = 0;
			prpu[i] = 0;
			return;
		}
	}

	// Store xcoord and W in one vector for eval_reaction
	vector<double> input(2, 0.0);

	// Left boundary
	if (theProblem->left_bdry_is_Dirichlet()) {
		R[0] = 0.0;	
		prpu[0] = 0.0;
	} else {
		input[0] = xcoord[0];
		input[1] = W[0];
		vector<double> leftbdr = theProblem->eval_reaction(input);
		R[0] = leftbdr[0] * Deltax[0];
		prpu[0] = leftbdr[1] * Deltax[0];
	}

	// Right boundary
	if (theProblem->right_bdry_is_Dirichlet()) {
		R[numSubInt] = 0.0;
		prpu[numSubInt] = 0.0;
	} else {
		input[0] = xcoord[numSubInt];
		input[1] = W[numSubInt];
		vector<double> rightbdr = theProblem->eval_reaction(input);
		R[numSubInt] = rightbdr[0] * Deltax[numSubInt];
		prpu[numSubInt] = rightbdr[1] * Deltax[numSubInt];
	}

	// Internal rows
	for (int i = 1; i < numSubInt; i++) {
		input[0] = xcoord[i];
		input[1] = W[i];
		vector<double> thereaci = theProblem->eval_reaction(input);
		R[i] = thereaci[0] * Deltax[i];
		prpu[i] = thereaci[1] * Deltax[i];
	}

	return;
			/* Old version

			// Default to zero
			for (int i = 0; i < numSubInt; i++) {
				r[i] = 0.0;
			}

			// Left boundary
			if (theProblem->left_bdry_is_Dirichlet()) {
				r[0] = 0;
			} else {
				r[0] = xcoord[0];
				r[0] = theProblem->eval_reaction() * Deltax[0];
			}

			// Internal rows
			for (int i = 1; i < numSubInt; i++) {
				par[0] = xcoord[i];
				reac[i] = theProblem->eval_reaction(par) * Deltax[i];
			}

			// Right boundary
			if (theProblem->right_bdry_is_Dirichlet()) {
				reac[numSubInt] = 0;
			} else {
				par[0] = xcoord[numSubInt];
				reac[numSubInt] = theProblem->eval_reaction(par) * Deltax[numSubInt];
			}

			return reac;

			*/
}

// Build the external force vector
vector <double> ApproxBVP::AssembleForce() {
	vector <double> force(numSubInt + 1);
	vector <double> par(1);

	// Default to zero
	for (int i = 0; i <= numSubInt; i++) {
		force[i] = 0.0;
	}

	// Get boundary values
	double * leftBound = theProblem->get_left_bdry();
	double * rightBound = theProblem->get_right_bdry();

	// Left boundary
	if (theProblem->left_bdry_is_Dirichlet()) {
		force[0] = leftBound[1];
	} else {
		force[0] = -leftBound[1];
		if (theProblem->extForce_present()) {
			par[0] = xcoord[0];
			force[0] += theProblem->eval_ext_force(par) * Deltax[0];
		}
	}

	// Internal rows
	if (theProblem->extForce_present()) {
		for (int i = 1; i < numSubInt; i++) {
			par[0] = xcoord[i];
			force[i] = theProblem->eval_ext_force(par) * Deltax[i];
		}
	}

	// Right boundary
	if (theProblem->right_bdry_is_Dirichlet()) {
		force[numSubInt] = rightBound[1];
	} else {
		force[numSubInt] = -rightBound[1];
		if (theProblem->extForce_present()) {
			par[0] = xcoord[numSubInt];
			force[numSubInt] += theProblem->eval_ext_force(par) * Deltax[numSubInt];
		}
	}

	return force;
}

// Solve the two-point BVP!!
vector <double> ApproxBVP::solve(int maxIter, double tol) {
	tridiagonal_matrix * Gp, *A;
	vector<double> R, Rp, F;
	A = new tridiagonal_matrix(numSubInt + 1);
	AssembleDiffusion(A);
	F = AssembleForce();

	vector<double> U(numSubInt + 1, 2.5);
	vector<double> h(numSubInt + 1, 0.0);
	vector<double> G(numSubInt + 1, 0.0);
	U[0] = F[0];
	U[numSubInt] = F[numSubInt];


	// The iterative solver method
	for (int iter = 0; iter <= maxIter; iter++) {
		// Discrepancy via L-infty norm
		double discrep(-5);

		// Copy A over to Gp
		Gp = new tridiagonal_matrix(A);
		if (theProblem->reaction_present()) {
			AssembleReaction(U, R, Rp);
			for (int i = 0; i < numSubInt + 1; i++) {
				Gp->add_to_diagonal_entry(i, Rp[i]);
			}
		}

			// Fill G to contain AU + R - F
			vector<double> AU = A->tridiagonal_multiply(U);


			for (int i = 0; i < numSubInt + 1; i++) {
				G[i] = (theProblem->reaction_present()) ? AU[i] + R[i] - F[i] : AU[i] - F[i];
			}

			// Solve linear system for -h
			h = Gp->solve_linear_system(G);

			// Update U based on h
			for (int i = 0; i < numSubInt + 1; i++) {
				U[i] -= h[i];
			}

			// Get maximum of the discrepancies
			for (int i = 0; i < numSubInt + 1; i++) {
				discrep = (abs(h[i]) > discrep) ? abs(h[i]) : discrep;
			}

			// Check for convergence
			if (discrep < tol) {
				cout << "Convergence attained after " << iter << " iterations." << endl;
				break;
			} else if (iter == maxIter) {
				cout << "Convergence not reached after " << maxIter << " iterations." << endl;
				break;
			}

			delete Gp;
	}

	delete A;
	return U;


	/* Old version
	vector <double> sol(numSubInt + 1);

	tridiagonal_matrix * A;
	vector <double> reac, force;

	A = new tridiagonal_matrix(numSubInt + 1);

	// Build tridiagonal matrix from diffusion coefficients
	AssembleDiffusion(A);

	// Assemble reaction and force vectors
	if (theProblem->reaction_present()) {
		reac = AssembleReaction();
		for (int i = 0; i <= numSubInt; i++) {
			A->add_to_diagonal_entry(i, reac[i]);
		}
	}
	
	force = AssembleForce();					// Needed for boundary conditions and initialized at zero anyway

	// Solve the tridiagonal matrix!
	sol = A->solve_linear_system(force);

	return sol;
	*/
}

// Empty destructor.
ApproxBVP::~ApproxBVP() {
	;
}