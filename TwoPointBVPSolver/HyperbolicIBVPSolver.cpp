#ifndef HYPERBOLIC_IBVP_SOLVER_H
#define HYPERBOLIC_IBVP_SOLVER_H
#include "HyperbolicIBVPSolver.h"
#endif

HyperbolicIBVPSolver::HyperbolicIBVPSolver(ApproxBVP * prob) {
	appr_bvp = prob;
	diffusmat = appr_bvp->calcDiffusion();
	theProblem = appr_bvp->getBVP();
}

HyperbolicIBVPSolver::~HyperbolicIBVPSolver() {
	delete diffusmat;
}

TwoStepImplicitRPR::TwoStepImplicitRPR(ApproxBVP * prob) : HyperbolicIBVPSolver(prob) {
	cout << "A two-step implicit Right Point Rule scheme is utilized." << endl;
}

vector<double> TwoStepImplicitRPR::two_step_march(vector<double> &timelevel,
	vector<double> &Unm1, vector<double> &Un) {
	
	// Basic parameters
	int k = Un.size();
	double tol = 1e-9;
	double maxIter = 20;

	double tnm1 = timelevel[0];
	double tn = timelevel[1];
	double tnp1 = timelevel[2];
	double tnphalf = (tn + tnp1) / 2;
	double dtnp1 = tnp1 - tn;
	double dtn = tn - tnm1;
	double dtav = (tnp1 - tnm1) / 2;

	// In case of Dirichlet conditions
	bool leftDirichlet = appr_bvp->left_bound_Dirichlet();
	bool rightDirichlet = appr_bvp->right_bound_Dirichlet();
	vector<double> Ftnp1 = appr_bvp->calcForce(tnp1);

	// Parameters that don't change each Newton iteration
	vector<double> M = appr_bvp->calcLumpedMass();					// Lumped mass
	vector<double> Mnp1(k);										// Divided by tnp1
	vector<double> Mnp1Un(k);										// Lumped mass * Un / tnp1
	vector<double> MnUn(k);											// Lumped mass * Un / tn
	vector<double> MnUnm1(k);										// Lumped mass * Unm1 / tn
	vector<double>Ftphalf = appr_bvp->calcForce(tnphalf);			// Force(tnphalf) * dtav
	for (int i = 0; i < k; i++) {
		Mnp1[i] = M[i] / dtnp1;
		Mnp1Un[i] = M[i] * Un[i] / tnp1;
		MnUn[i] = M[i] * Un[i] / tn;
		MnUnm1[i] = M[i] * Unm1[i] / tn;
		Ftphalf[i] *= dtav;
	}

	tridiagonal_matrix Atav(diffusmat);								// Diffusion * dtav
	for (int i = 0; i < k - 1; i++) {
		Atav.scalarmult_diagonal_entry(i, dtav);
		Atav.scalarmult_upper_diagonal_entry(i, dtav);
		Atav.scalarmult_lower_diagonal_entry(i, dtav);
	}
	Atav.scalarmult_diagonal_entry(k - 1, dtav);

	// Initial guess for output
	vector<double> Unp1(k);
	for (int i = 0; i < k; i++) {
		Unp1[i] = Un[i];
	}

	// Parameters that do change each Newton iteration
	vector<double> r, prpu;											// Reaction and partialreactionpartialu
	vector<double> Uav(k);											// (Un + Unp1) / 2
	vector<double> AtavUav(k);										// Atav * Uav
	vector<double> Mnp1Unp1(k);										// Mtnp1 * Unp1
	vector<double> W(k);											// The overall Newton function
	tridiagonal_matrix Wp(Atav);									// Its derivative
	vector<double> h(k);											// The Newton step
	double err = -1000;												// L-infty norm of h

	// The Newton iterator
	for (int i = 1; i <= maxIter; i++) {
		// Dirichlet correction for Unp1
		if (leftDirichlet) Unp1[0] = Ftnp1[0];
		if (rightDirichlet) Unp1[k] = Ftnp1[k];

		for (int j = 0; j < k; j++) {
			Uav[j] = (Un[j] + Unp1[j]) / 2;
			Mnp1Unp1[j] = Mnp1[j] * Unp1[j];
		}

		AtavUav = Atav.tridiagonal_multiply(Uav);

		appr_bvp->calcReaction(Uav, r, prpu);
		for (int j = 0; j < k; j++) {
			r[j] *= dtav;
			prpu[j] *= dtav;
		}

		// Set up the Newton functions to solve
		for (int j = 0; j < k; j++) {
			W[j] = Mnp1Unp1[j] - Mnp1Un[j] - MnUn[j] + MnUnm1[j] + AtavUav[j] + r[j] - Ftphalf[j];
			W[j] *= -1;
			Wp.add_to_diagonal_entry(j, (Mnp1[j] + prpu[j]));
		}

		// Dirichlet correction to Wp
		if (leftDirichlet) {
			Wp.set_diagonal_entry(0, 1.0);
			Wp.set_upper_diagonal_entry(0, 0.0);
		}
		if (rightDirichlet) {
			Wp.set_diagonal_entry(k - 1, 1.0);
			Wp.set_lower_diagonal_entry(k - 1, 0.0);
		}

		h = Wp.solve_linear_system(W);

		for (int j = 0; j < k; j++) {
			Unp1[j] += h[j];
			if (abs(h[j]) > err) err = h[j];
		}

		if (abs(err) < tol) {
			cout << "Convergence after " << i << "iterations" << endl;
			break;
		}
		if (i == maxIter) cout << "No convergence, " << i << " iterations, err = " << err << endl;

	}

	return Unp1;
}

vector<double> TwoStepImplicitRPR::first_march(vector<double> &timelevel,
	vector<double> &U0, vector<double> & didt, vector<double> & d2idt2) {

	// Basic parameters
	int k = U0.size();

	double t1 = timelevel[1];
	double t0 = timelevel[0];
	double dt = t1 - t0;

	// The eventual return value
	vector<double> U1(k);

	// Dirichlet parameters
	bool leftDirichlet = appr_bvp->left_bound_Dirichlet();
	bool rightDirichlet = appr_bvp->right_bound_Dirichlet();
	int left = 1;
	int right = 0;
	vector<double> F1 = appr_bvp->calcForce(t1);
	if (leftDirichlet) {
		U1[0] = F1[0];
		left = 0;
	}
	if (rightDirichlet) {
		U1[k - 1] = F1[k - 1];
		right = 1;
	}

	// The initial conditions
	vector<double> Q0 = U0;
	vector<double> Q1 = didt;
	vector<double> Q2 = d2idt2;

	// The initial force
	vector<double> F0 = appr_bvp->calcForce(t0);

	// The reaction
	vector<double> r, prpu;
	appr_bvp->calcReaction(Q0, r, prpu);

	// The diffusion
	vector<double> kappa(k);
	vector<double> pkappa(k);
	vector<double> par(1);
	for (int i = 0; i < k; i++) {
		par[0] = Q0[i];
		kappa[i] = theProblem->eval_diffusion(par);
		pkappa[i] = theProblem->eval_partial_diffusion(par);
	}

	// The terms for our overall equation
	vector<double> Q1dt(k);
	vector<double> Fdt(k);
	vector<double> Rdt(k);
	vector<double> ktQ2(k);
	vector<double> kptQ1(k);
	for (int i = 0; i < k; i++) {
		Q1dt[i] = Q1[i] * dt;
		Fdt[i] = F0[i] * dt * dt / 2.0;
		Rdt[i] = r[i] * dt * dt / 2.0;
		ktQ2[i] = kappa[i] * Q2[i] * dt * dt / 2.0;
		kptQ1[i] = pkappa[i] * Q1[i] * dt * dt / 2.0;
	}

	// The overall equation
	for (int i = 1 - left; i < (k - 1) - right; i++) {
		U1[i] = Q0[i] + Q1dt[i] + Fdt[i] - Rdt[i] + ktQ2[i] + kptQ1[i];
	}

	return U1;
}