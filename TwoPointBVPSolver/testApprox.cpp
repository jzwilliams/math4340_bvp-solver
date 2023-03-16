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

#ifndef PARABOLIC_IBVP_SOLVER_H
#define PARABOLIC_IBVP_SOLVER_H
#include "ParabolicIBVPSolver.h"
#endif

#ifndef HYPERBOLIC_IBVP_SOLVER_H
#define HYPERBOLIC_IBVP_SOLVER_H
#include "HyperbolicIBVPSolver.h"
#endif

// Time-testing
#include <chrono>
using std::chrono::high_resolution_clock;
using std::chrono::duration;

using namespace std;

double diffusioncoeff(vector<double> &par);
double partialdiffusioncoeff(vector<double> &par);
double reactcoeff(vector<double> &par);
double partialrcoeff(vector<double> &par);
double forcecoeff(vector<double> &par);
bool true_sol_is_present = true;			// Is there a true solution?

// Time-dependent parameters
double leftbdryvalues(vector<double> &par);
double rightbdryvalues(vector<double> &par);

double initialcondition(double x);
double dinitdt(double x);
double d2initdt2(double x);

double truesol(double x);
double truesol(double x, double t);
double newton(double guess);				// Newton's method, 1-D, for true soln


// -------------------------------------------------------------------------------- //
// -------------------------------------------------------------------------------- //
//							Two-point BVP main() file								//
// -------------------------------------------------------------------------------- //
// -------------------------------------------------------------------------------- //

/* 
int main() {
	// ------------------------------------------------------------ //
	//					Set up the two-point BVP					//
	// ------------------------------------------------------------ //

	double *dom = new double[2];
	dom[0] = 0.0;
	dom[1] = 1.0;

	// Reaction and external force present?
	bool reactPresent = true;
	bool extForcePresent = false;

	// Boundary conditions Dirichlet?
	bool leftDirichlet = true;
	bool rightDirichlet = true;

	// Iteration parameters
	int maxIter = 100;
	double tol = pow(10, -9);


	TwoPointBVP *prob = new TwoPointBVP(dom, diffusioncoeff);

	// gamma_L, g_L
	double *Lbval = new double[2];
	Lbval[0] = 0.0;
	Lbval[1] = 0.0;
	prob->set_left_bdry(leftDirichlet, Lbval);

	double *Rbval = new double[2];
	Rbval[0] = 0.0;
	Rbval[1] = 0.0;
	prob->set_right_bdry(rightDirichlet, Rbval);

	// Set reaction
	if (reactPresent) prob->set_reaction(reactcoeff,partialrcoeff);

	// Set external force
	if (extForcePresent) prob->set_ext_force(forcecoeff);


	// ------------------------------------------------------------ //
	//				Basic information about the BVP					//
	// ------------------------------------------------------------ //

	cout << "Basic information about the two-point BVP: " << endl;

	double *domain = prob->get_domain();
	cout << "Domain: (" << domain[0] << ", " << domain[1] << "). " << endl;

	double *val;
	val = prob->get_left_bdry();
	if (prob->left_bdry_is_Dirichlet()) {
		cout << "Left boundary is Dirichlet. Value: " << val[1] << "." << endl;
	} else {
		cout << "Left boundary not Dirichlet. Values: gamma = " << val[0] << ", g = " << val[1] << "." << endl;
	}

	val = prob->get_right_bdry();
	if (prob->right_bdry_is_Dirichlet()) {
		cout << "Right boundary is Dirichlet. Value: " << val[1] << "." << endl;
	} else {
		cout << "Right boundary not Dirichlet. Values: gamma = " << val[0] << ", g = " << val[1] << "." << endl;
	}

	if (prob->reaction_present()) {
		cout << "Reaction is present." << endl;
	} else {
		cout << "No reaction present." << endl;
	}

	if (prob->extForce_present()) {
		cout << "External force is present." << endl;
	} else {
		cout << "No external force." << endl;
	}

	// ------------------------------------------------------------ //
	//					Approximate the solution					//
	// ------------------------------------------------------------ //

	int numSubInt = 256;
	double *subintervals = new double[numSubInt];

	// Equal-length sub-intervals
	for (int i = 0; i < numSubInt; i++) {
		subintervals[i] = (dom[1] - dom[0]) / numSubInt;
	}

	ApproxBVP *method = new ApproxBVP(numSubInt, subintervals, prob);

	// Put all the info into the approximate solution, and solve it
	vector<double> approxSol = method->solve(maxIter,tol);
	
	vector<double> xcoord = method->getXCoord();

	// Error analysis: L-infty norm
	if (true_sol_is_present) {
		double err = abs(approxSol[0] - truesol(xcoord[0]));
		for (int i = 1; i <= numSubInt; i++) {
			if (abs(approxSol[i] - truesol(xcoord[i])) > err) {
				err = abs(approxSol[i] - truesol(xcoord[i]));
			}
		}
		cout << "Error: " << endl;
		cout << "h" << "\t" << "L-infty norm" << "\t" << endl;
		cout << subintervals[0] << "\t" << err << "\t" << endl;	// pick these up by hand
	}

	ofstream fileout;
	fileout.open("approximatesol.txt");
	for (int i = 0; i <= numSubInt; i++) {
		fileout << xcoord[i] << "\t" << approxSol[i] << " \t" << endl;
	}
	fileout.close();

	if (true_sol_is_present) {
		fileout.open("truesol.txt");
		//int nres = numSubInt;
		//double s = (dom[1] - dom[0]) / nres;
		//double x = dom[0];
		for (int i = 0; i <= numSubInt; i++) {
			fileout << xcoord[i] << "\t" << truesol(xcoord[i]) << "\t" << endl;
			//x += s;
		}
		fileout.close();
	}

	// Plot solution(s)
	ifstream filein;



	// Clean up
	delete method;
	delete[]subintervals;
	delete[]Lbval;
	delete[]Rbval;
	delete prob;
	delete[]dom;

	return 0;
}
*/


// -------------------------------------------------------------------------------- //
// -------------------------------------------------------------------------------- //
//							Parabolic IBVP main() file								//
// -------------------------------------------------------------------------------- //
// -------------------------------------------------------------------------------- //

int main() {
	// Starting time
	auto startTime = high_resolution_clock::now();

	// ------------------------------------------------------------ //
	//						Set up the IBVP							//
	// ------------------------------------------------------------ //

	// Set these manually
	bool truesolpresent = true;
	bool reactPresent = false;
	bool extForcePresent = true;

	double * dom = new double[2];
	dom[0] = -1.0;
	dom[1] = 1.0;

	TwoPointBVP * prob = new TwoPointBVP(dom, diffusioncoeff);

	bool leftBdDirichlet = true;
	bool rightBdDirichlet = true;
	double gammaa = 0;
	double gammab = 0;
	prob->set_left_bdry(leftBdDirichlet, gammaa, leftbdryvalues);
	prob->set_right_bdry(rightBdDirichlet, gammab, rightbdryvalues);

	if (reactPresent) prob->set_reaction(reactcoeff, partialrcoeff);
	if (extForcePresent) prob->set_ext_force(forcecoeff);


	// ------------------------------------------------------------ //
	//					Approximate the solution					//
	// ------------------------------------------------------------ //
	int numSubInt = 8;

	double *subintervals = new double[numSubInt];
	for (int i = 0; i < numSubInt; i++) {
			subintervals[i] = (dom[1] - dom[0]) / numSubInt;
	}
	ApproxBVP * apprbvp = new ApproxBVP(numSubInt, subintervals, prob);

	// The x-coordinates
	vector<double> xcoord = apprbvp->getXCoord();

	// The time-coordinates
	int numTimeSteps = 10000;
	double finalTime = 10;
	double timeStep = static_cast<double>(finalTime) / numTimeSteps;
	vector<double> timeLevel(numTimeSteps + 1);
	timeLevel[0] = 0.0;
	for (int i = 1; i < numTimeSteps + 1; i++) {
		timeLevel[i] = timeLevel[i - 1] + timeStep;
	}

	// Solve the IBVP
	ParabolicIBVPSolver * fdmethod;
	//fdmethod = new ForwardEuler(apprbvp);
	fdmethod = new BackwardEuler(apprbvp);

	vector<double> timeInterval(2);
	vector<double> Ul(numSubInt + 1);
	for (int i = 0; i < numSubInt + 1; i++) {
		Ul[i] = initialcondition(xcoord[i]);
	}

	vector<vector<double>> approxsol(numTimeSteps + 1);
	for (int i = 0; i < numTimeSteps + 1; i++) {
		approxsol[i].resize(numSubInt + 1);
	}

	// First row of approximate solution: initial conditions
	for (int i = 0; i < numSubInt + 1; i++) {
		approxsol[0][i] = Ul[i];
	}

	vector<double> Ur;
	for (int i = 1; i <= numTimeSteps; i++) {
		timeInterval[0] = timeLevel[i - 1];
		timeInterval[1] = timeLevel[i];
		Ur = fdmethod->one_step_march(timeInterval, Ul);
		approxsol[i] = Ur;
		Ul = Ur;
	}

	// Write to file
	ofstream timef, spacef, approxsolf;

	timef.open("timecoord.txt");
	spacef.open("spacecoord.txt");
	approxsolf.open("approximatesol.txt");

	for (int i = 0; i < numTimeSteps + 1; i++) {
		timef << timeLevel[i] << endl;
	}

	for (int i = 0; i < numSubInt + 1; i++) {
		spacef << xcoord[i] << endl;
	}

	for (int i = 0; i < numTimeSteps + 1; i++) {
		for (int j = 0; j < numSubInt + 1; j++) {
			approxsolf << approxsol[i] [j] << " ";
		}
		approxsolf << endl;
	}

	timef.close();
	spacef.close();
	approxsolf.close();

	// The true solution
	if (truesolpresent) {
		vector<vector<double>> truesolvec(numTimeSteps + 1);
		for (int i = 0; i < numTimeSteps + 1; i++) {
			truesolvec[i].resize(numSubInt + 1);
		}
		for (int i = 0; i < numTimeSteps + 1; i++) {
			for (int j = 0; j < numSubInt +1; j++) {
				truesolvec[i] [j] = truesol(xcoord[j], timeLevel[i]);
			}
		}

		ofstream truesolf;
		truesolf.open("truesol.txt");

		for (int i = 0; i < numTimeSteps + 1; i++) {
			for (int j = 0; j < numSubInt + 1; j++) {
				truesolf << truesolvec[i] [j] << " ";
			}
			truesolf << endl;
		}

		truesolf.close();

		/*
		// Ext. force, true
		vector<double> par(2);
		par[1] = timeLevel[numTimeSteps];
		cout << "True force " << endl;
		for (int i = 0; i < numSubInt + 1; i++) {
			par[0] = xcoord[i];
			cout << forcecoeff(par) << " ";
		}
		cout << endl;
		*/

		// Closing time
		auto endTime = high_resolution_clock::now();

		// Time elapsed
		duration<double> elapsedTime = endTime - startTime;

		// Output some basic information about what happened
		cout << "Number of space points: " << numSubInt + 1 << endl;
		cout << "Number of time points: " << numTimeSteps << endl;
		cout << "Time taken: " << elapsedTime.count() << endl;

		// Quantify error at final timestep
		double theErr = -1000;
		for (int i = 0; i < numSubInt + 1; i++) {
			double lifty = abs(approxsol[numTimeSteps][i] - truesolvec[numTimeSteps][i]);
			theErr = (lifty > theErr) ? lifty : theErr;
		}
		cout << "h = " << (dom[1] - dom[0]) / numSubInt << endl;
		cout << "Max. error at final timestep (in L-infty norm) = " << theErr << endl;

		cout << "h \t L_infty" << endl;
		cout << (dom[1] - dom[0]) / numSubInt << "\t" << theErr << endl;
		/*
		cout << "True solution" << endl;
		for (int i = 0; i < numSubInt + 1; i++) {
			cout << truesolvec[numTimeSteps][i] << " ";
		}
		cout << endl;
		cout << "Approximate solution" << endl;
		for (int i = 0; i < numSubInt + 1; i++) {
			cout << approxsol[numTimeSteps][i] << " ";
		}
		cout << endl;
		*/
	}
	


	// ------------------------------------------------------------ //
	//							Clean up							//
	// ------------------------------------------------------------ //
	delete[] dom;
	delete prob;
	delete[] subintervals;
	delete apprbvp;
	delete fdmethod;

	return 0;
}

// -------------------------------------------------------------------------------- //
// -------------------------------------------------------------------------------- //
//							Hyperbolic IBVP main() file								//
// -------------------------------------------------------------------------------- //
// -------------------------------------------------------------------------------- //

double diffusioncoeff(vector<double> &par) {
	return 1;
}

double partialdiffusioncoeff(vector<double> &par) {
	return 0;
}

double reactcoeff(vector<double> &par) {
	return 0;
}

double partialrcoeff(vector<double> &par) {
	return 0;
}

double forcecoeff(vector<double> &par) {
	double x = par[0];
	double t = par[1];
	return exp(-t)*(0.5*(1.0 - t)*(1.0 - x * x) + t);
}

// The true solution, non-time-dependent
double truesol(double x) {
	return 0;
}


// Time-dependent parameters

double leftbdryvalues(vector<double> &par) {
	// Input vector: time
	double t = par[0];
	return exp(-t)*cos(1.0) - 1.0;
}

double rightbdryvalues(vector<double> &par) {
	// Input vector: time
	double t = par[0];
	return exp(-t)*cos(1.0) + 1.0;
}

double initialcondition(double x) {
	return x + cos(x);
}

double dinitdt(double x) {
	return 0;
}

double d2initdt2(double x) {
	return 0;
}

// The true solution, time-dependent
double truesol(double x, double t) {
	return exp(-t)*(cos(x) + 0.5*t*(1.0 - x * x)) + x;


}

double newton(double guess) {
	double theta = guess;
	double lambda = 0.1;		// 0.1, 0.5, 1, 2
	double f = sqrt(2 * lambda)*cosh(0.25*theta) - theta;
	double fp = -0.25*sqrt(2 * lambda)*sinh(0.25*theta) - 1;
	double h = 0;

	for (int i = 0; i < 10000; i++) {
		double temp = f;
		h = -f / fp;
		theta += h;
		f = sqrt(2 * lambda)*cosh(0.25*theta) - theta;
		fp = -0.25*sqrt(2 * lambda)*sinh(0.25*theta) - 1;
		if (abs(temp - f) < 0.00000001) {
			break;
		}
	}

	return theta;
}