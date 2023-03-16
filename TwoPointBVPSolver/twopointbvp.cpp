#ifndef TWOPOINTBVP_H
#define TWOPOINTBVP_H
#include "twopointbvp.h"
#endif

// Jacob Williams
// MATH 4340
// This file implements all the functions listed in the header!
	// This should be compiled into an object file, not an executable btw

// Constructor for BVP or parabolic IBVP
TwoPointBVP::TwoPointBVP(double *dom, double(*dfunc)(vector<double> &)) {
	// Default: no reaction or external force
	domain = dom;
	diffusion = dfunc;
	containsReaction = false;
	containsExtForce = false;
}

// Constructor for hyperbolic IBVP
TwoPointBVP::TwoPointBVP(double *dom, double(*dfunc)(vector<double> &), double(*pdfunc)(vector<double> &)) {
	// Default: no reaction or external force
	domain = dom;
	diffusion = dfunc;
	partialdiffusion = pdfunc;
	containsReaction = false;
	containsExtForce = false;
}

// The left boundary conditions
void TwoPointBVP::set_left_bdry(bool _leftIsDirichlet, double *val) {
	leftBoundIsDirichlet = _leftIsDirichlet;
	leftBoundVals = val;
}

// The right boundary conditions
void TwoPointBVP::set_right_bdry(bool _rightIsDirichlet, double *val) {
	rightBoundIsDirichlet = _rightIsDirichlet;
	rightBoundVals = val;
}

// Time-dependent boundary conditions
void TwoPointBVP::set_left_bdry(bool _leftIsDirichlet, double gamleft, double(*gleft)(vector<double>&)) {
	leftBoundIsDirichlet = _leftIsDirichlet;
	gammaa = gamleft;
	ga = gleft;
}

void TwoPointBVP::set_right_bdry(bool _rightIsDirichlet, double gamright, double(*gright)(vector<double> &par)) {
	rightBoundIsDirichlet = _rightIsDirichlet;
	gammab = gamright;
	gb = gright;
}

vector<double> TwoPointBVP::calcLeftBdry(double t) const {
	vector<double> par(1);
	par[0] = t;

	vector<double> leftbound(2);
	leftbound[0] = gammaa;
	leftbound[1] = ga(par);
	return leftbound;
}

vector<double> TwoPointBVP::calcRightBdry(double t) const {
	vector<double> par(1);
	par[0] = t;

	vector<double> rightbound(2);
	rightbound[0] = gammaa;
	rightbound[1] = gb(par);

	return rightbound;
}



// The reaction term
void TwoPointBVP::set_reaction(double(*funcone) (vector<double> &),
	double(*functwo) (vector<double> &)) {
	reaction = funcone;
	partialreactionpartialu = functwo;
	containsReaction = true;
}

// The external force
void TwoPointBVP::set_ext_force(double(*func)(vector <double> &)) {
	extForce = func;
	containsExtForce = true;
}

/* 
	Get the properties of the BVP
*/

double * TwoPointBVP::get_domain()  const {
	return domain;
}

bool TwoPointBVP::left_bdry_is_Dirichlet() const {
	return leftBoundIsDirichlet;
}

bool TwoPointBVP::right_bdry_is_Dirichlet() const {
	return rightBoundIsDirichlet;
}

double * TwoPointBVP::get_left_bdry() const {
	return leftBoundVals;
}


double * TwoPointBVP::get_right_bdry() const {
	return rightBoundVals;
}


double TwoPointBVP::eval_diffusion(vector<double> &x) const {
	return diffusion(x);
}

double TwoPointBVP::eval_partial_diffusion(vector<double> &x) const {
	return partialdiffusion(x);
}

bool TwoPointBVP::reaction_present() const {
	return containsReaction;
}

vector<double> TwoPointBVP::eval_reaction(vector<double> & reac) const {
	vector<double> evalreac(2, 0);

	evalreac[0] = reaction(reac); 
	evalreac[1] = partialreactionpartialu(reac);

	return evalreac;
}

bool TwoPointBVP::extForce_present() const {
	return containsExtForce;
}

double TwoPointBVP::eval_ext_force(vector<double> &x) const {
	return extForce(x);
}

TwoPointBVP::~TwoPointBVP() {
	;
}