/*
 * FermiSolver.cpp
 *
 *  Created on: Jun 7, 2012
 *      Author: Karl Lopker
 */

#include "FermiSolver.h"
using namespace arma;
using namespace std;

void fermiSolver(const vec* INTa_, const double v_bias, const bool isElectronLevel, vec* Ef_) {
	int N = Ef_->n_elem;

	mat Y = zeros<mat>(N,N);
	vec delta = zeros<vec>(N,1);

	Y.diag(1) = INTa_->rows(1,N-1);
	Y.diag(-1) = INTa_->rows(1,N-1);
	Y.diag() = -(INTa_->rows(0,N-1) + INTa_->rows(1,N));

	double p = -exp((0.5*q*v_bias)/(kb*temperature));
	double n = -exp((-0.5*q*v_bias)/(kb*temperature));

	if(!isElectronLevel){
		double t = n;
		n = p;
		p = t;
	}

	delta.at(0) = INTa_->at(0)*n;
	delta.at(N-1) = INTa_->at(N)*p;

	double d = (isElectronLevel)?1:-1;
	*Ef_ = solve(Y,delta);
	*Ef_ = d*kb*temperature*log(*Ef_);
}

