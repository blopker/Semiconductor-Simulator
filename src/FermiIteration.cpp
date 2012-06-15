/*
 * FermiIteration.cpp
 *
 *  Created on: Jun 7, 2012
 *      Author: Karl Lopker
 */

#include "FermiIteration.h"
using namespace arma;
using namespace std;

double makeINTa(double mu, double ni, double e) {
	return mu * ni * exp(e / (kb * temperature)) * kb * temperature;
}

void interpolate(const vec* in, vec* out) {
	int N = in->n_elem;
	out->at(0) = in->at(0);
	out->at(1) = 0.25 * (2 * in->at(0) + in->at(2) + in->at(1));
	out->at(N - 1) = 0.25 * (2 * in->at(N - 1) + in->at(N - 2) + in->at(N - 3));
	out->at(N) = in->at(N - 1);
	for (int i = 2; i < N - 1; ++i) {
		out->at(i) = 0.25
				* (in->at(i - 1) + in->at(i) + in->at(i + 1) + in->at(i - 2));
	}
}

void fermiIteration(int N, double dx, const vec* ni, const vec* n, const vec* p,
		vec* V, double V_left_B, double V_right_B, double v_bias,
		const vec* mu_p, const vec* mu_n, vec* Efn, vec* Efp) {

	vec INTa_n = zeros<vec>(N + 1, 1);
	vec INTa_p = zeros<vec>(N + 1, 1);

	double e;
	double n_i;
	double mu;

	// initalize first two and last two elements.
	mu = mu_n->at(0);
	e = 0.5 * q * (V_left_B + V->at(0));
	n_i = ni->at(0);
	INTa_n.at(0) = makeINTa(mu, n_i, -e);

	mu = mu_p->at(0);
	INTa_p.at(0) = makeINTa(-mu, n_i, e);

	mu = 0.25 * (2 * mu_n->at(0) + mu_n->at(2) + mu_n->at(1));
	e = 0.25 * q * (V_left_B + V->at(0) + V->at(1) + V->at(2));
	n_i = 0.25 * (2 * ni->at(0) + ni->at(2) + ni->at(1));
	INTa_n.at(1) = makeINTa(mu, n_i, -e);

	mu = 0.25 * (2 * mu_p->at(0) + mu_p->at(2) + mu_p->at(1));
	INTa_p.at(1) = makeINTa(-mu, n_i, e);

	mu = 0.25 * (2 * mu_n->at(N - 1) + mu_n->at(N - 2) + mu_n->at(N - 3));
	e = 0.25 * q * (V_right_B + V->at(N - 1) + V->at(N - 2) + V->at(N - 3));
	n_i = 0.25 * (2 * ni->at(N - 1) + ni->at(N - 2) + ni->at(N - 3));
	INTa_n.at(N - 1) = makeINTa(mu, n_i, -e);

	mu = 0.25 * (2 * mu_p->at(N - 1) + mu_p->at(N - 2) + mu_p->at(N - 3));
	INTa_p.at(N - 1) = makeINTa(-mu, n_i, e);

	mu = mu_n->at(N - 1);
	e = 0.5 * q * (V_right_B + V->at(N - 1));
	n_i = ni->at(N - 1);
	INTa_n.at(N) = makeINTa(mu, n_i, -e);

	mu = mu_p->at(N - 1);
	INTa_p.at(N) = makeINTa(-mu, n_i, e);

//	Generating Interleaved Mesh Coefficients
//	Energy & intrinsic concentration on the Interleaved Mesh
	for (int j = 2; j < (N - 1); ++j) {
		double e;
		double n_i;
		double mu;
		mu = 0.25
				* (mu_n->at(j) + mu_n->at(j - 1) + mu_n->at(j - 2)
						+ mu_n->at(j + 1));
		e = 0.25 * q * (V->at(j) + V->at(j - 1) + V->at(j - 2) + V->at(j + 1));
		n_i = 0.25
				* (ni->at(j) + ni->at(j - 1) + ni->at(j - 2) + ni->at(j + 1));

		INTa_n.at(j) = makeINTa(mu, n_i, -e);

		mu = 0.25
				* (mu_p->at(j) + mu_p->at(j - 1) + mu_p->at(j - 2)
						+ mu_p->at(j + 1));
		INTa_p.at(j) = makeINTa(-mu, n_i, e);
	}

	fermiSolver(&INTa_n, v_bias, true, Efn);

	fermiSolver(&INTa_p, v_bias, false, Efp);

}
