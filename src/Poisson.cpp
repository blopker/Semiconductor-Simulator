/*
 * Poisson.cpp
 *
 *  Created on: Jun 7, 2012
 *      Author: Karl Lopker
 */
#include "Poisson.h"
#include "Laplacian.h"
using namespace arma;
using namespace std;

void poisson(int N, double dx, int lap_kind, int itertime, const vec* ni,
		const vec* Efn, const vec* Efp, const vec* epsilon, const vec* Nb,
		double V_left_B, double V_right_B, vec* n, vec* p, vec* V) {
	vec dV = zeros<vec>(N, 1);
	mat lap = zeros<mat>(N, N);
	mat A = zeros<mat>(N, N);
	laplacian(N, dx, lap_kind, &lap);

	for (int i = 0; i < itertime; ++i) {
		double max = (dV.max() > abs(dV.min())) ? dV.max() : abs(dV.min());
		if (max < 1) {
			*V += dV;
		} else {
			*V += dV / max / 20;
		}

		*n = *ni % exp(*Efn / (kb * temperature))
				% exp(*V * -q / (kb * temperature));
		*p = *ni % exp(*Efp * -1 / (kb * temperature))
				% exp(*V * q / (kb * temperature));

		vec b = -(lap * *V) + q / *epsilon % (-*n + *p - *Nb);
		b.at(0) -= lap.at(1, 0) * (V_left_B);
		b.at(N - 1) -= lap(1, 2) * (V_right_B);

		A = lap;
		A.diag() += pow(q, 2) / (*epsilon * kb * temperature) % (-*n - *p);

		dV = solve(A, b);
		max = (b.max() > abs(b.min())) ? b.max() : abs(b.min());

		if (max < 100) {
			break;
		}
	}
}

