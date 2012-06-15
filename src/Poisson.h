/*
 * Poisson.h
 *
 *  Created on: Jun 7, 2012
 *      Author: Karl Lopker
 */

#ifndef POISSON_H_
#define POISSON_H_
#include "armadillo"
#include "constants.h"
using namespace arma;

void poisson(int N, double dx, int lap_kind, int itertime, const vec* ni,
		const vec* Efn, const vec* Efp, const vec* epsilon, const vec* Nb,
		double V_left_B, double V_right_B, vec* n, vec* p, vec* V);

#endif /* POISSON_H_ */
