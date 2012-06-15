/*
 * FermiIteration.h
 *
 *  Created on: Jun 7, 2012
 *      Author: Karl Lopker
 */

#ifndef FERMIITERATION_H_
#define FERMIITERATION_H_
#include "armadillo"
#include "constants.h"
#include "FermiSolver.h"
using namespace arma;

void fermiIteration(int N, double dx, const vec* ni, const vec* n,
		const vec* p, vec* V, double V_left_B, double V_right_B,
		double v_bias, const vec* mu_p, const vec* mu_n, vec* Efn,
		vec* Efp);

#endif /* FERMIITERATION_H_ */
