/*
 * FermiSolver.h
 *
 *  Created on: Jun 7, 2012
 *      Author: Karl Lopker
 */

#ifndef FERMISOLVER_H_
#define FERMISOLVER_H_
#include "armadillo"
#include "constants.h"
using namespace arma;

void fermiSolver(const vec* INTa_, const double v_bias, const bool isElectronLevel, vec* Ef_);

#endif /* FERMISOLVER_H_ */
