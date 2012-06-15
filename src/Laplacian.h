/*
 * Laplacian.h
 *
 *  Created on: Jun 7, 2012
 *      Author: Karl Lopker
 */

#ifndef LAPLACIAN_H_
#define LAPLACIAN_H_
#include "armadillo"
#include "constants.h"

void laplacian(int N, double dx, int lap_kind, arma::mat* lap);

#endif /* LAPLACIAN_H_ */
